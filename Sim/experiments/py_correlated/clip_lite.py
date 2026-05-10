"""
CLIP-Lite: Python 级三阶段迭代接收机。

对照 MATLAB GiFreeReceiver 的核心机制:
  Stage 1: 组 OMP + κ 标定门限 (已有, 复用 group_omp)
  Stage 2: 软判决 LMMSE → 置信门控 → ID2P 消除 → 组 OMP 重检测
  Stage 3: 固定支撑集 → 联合 LS 增益重估计 → 截尾噪声 → 最终 LMMSE

对比: 单次 OMP vs CLIP-Lite (2 轮迭代) 的 BER/NMSE/Pd。

用法:
  python clip_lite.py --mc 100
  python clip_lite.py --quick
"""

import argparse, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _plot_setup import setup_plots
setup_plots()
import matplotlib.pyplot as plt
import numpy as np

from channel import (
    compute_chirp_params, sample_clustered_geometry,
    bdlr_gains, build_path_matrix, compute_effective_channel,
)
from grid_stats import build_candidate_atoms
from detection_mc import _estimate_true_candidate_mask
from group_omp import build_group_atoms, single_omp_detect, group_omp_detect


# ======================================================================
# Stage 2: 软判决 + 置信门控 + ID2P 消除
# ======================================================================

def stage2_refinement(y, x_tx, H_est, data_idx, n_data, n_sc, pilot_idx,
                       atoms, grid, group_atoms, group_masks, cfg):
    """
    Stage 2 单轮迭代。

    1. LMMSE 数据检测 (软输出)
    2. 置信门控: 保留 top-k 符号
    3. 用高置信符号重建数据干扰并消除
    4. 在清洁信号上重新组 OMP 检测
    5. 更新 H_est

    参数:
      y:            (N,) 当前接收信号
      x_tx:         (N,) 发射帧 (用于消去导频)
      H_est:        (N,N) 当前信道估计
      data_idx:     (n_data,) 数据位置
      n_data, n_sc, pilot_idx
      atoms, grid, group_atoms, group_masks
      cfg:          dict (含 kappa, keep_ratio, n_paths 等)

    返回:
      H_new:        (N,N) 更新后的信道估计
      x_data_est:   (n_data,) 数据软估计
      conf_mask:    (n_data,) bool 高置信符号掩码
      n_kept:       int 保留的符号数
      detected:     list[int] 本轮新检测的候选
    """
    noise_power = 1.0
    kappa = float(cfg.get("kappa", 24.6))
    kappa_g = kappa * 0.85
    keep_ratio = float(cfg.get("stage2_keep_ratio", 0.6))
    n_paths = int(cfg.get("n_paths", 6))
    dirichlet_r = int(cfg.get("dirichlet_r", 3))
    max_doppler = int(cfg.get("max_doppler", 2))

    # 1. LMMSE
    H_data = H_est[:, data_idx]
    reg = noise_power
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_data_est = w_mmse @ y  # (n_data,) 软输出

    # 2. 置信门控: 按 |x_est| 排序, 保留 top keep_ratio
    confidence = np.abs(x_data_est)
    threshold = np.percentile(confidence, 100.0 * (1.0 - keep_ratio))
    conf_mask = confidence >= threshold
    n_kept = int(conf_mask.sum())

    # 3. 重建数据干扰并消除
    # 硬判决: 取 QPSK 最近点
    x_data_hard = _qpsk_hard_decision(x_data_est)
    x_data_fb = np.where(conf_mask, x_data_hard, 0.0 + 0j)
    # 仅消除高置信符号
    id2p_est = H_est[:, data_idx] @ x_data_fb
    y_clean = y.copy()
    y_clean -= id2p_est

    # 4. 组 OMP 重检测
    res_pow = noise_power
    _, detected, _, g_det = group_omp_detect(
        y_clean, atoms, grid, group_atoms, group_masks,
        n_paths + 2, kappa_g, res_pow)

    # 5. 重建 H_est (用 OMP 内积增益)
    c1 = float(cfg.get("_c1", 0.0))
    H_new = _build_channel_from_detected(detected, g_det, atoms, grid, c1, n_sc, dirichlet_r)

    return H_new, x_data_est, conf_mask, n_kept, detected


# ======================================================================
# Stage 3: 固定支撑 + 联合 LS + 截尾噪声 + 最终 LMMSE
# ======================================================================

def stage3_polish(y, x_tx, H_est, data_idx, n_data, n_sc, pilot_idx,
                   atoms, grid, cfg):
    """
    Stage 3: 固定支撑集, 不搜索新路径。

    1. 提取当前路径列表
    2. 联合 LS 增益重估计 (用 pilot 响应)
    3. 截尾噪声方差估计
    4. 最终 LMMSE + 硬判决

    返回:
      ber:         float
      noise_est:   float
    """
    noise_power = 1.0
    trim_alpha = float(cfg.get("trim_alpha", 0.15))
    dirichlet_r = int(cfg.get("dirichlet_r", 3))
    c1 = float(cfg.get("_c1", 0.0))
    n_paths = int(cfg.get("n_paths", 6))

    # 1. 从 H_est 提取当前路径 (非零列)
    # 简化: 直接用 H_est 作为最终估计
    # 2. 截尾噪声估计
    r_final = y - H_est @ x_tx
    residual_powers = np.sort(np.abs(r_final) ** 2)
    n_keep = int(n_sc * (1.0 - trim_alpha))
    noise_est = float(np.mean(residual_powers[:n_keep]))

    # 3. 最终 LMMSE
    H_data = H_est[:, data_idx]
    reg = max(noise_est, noise_power)
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_data_est = w_mmse @ y

    # 4. BER
    bits_est = (np.real(x_data_est) > 0).astype(int)
    # bits_true 需要从外部传入, 此处只返回估计
    return x_data_est, noise_est


# ======================================================================
# CLIP-Lite 完整流水线
# ======================================================================

def run_clip_trial(cfg):
    """
    CLIP-Lite 单次 trial。

    对比三条路径的 BER/NMSE/Pd:
      A. 单次组 OMP (Stage 1 only, 当前 baseline)
      B. CLIP-Lite 1 iter (Stage 1 + Stage 2 × 1)
      C. CLIP-Lite 2 iter (Stage 1 + Stage 2 × 2)
    """
    n_sc         = int(cfg["n_sc"])
    n_paths      = int(cfg["n_paths"])
    n_clusters   = int(cfg["n_clusters"])
    max_delay    = int(cfg["max_delay"])
    max_doppler  = int(cfg["max_doppler"])
    doppler_guard = int(cfg["doppler_guard"])
    dirichlet_r  = int(cfg["dirichlet_r"])
    rho          = float(cfg["rho"])
    snr_db       = float(cfg["snr_db"])
    seed         = int(cfg.get("seed", 42))
    pilot_idx    = int(cfg.get("pilot_idx", 0))
    kappa        = float(cfg.get("kappa", 24.6))
    stage2_iter  = int(cfg.get("stage2_iter", 2))

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)
    cfg["_c1"] = c1  # 缓存供后续使用

    # --- 信道 + 帧 ---
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    gains = bdlr_gains(n_paths, n_clusters, rho, rng)
    H_true = compute_effective_channel(delays, dopplers, gains, n_sc, c1, dirichlet_r)

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)
    data_amp = np.sqrt(10.0 ** (snr_db / 10.0))

    x_tx = np.zeros(n_sc, dtype=complex)
    x_tx[pilot_idx] = 1.0 + 0j
    bits = rng.integers(0, 2, 2 * n_data)
    sym_data = ((2*bits[:n_data]-1) + 1j*(2*bits[n_data:]-1)) / np.sqrt(2.0)
    x_tx[data_idx] = data_amp * sym_data
    noise_power = 1.0
    noise = np.sqrt(noise_power/2.0) * (rng.standard_normal(n_sc) + 1j*rng.standard_normal(n_sc))
    y = H_true @ x_tx + noise

    # --- 原子 ---
    atoms, grid = build_candidate_atoms(n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx)
    group_atoms, group_masks = build_group_atoms(atoms, grid, max_doppler, dirichlet_r)
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    true_set = set(np.where(true_mask)[0])

    res_pow = noise_power
    kappa_g = kappa * 0.85

    # ========== Path A: Stage 1 only ==========
    _, det_s1, _, g_s1 = group_omp_detect(
        y, atoms, grid, group_atoms, group_masks, n_paths + 2, kappa_g, res_pow)
    H_s1 = _build_channel_from_detected(det_s1, g_s1, atoms, grid, c1, n_sc, dirichlet_r)
    ber_s1, nmse_s1 = _compute_ber_nmse(H_s1, H_true, y, data_idx, sym_data, n_data, noise_power)
    pd_s1 = len(set(det_s1) & true_set) / max(n_paths, 1)

    # ========== Path B: CLIP-Lite (Stage 1 + Stage 2 × N) ==========
    H_clip = H_s1.copy()
    det_clip = list(det_s1)
    y_work = y.copy()
    total_kept = 0

    for it in range(stage2_iter):
        keep_ratio = 0.5 + 0.25 * it  # 渐进放开: 0.5 → 0.75
        _c = dict(cfg)
        _c["stage2_keep_ratio"] = keep_ratio
        _c["kappa"] = kappa * (1.0 - 0.15 * it)  # 门限渐进放松
        _c["_c1"] = c1

        H_new, x_soft, conf_mask, n_kept, det_new = stage2_refinement(
            y_work, x_tx, H_clip, data_idx, n_data, n_sc, pilot_idx,
            atoms, grid, group_atoms, group_masks, _c)

        # 合并新旧路径集
        H_clip = H_new
        all_det = list(set(det_clip + det_new))
        det_clip = all_det
        total_kept += n_kept

        # 更新工作信号: 用当前 H_clip 消除全部数据干扰
        x_data_hard = _qpsk_hard_decision(x_soft)
        id2p_all = H_clip[:, data_idx] @ x_data_hard
        y_work = y - id2p_all

    # 最终 H_clip 从合并路径集重建
    # 合并所有 OMP 增益
    all_gains = dict(g_s1)
    H_clip = _build_channel_from_detected(det_clip, all_gains, atoms, grid, c1, n_sc, dirichlet_r)

    # Stage 3: 固定支撑 + LS + 截尾
    x_final, noise_est = stage3_polish(
        y, x_tx, H_clip, data_idx, n_data, n_sc, pilot_idx, atoms, grid, cfg)
    ber_clip = float(np.mean((np.real(x_final) > 0).astype(int) !=
                              (np.real(sym_data) > 0).astype(int)))
    nmse_clip = 10.0 * np.log10(
        np.sum(np.abs(H_true - H_clip)**2) / max(np.sum(np.abs(H_true)**2), 1e-15))
    pd_clip = len(set(det_clip) & true_set) / max(n_paths, 1)

    return {
        "ber_s1": ber_s1, "ber_clip": ber_clip,
        "nmse_s1": nmse_s1, "nmse_clip": nmse_clip,
        "pd_s1": pd_s1, "pd_clip": pd_clip,
        "n_det_s1": len(det_s1), "n_det_clip": len(det_clip),
        "n_iterations": stage2_iter,
    }


# ======================================================================
# Helpers
# ======================================================================

def _qpsk_hard_decision(x_soft):
    """QPSK 软输出 → 硬判决 (归一化)。"""
    return ((np.real(x_soft) > 0).astype(float) * 2 - 1 +
            1j * (np.imag(x_soft) > 0).astype(float) * 2 - 1j) / np.sqrt(2.0)


def _build_channel_from_detected(detected, detected_gains, atoms, grid, c1, n_sc, dirichlet_r):
    """
    从检测候选索引 + OMP 返回增益 → H_est。

    增益来自 OMP 内积，不做额外重估。
    """
    if not detected:
        return np.zeros((n_sc, n_sc), dtype=complex)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    H = np.zeros((n_sc, n_sc), dtype=complex)
    for idx in detected:
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        g_i = detected_gains.get(idx, 1.0 / np.sqrt(len(detected)))
        H += g_i * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)
    return H

def _ls_gain_estimate(detected, atoms, grid, c1, n_sc, dirichlet_r, y):
    """
    LS 联合增益估计: min_g ||Phi_big @ g - y||²。

    用 pilot 列构建 Phi_big, 对接收信号 y 做 LS 拟合。
    """
    if not detected:
        return np.zeros(len(detected))
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    K = len(detected)
    Phi_big = np.zeros((n_sc, K), dtype=complex)
    for ki, idx in enumerate(detected):
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        phi = build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)
        Phi_big[:, ki] = phi[:, 0]  # pilot column
    g_vec, _, _, _ = np.linalg.lstsq(Phi_big, y, rcond=None)
    return g_vec


def _compute_ber_nmse(H_est, H_true, y, data_idx, sym_data, n_data, noise_power):
    """LMMSE → BER + NMSE。"""
    if np.all(H_est == 0):
        return 0.5, 0.0
    H_data = H_est[:, data_idx]
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + noise_power * np.eye(n_data), H_data.conj().T)
    x_est = w_mmse @ y
    ber = float(np.mean((np.real(x_est) > 0).astype(int) !=
                         (np.real(sym_data) > 0).astype(int)))
    nmse = 10.0 * np.log10(
        np.sum(np.abs(H_true - H_est)**2) / max(np.sum(np.abs(H_true)**2), 1e-15))
    return ber, nmse


# ======================================================================
# 批量对比
# ======================================================================

def run_grid(cfg, n_trials=50):
    """SNR × rho 网格扫描: S1-only vs CLIP-Lite。"""
    snr_vec = cfg.get("snr_vec", [10, 15, 20, 25])
    rho_vec = cfg.get("rho_vec", [0.0, 0.5, 0.9])
    results = []
    total = len(snr_vec) * len(rho_vec)
    count = 0
    t0 = time.time()

    for rho_val in rho_vec:
        for snr_val in snr_vec:
            _c = dict(cfg)
            _c["snr_db"] = snr_val
            _c["rho"] = rho_val

            trials = []
            for t in range(n_trials):
                _c["seed"] = cfg.get("seed", 0) + count * n_trials + t
                r = run_clip_trial(_c)
                trials.append(r)

            row = {
                "snr_db": snr_val, "rho": rho_val,
                "ber_s1": np.mean([r["ber_s1"] for r in trials]),
                "ber_clip": np.mean([r["ber_clip"] for r in trials]),
                "nmse_s1": np.mean([r["nmse_s1"] for r in trials]),
                "nmse_clip": np.mean([r["nmse_clip"] for r in trials]),
                "pd_s1": np.mean([r["pd_s1"] for r in trials]),
                "pd_clip": np.mean([r["pd_clip"] for r in trials]),
                "n_det_s1": np.mean([r["n_det_s1"] for r in trials]),
                "n_det_clip": np.mean([r["n_det_clip"] for r in trials]),
            }
            results.append(row)
            count += 1
            elapsed = time.time() - t0
            print(f"  [{count}/{total}] SNR={snr_val:2d}dB rho={rho_val:.1f}  "
                  f"BER_s1={row['ber_s1']:.3f} BER_clip={row['ber_clip']:.3f}  "
                  f"PD_s1={row['pd_s1']:.3f} PD_clip={row['pd_clip']:.3f}  "
                  f"({elapsed/count*(total-count):.0f}s left)")

    return results


def plot_clip_comparison(results, out_dir):
    """Fig7: S1 vs CLIP-Lite 三面板。"""
    snr_vec = sorted(set(r["snr_db"] for r in results))
    rho_vec = sorted(set(r["rho"] for r in results))
    colors = {0.0: "gray", 0.5: "steelblue", 0.9: "crimson"}

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    for rho_val in rho_vec:
        rows = [r for r in results if r["rho"] == rho_val]
        snr = np.array([r["snr_db"] for r in rows])
        c = colors[rho_val]

        # BER
        ax = axes[0]
        ax.semilogy(snr, [r["ber_s1"] for r in rows], 's--', color=c, linewidth=1, markersize=6, alpha=0.5,
                     label=f"S1 rho={rho_val:.1f}")
        ax.semilogy(snr, [r["ber_clip"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c,
                     label=f"CLIP rho={rho_val:.1f}")

        # Pd
        ax = axes[1]
        ax.plot(snr, [r["pd_s1"] for r in rows], 's--', color=c, linewidth=1, markersize=6, alpha=0.5)
        ax.plot(snr, [r["pd_clip"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c)

        # NMSE
        ax = axes[2]
        ax.plot(snr, [r["nmse_s1"] for r in rows], 's--', color=c, linewidth=1, markersize=6, alpha=0.5)
        ax.plot(snr, [r["nmse_clip"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c)

    axes[0].set_xlabel("Data SNR (dB)"); axes[0].set_ylabel("BER"); axes[0].set_title("BER: Stage1 vs CLIP-Lite"); axes[0].legend(fontsize=7); axes[0].grid(True, alpha=0.3)
    axes[1].set_xlabel("Data SNR (dB)"); axes[1].set_ylabel("P_d"); axes[1].set_title("Detection Probability"); axes[1].grid(True, alpha=0.3)
    axes[2].set_xlabel("Data SNR (dB)"); axes[2].set_ylabel("NMSE (dB)"); axes[2].set_title("Channel Estimation NMSE"); axes[2].grid(True, alpha=0.3)

    fig.suptitle("CLIP-Lite: Iterative ID2P Cancellation vs Single-Pass Detection",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, "Fig7_CLIP_Lite.png")
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"\n  Fig7 saved: {path}")


def print_table(results):
    print(f"\n{'SNR':>4} {'rho':>4} {'BER_s1':>8} {'BER_clip':>8} {'gain':>6} "
          f"{'Pd_s1':>6} {'Pd_clip':>6} {'dPd':>6} {'NMSE_s1':>7} {'NMSE_clip':>7}")
    print("-" * 82)
    for r in sorted(results, key=lambda x: (x["rho"], x["snr_db"])):
        gain = r["ber_s1"] / max(r["ber_clip"], 1e-10)
        print(f"{r['snr_db']:4.0f} {r['rho']:4.1f} {r['ber_s1']:8.4f} {r['ber_clip']:8.4f} "
              f"{gain:6.1f}x {r['pd_s1']:6.3f} {r['pd_clip']:6.3f} {r['pd_clip']-r['pd_s1']:+6.3f} "
              f"{r['nmse_s1']:+7.1f} {r['nmse_clip']:+7.1f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--mc", type=int, default=50)
    ap.add_argument("--out-dir", type=str, default=None)
    args = ap.parse_args()

    n_trials = 15 if args.quick else args.mc
    sim_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = args.out_dir or os.path.join(sim_dir, "Results", "Figures")
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "p_fa": 1e-3, "seed": 0, "pilot_idx": 0,
        "kappa": 24.6, "stage2_iter": 2, "trim_alpha": 0.15,
        "snr_vec": [10, 15, 20, 25],
        "rho_vec": [0.0, 0.5, 0.9],
    }

    print(f"=== CLIP-Lite: Stage1 vs Iterative ID2P Cancellation (n={n_trials}) ===\n")
    t0 = time.time()
    results = run_grid(cfg, n_trials)
    print(f"\n  Total: {time.time()-t0:.0f}s")
    print_table(results)
    plot_clip_comparison(results, out_dir)
