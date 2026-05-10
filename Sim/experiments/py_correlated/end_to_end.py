#!/usr/bin/env python3
"""
全链路仿真：BDLR 信道 → 组/单 OMP → LMMSE → BER/NMSE/Pd。

对比维度:
  - 检测器: 单径 OMP vs 组级 OMP
  - 信道:    i.i.d. (ρ=0) vs BDLR 相关 (ρ=0.5, 0.9)
  - SNR:     [0, 5, 10, 15, 20, 25] dB

用法:
  python end_to_end.py                  # 默认参数, 100 trials/point
  python end_to_end.py --quick          # 小规模快速验证 (20 trials)
  python end_to_end.py --mc 200         # 200 trials/point
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
from group_omp import (
    build_group_atoms,
    single_omp_detect,
    group_omp_detect,
)


# ======================================================================
# 单次 trial
# ======================================================================

def run_one_trial(cfg):
    """
    一次完整全链路 trial。

    返回 dict: ber_s, ber_g, nmse_s, nmse_g, pd_s, pd_g, n_det_s, n_det_g
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
    p_fa         = float(cfg.get("p_fa", 1e-3))
    seed         = int(cfg.get("seed", 42))
    pilot_idx    = int(cfg.get("pilot_idx", 0))
    kappa        = float(cfg.get("kappa", 24.6))

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)

    # --- 信道 ---
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    gains = bdlr_gains(n_paths, n_clusters, rho, rng)
    H_true = compute_effective_channel(delays, dopplers, gains, n_sc, c1, dirichlet_r)

    # --- 帧 ---
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)
    data_amp = np.sqrt(10.0 ** (snr_db / 10.0))

    x_tx = np.zeros(n_sc, dtype=complex)
    x_tx[pilot_idx] = 1.0 + 0j
    bits = rng.integers(0, 2, 2 * n_data)
    sym_data = ((2*bits[:n_data]-1) + 1j*(2*bits[n_data:]-1)) / np.sqrt(2.0)
    x_tx[data_idx] = data_amp * sym_data
    noise_power = 1.0
    noise = np.sqrt(noise_power/2.0)*(rng.standard_normal(n_sc) + 1j*rng.standard_normal(n_sc))
    y = H_true @ x_tx + noise

    # --- 原子 + 组原子 ---
    atoms, grid = build_candidate_atoms(n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx)
    group_atoms, group_masks = build_group_atoms(atoms, grid, max_doppler, dirichlet_r)
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    true_set = set(np.where(true_mask)[0])

    res_pow = 1.0  # noise floor
    kappa_g = kappa * 0.85

    # --- 单径 OMP ---
    det_s, g_s = single_omp_detect(y, atoms, n_paths + 2, kappa, res_pow)
    H_est_s = _build_channel_from_paths(det_s, g_s, atoms, grid, c1, n_sc, dirichlet_r)
    ber_s, nmse_s = _eval_ber_nmse(H_est_s, H_true, y, x_tx, data_idx, sym_data, n_data, noise_power)
    pd_s = len(set(det_s) & true_set) / max(n_paths, 1)

    # --- 组级 OMP ---
    groups, det_g, n_iter, g_g = group_omp_detect(
        y, atoms, grid, group_atoms, group_masks, n_paths + 2, kappa_g, res_pow)
    H_est_g = _build_channel_from_paths(det_g, g_g, atoms, grid, c1, n_sc, dirichlet_r)
    ber_g, nmse_g = _eval_ber_nmse(H_est_g, H_true, y, x_tx, data_idx, sym_data, n_data, noise_power)
    pd_g = len(set(det_g) & true_set) / max(n_paths, 1)

    return {
        "ber_s": ber_s, "ber_g": ber_g,
        "nmse_s": nmse_s, "nmse_g": nmse_g,
        "pd_s": pd_s, "pd_g": pd_g,
        "n_det_s": len(det_s), "n_det_g": len(det_g),
    }


def _build_channel_from_paths(detected, detected_gains, atoms, grid, c1, n_sc, dirichlet_r):
    """从检测候选索引 + OMP 内积增益 → H_est。"""
    if not detected:
        return np.zeros((n_sc, n_sc), dtype=complex)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    H = np.zeros((n_sc, n_sc), dtype=complex)
    for idx in detected:
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        g_i = detected_gains.get(idx, 0.0)
        H += g_i * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)
    return H


def _eval_ber_nmse(H_est, H_true, y, x_tx, data_idx, sym_data, n_data, noise_power):
    """LMMSE 检测 → BER + NMSE。"""
    if np.all(H_est == 0):
        return 0.5, 0.0

    H_data = H_est[:, data_idx]
    reg = noise_power
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_est = w_mmse @ y

    bits_est = (np.real(x_est) > 0).astype(int)
    bits_true = (np.real(sym_data) > 0).astype(int)
    ber = float(np.mean(bits_est != bits_true))

    nmse_num = np.sum(np.abs(H_true - H_est) ** 2)
    nmse_den = max(np.sum(np.abs(H_true) ** 2), 1e-15)
    nmse_db = 10.0 * np.log10(nmse_num / nmse_den)
    return ber, nmse_db


# ======================================================================
# MC 批量
# ======================================================================

def run_grid(cfg, n_trials=100):
    """
    SNR × rho 网格扫描。

    返回:
      rows: list[dict] 每行一个 (snr, rho) 组合
    """
    snr_vec = cfg.get("snr_vec", [0, 5, 10, 15, 20, 25])
    rho_vec = cfg.get("rho_vec", [0.0, 0.5, 0.9])
    results = []

    total = len(snr_vec) * len(rho_vec)
    count = 0
    t_start = time.time()

    for rho_val in rho_vec:
        for snr_val in snr_vec:
            _c = dict(cfg)
            _c["snr_db"] = snr_val
            _c["rho"] = rho_val

            trials = []
            for t in range(n_trials):
                _c["seed"] = cfg.get("seed", 0) + count * n_trials + t
                r = run_one_trial(_c)
                trials.append(r)

            row = {
                "snr_db": snr_val,
                "rho": rho_val,
                "ber_s": np.mean([r["ber_s"] for r in trials]),
                "ber_g": np.mean([r["ber_g"] for r in trials]),
                "nmse_s": np.mean([r["nmse_s"] for r in trials]),
                "nmse_g": np.mean([r["nmse_g"] for r in trials]),
                "pd_s": np.mean([r["pd_s"] for r in trials]),
                "pd_g": np.mean([r["pd_g"] for r in trials]),
                "ber_s_std": np.std([r["ber_s"] for r in trials]),
                "ber_g_std": np.std([r["ber_g"] for r in trials]),
                "n_det_s": np.mean([r["n_det_s"] for r in trials]),
                "n_det_g": np.mean([r["n_det_g"] for r in trials]),
            }
            results.append(row)
            count += 1
            elapsed = time.time() - t_start
            print(f"  [{count}/{total}] SNR={snr_val:2d}dB rho={rho_val:.1f}  "
                  f"BER_s={row['ber_s']:.2e} BER_g={row['ber_g']:.2e}  "
                  f"Pd_s={row['pd_s']:.3f} Pd_g={row['pd_g']:.3f}  "
                  f"NMSE_s={row['nmse_s']:+.1f}dB NMSE_g={row['nmse_g']:+.1f}dB  "
                  f"({elapsed/count*(total-count):.0f}s left)")

    return results


# ======================================================================
# 绘图
# ======================================================================

def plot_end_to_end(results, out_dir):
    """BER / Pd / NMSE 三面板 × 单径/组双曲线。"""
    snr_vec = sorted(set(r["snr_db"] for r in results))
    rho_vec = sorted(set(r["rho"] for r in results))
    colors  = {0.0: "gray", 0.5: "steelblue", 0.9: "crimson"}
    markers = {0.0: "s", 0.5: "D", 0.9: "o"}

    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    for ri, rho_val in enumerate(rho_vec):
        rows = [r for r in results if r["rho"] == rho_val]
        snr = np.array([r["snr_db"] for r in rows])
        c = colors[rho_val]
        m = markers[rho_val]
        label_s = f"single rho={rho_val:.1f}"
        label_g = f"group  rho={rho_val:.1f}"

        # BER (row 0)
        ax = axes[0, 0]
        ax.semilogy(snr, [r["ber_s"] for r in rows], f'{m}-', color=c, linewidth=1.5, markersize=7, alpha=0.6, label=label_s)
        ax.semilogy(snr, [r["ber_g"] for r in rows], f'{m}--', color=c, linewidth=2, markersize=9, markerfacecolor=c, label=label_g)

        # Pd (row 0)
        ax = axes[0, 1]
        ax.plot(snr, [r["pd_s"] for r in rows], f'{m}-', color=c, linewidth=1.5, markersize=7, alpha=0.6)
        ax.plot(snr, [r["pd_g"] for r in rows], f'{m}--', color=c, linewidth=2, markersize=9, markerfacecolor=c)

        # NMSE (row 0)
        ax = axes[0, 2]
        ax.plot(snr, [r["nmse_s"] for r in rows], f'{m}-', color=c, linewidth=1.5, markersize=7, alpha=0.6)
        ax.plot(snr, [r["nmse_g"] for r in rows], f'{m}--', color=c, linewidth=2, markersize=9, markerfacecolor=c)

        # dPd (row 1, left)
        ax = axes[1, 0]
        dPd = np.array([r["pd_g"] - r["pd_s"] for r in rows])
        ax.bar(snr + (ri-1)*1.2, dPd, width=1.0, color=c, alpha=0.8, label=f"rho={rho_val:.1f}")

        # BER gain (row 1, middle)
        ax = axes[1, 1]
        bgain = np.array([r["ber_s"]/max(r["ber_g"], 1e-10) for r in rows])
        ax.plot(snr, bgain, f'{m}-', color=c, linewidth=1.5, markersize=7, label=f"rho={rho_val:.1f}")

        # dNMSE (row 1, right)
        ax = axes[1, 2]
        dN = np.array([r["nmse_g"] - r["nmse_s"] for r in rows])
        ax.bar(snr + (ri-1)*1.2, dN, width=1.0, color=c, alpha=0.8, label=f"rho={rho_val:.1f}")

    # 格式化
    axes[0,0].set_xlabel("Data SNR (dB)"); axes[0,0].set_ylabel("BER"); axes[0,0].set_title("BER"); axes[0,0].legend(fontsize=7); axes[0,0].grid(True, alpha=0.3)
    axes[0,1].set_xlabel("Data SNR (dB)"); axes[0,1].set_ylabel("P_d"); axes[0,1].set_title("Detection Probability"); axes[0,1].grid(True, alpha=0.3)
    axes[0,2].set_xlabel("Data SNR (dB)"); axes[0,2].set_ylabel("NMSE (dB)"); axes[0,2].set_title("Channel Estimation NMSE"); axes[0,2].grid(True, alpha=0.3)

    axes[1,0].set_xlabel("Data SNR (dB)"); axes[1,0].set_ylabel("dPd"); axes[1,0].set_title("Group OMP Pd Gain (dPd)");
    axes[1,0].axhline(0, color='k', linewidth=0.5); axes[1,0].legend(fontsize=7); axes[1,0].grid(True, alpha=0.3, axis='y')
    axes[1,1].set_xlabel("Data SNR (dB)"); axes[1,1].set_ylabel("BER_s / BER_g"); axes[1,1].set_title("BER Ratio (single/group, >1 = group better)");
    axes[1,1].axhline(1, color='k', linewidth=0.5); axes[1,1].legend(fontsize=7); axes[1,1].grid(True, alpha=0.3)
    axes[1,2].set_xlabel("Data SNR (dB)"); axes[1,2].set_ylabel("dNMSE (dB)"); axes[1,2].set_title("NMSE Improvement (negative = better)");
    axes[1,2].axhline(0, color='k', linewidth=0.5); axes[1,2].legend(fontsize=7); axes[1,2].grid(True, alpha=0.3, axis='y')

    fig.suptitle("End-to-End: Single-OMP vs Group-OMP under BDLR Correlated Channel",
                 fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(out_dir, "Fig6_EndToEnd.png")
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"\n  Fig6 saved: {path}")


def print_summary_table(results):
    """打印汇总表。"""
    snr_list = sorted(set(r["snr_db"] for r in results))
    print(f"\n{'SNR':>4} {'rho':>4} {'BER_s':>9} {'BER_g':>9} {'BER_imp':>8} "
          f"{'Pd_s':>6} {'Pd_g':>6} {'dPd':>6} {'NMSE_s':>7} {'NMSE_g':>7}")
    print("-" * 85)
    for r in sorted(results, key=lambda x: (x["rho"], x["snr_db"])):
        ber_ratio = r["ber_s"] / max(r["ber_g"], 1e-10)
        print(f"{r['snr_db']:4.0f} {r['rho']:4.1f} {r['ber_s']:9.2e} {r['ber_g']:9.2e} "
              f"{ber_ratio:8.1f}x {r['pd_s']:6.3f} {r['pd_g']:6.3f} {r['pd_g']-r['pd_s']:+6.3f} "
              f"{r['nmse_s']:+7.1f} {r['nmse_g']:+7.1f}")


# ======================================================================
# main
# ======================================================================

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="小规模验证 (20 trials)")
    ap.add_argument("--mc", type=int, default=100)
    ap.add_argument("--out-dir", type=str, default=None)
    args = ap.parse_args()

    n_trials = 20 if args.quick else args.mc
    sim_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = args.out_dir or os.path.join(sim_dir, "Results", "Figures")
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "p_fa": 1e-3, "seed": 0, "pilot_idx": 0,
        "kappa": 24.6,
        "snr_vec": [0, 5, 10, 15, 20, 25],
        "rho_vec": [0.0, 0.5, 0.9],
    }

    print(f"=== End-to-End: Single vs Group OMP (n_trials={n_trials}) ===\n")
    t0 = time.time()
    results = run_grid(cfg, n_trials)
    elapsed = time.time() - t0
    print(f"\n  Total: {elapsed:.0f}s ({elapsed/n_trials/len(results):.1f}s/trial)")

    print_summary_table(results)
    plot_end_to_end(results, out_dir)

    # CSV
    csv_path = os.path.join(out_dir, "end_to_end_results.csv")
    with open(csv_path, "w") as f:
        keys = ["snr_db","rho","ber_s","ber_g","pd_s","pd_g","nmse_s","nmse_g","ber_s_std","ber_g_std","n_det_s","n_det_g"]
        f.write(",".join(keys) + "\n")
        for r in results:
            f.write(",".join(str(r[k]) for k in keys) + "\n")
    print(f"  CSV saved: {csv_path}")
