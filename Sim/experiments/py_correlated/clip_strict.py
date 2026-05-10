"""
CLIP 严格对齐实现 — 三阶段独立函数，逐一对应 MATLAB GiFreeReceiver。

对照关系:
  stage1_coarse_lock   ↔ GiFreeReceiver Phase 1 (bootstrap + CLIP coarse)
  stage2_iteration     ↔ GiFreeReceiver Phase 2 (iterative refinement)
  stage3_polish        ↔ GiFreeReceiver Phase 3 (final polishing)

每阶段内:
  - 渐进 CFAR 严格度 (progressiveStrictness: s_max → s_min 线性衰减)
  - 路径稳定度门控 (pathStabilityGate: τ_s = 2 次连续检测)
  - 置信门控 (confidenceGate: top-ρ 符号用于 ID2P 消除)
  - 残差功率跟踪 (residualPower: 每次迭代更新)

用法:
  python clip_strict.py --mc 30
  python clip_strict.py --quick
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
from group_omp import build_group_atoms


# ======================================================================
# 公共原子: LMMSE, 增益估计, 残差功率
# ======================================================================

def lmmse_detect_with_confidence(y, H_est, data_idx, n_data, noise_power,
                                  pilot_idx=0, pilot_amp=None,
                                  residual_interf=0.0, data_amp_sq=1.0):
    """
    LMMSE 数据检测, 返回软符号 + 逐符号可靠性。

    先减去 pilot 响应, 再做 LMMSE。
    正则化系数 reg = (noise_power + residual_interf) / data_amp² 对齐 MATLAB
    GiFreeReceiver 的 residualInterf 逻辑。
    """
    # 消除 pilot 响应
    if pilot_amp is not None and pilot_amp > 0:
        y_eff = y - H_est[:, pilot_idx] * pilot_amp
    else:
        y_eff = y

    H_data = H_est[:, data_idx]
    # reg: 对齐 MATLAB ompRegParam = residualInterf / pilotAmpSq
    # 对 LMMSE data detection: reg = (σ² + ID2P_residual) / data_amp²
    reg = max(noise_power + residual_interf, noise_power) / max(data_amp_sq, 1.0)
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_soft = w_mmse @ y_eff

    # 可靠性: 距 QPSK 星座点的归一化距离的倒数
    # 高 SNR → x_soft 接近星座点 → reliability 大
    dist_to_const = np.minimum(
        np.abs(x_soft - (1+1j)/np.sqrt(2)),
        np.abs(x_soft - (-1+1j)/np.sqrt(2))
    )
    dist_to_const = np.minimum(dist_to_const,
        np.minimum(np.abs(x_soft - (1-1j)/np.sqrt(2)),
                   np.abs(x_soft - (-1-1j)/np.sqrt(2))))
    reliability = 1.0 / (dist_to_const + 0.01)

    return x_soft, reliability


def estimate_gains_pilot_constrained(detected_indices, atoms, grid, y, n_sc, c1):
    """
    Pilot 列联合 LS 增益估计。

    用每个候选的 pilot 列 phi[:,0]*pilot_amp 做联合 LS:
      [phi_1[:,0], phi_2[:,0], ...] @ g ≈ y

    比 pilot 行内积正确: pilot 行内积对整数 Doppler 候选退化为
    单个数据子载波的值, 与 pilot 响应无关。
    """
    if not detected_indices:
        return {}
    pilot_amp_est = np.sqrt(10.0 ** 3.5)  # ≈ 56.2
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    K = len(detected_indices)
    A = np.zeros((n_sc, K), dtype=complex)
    for ki, idx in enumerate(detected_indices):
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        from channel import build_path_matrix as _bpm
        phi = _bpm(d_est, nu_est, n_sc, loc_step, 3)  # dirichlet_r=3
        A[:, ki] = phi[:, 0] * pilot_amp_est
    g_vec, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    return {int(idx): g_vec[ki] for ki, idx in enumerate(detected_indices)}


def estimate_residual_power(r):
    """残差功率: mean(|r|²)。"""
    return float(np.mean(np.abs(r) ** 2))


def _build_H_from_detected(detected, gains, grid, c1, n_sc, dirichlet_r):
    """从 detected + gains → H_est。"""
    if not detected:
        return np.zeros((n_sc, n_sc), dtype=complex)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    H = np.zeros((n_sc, n_sc), dtype=complex)
    for idx in detected:
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        g_i = gains.get(int(idx), 0.0)
        H += g_i * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)
    return H


# ======================================================================
# Stage 1: 粗锁定 (Coarse-Lock)
# ======================================================================

def stage1_coarse_lock(y, x_tx, H_true_shape, atoms, grid,
                       group_atoms, group_masks, cfg):
    """
    CLIP Stage 1: 渐进 CFAR 阈值 + 路径稳定度门控 + 置信门控。

    参数 (对齐 MATLAB GiFreeConfig):
      cfg:
        progressive_init_scale:  float 渐进严格度起始值 (默认 3.0)
        progressive_final_scale: float 渐进严格度最终值 (默认 1.0)
        num_phase1_iter:         int   阶段1迭代次数 (默认 5)
        stability_threshold:     int   路径稳定度门控阈值 (默认 2次)
        stage1_keep_ratio:       float 置信门控保真比例 (默认 0.30)
        kappa_ind:               float 公式型 κ_ind = ln(M/P_fa)
        max_paths:               int   最大路径数

    返回:
      H_est:          (N,N) 初始信道估计
      stable_paths:   list[int] 稳定路径集
      data_soft:      (n_data,) 数据软估计
    """
    n_sc         = cfg["n_sc"]
    n_paths      = cfg["n_paths"]
    max_delay    = cfg["max_delay"]
    max_doppler  = cfg["max_doppler"]
    dirichlet_r  = cfg["dirichlet_r"]
    pilot_idx    = cfg["pilot_idx"]
    c1           = cfg["_c1"]
    noise_power  = 1.0

    s_max        = float(cfg.get("progressive_init_scale", 3.0))
    s_min        = float(cfg.get("progressive_final_scale", 1.0))
    n_iter       = int(cfg.get("num_phase1_iter", 5))
    tau_s        = int(cfg.get("stability_threshold", 2))
    keep_ratio   = float(cfg.get("stage1_keep_ratio", 0.30))
    kappa_ind    = float(cfg.get("kappa_ind", 9.9))

    data_amp_sq = float(cfg.get("data_amp_sq", 1.0))

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)

    # 初始化
    M = grid.shape[0]
    all_detected = set()
    path_histogram = {}   # idx → 连续检测次数
    stable_paths = set()
    H_est = np.zeros((n_sc, n_sc), dtype=complex)
    r_work = y.copy()
    gains = {}

    for iteration in range(n_iter):
        # 渐进严格度: s = s_max + (s_min - s_max) * iter / (n_iter - 1)
        progress = iteration / max(n_iter - 1, 1)
        s = s_max + (s_min - s_max) * progress
        kappa_eff = s * kappa_ind

        # --- 残差功率: 渐进方式 (早期用噪声底, 后期用实测 residualInterf) ---
        if iteration <= 1:
            res_pow = noise_power
        else:
            # 迭代 ≥2: 用稳定路径估计的 H_est 消去 pilot 响应, 余下为 residualInterf + noise
            x_pilot_only = np.zeros(n_sc, dtype=complex)
            x_pilot_only[pilot_idx] = x_tx[pilot_idx]
            pilot_residual = r_work - H_est @ x_pilot_only
            res_pow_measured = float(np.mean(np.abs(pilot_residual)**2))
            # 防止过大: 截断到 [noise_power, 10*noise_power]
            res_pow = float(np.clip(res_pow_measured, noise_power, 10.0 * noise_power))
        threshold = kappa_eff * res_pow

        # --- OMP 检测 (group-OMP with current threshold) ---
        detected_this_round = []
        r_local = r_work.copy()
        for _ in range(n_paths + 2):
            scores = np.abs(group_atoms @ r_local.conj()) ** 2
            best = int(np.argmax(scores))
            if float(scores[best]) < threshold:
                break
            members = list(group_masks[best])
            A_g = atoms[members, :].T
            gv, _, _, _ = np.linalg.lstsq(A_g, r_local, rcond=None)
            r_local -= A_g @ gv
            detected_this_round.append(best)
            for mi in members:
                all_detected.add(int(mi))

        # --- 路径稳定度门控 ---
        for idx in detected_this_round:
            path_histogram[idx] = path_histogram.get(idx, 0) + 1
            if path_histogram[idx] >= tau_s:
                stable_paths.add(idx)

        # --- 重建 H_est (pilot 列联合 LS 增益估计) ---
        stable_list = list(stable_paths) if stable_paths else detected_this_round
        gains = estimate_gains_pilot_constrained(stable_list, atoms, grid, y, n_sc, c1)
        H_est = _build_H_from_detected(stable_list, gains, grid, c1, n_sc, dirichlet_r)

        # --- 置信门控 + ID2P 消除 (Stage 1 只做一次, 在最后) ---
        if iteration == n_iter - 1:
            # LMMSE 在原始 y 上做 (减去 pilot 响应), r_work 仅用于 OMP
            # residual_interf: 实测 pilot 位置残差 - noise
            y_pilot_only = np.zeros(n_sc, dtype=complex)
            y_pilot_only[pilot_idx] = x_tx[pilot_idx]
            pilot_resid_local = y - H_est @ y_pilot_only
            res_interf = max(float(np.mean(np.abs(pilot_resid_local)**2)) - noise_power, 0.0)
            x_soft, reliability = lmmse_detect_with_confidence(
                y, H_est, data_idx, n_data, noise_power,
                pilot_idx=pilot_idx, pilot_amp=np.sqrt(10.0**3.5),
                residual_interf=res_interf, data_amp_sq=data_amp_sq)
            # 按可靠性排序, 保留 top keep_ratio
            rank = np.argsort(reliability)[::-1]
            n_keep = int(n_data * keep_ratio)
            keep_mask = np.zeros(n_data, dtype=bool)
            keep_mask[rank[:n_keep]] = True

            # QPSK 硬判决 (仅对保留符号)
            x_hard = _qpsk_hard(x_soft)
            x_fb = np.zeros(n_data, dtype=complex)
            x_fb[keep_mask] = x_hard[keep_mask]

            # ID2P 消除
            id2p = H_est[:, data_idx] @ x_fb
            r_work = y - id2p
        else:
            # 中间迭代: LMMSE 在原始 y 上做
            x_soft, _ = lmmse_detect_with_confidence(
                y, H_est, data_idx, n_data, noise_power,
                pilot_idx=pilot_idx, pilot_amp=np.sqrt(10.0**3.5),
                data_amp_sq=data_amp_sq)
            x_hard = _qpsk_hard(x_soft)
            id2p = H_est[:, data_idx] @ x_hard
            r_work = y - id2p

    # 最终稳定路径集
    final_paths = list(stable_paths) if stable_paths else list(all_detected)[:n_paths]
    # 从初始 y 重估增益
    gains_final = estimate_gains_pilot_constrained(final_paths, atoms, grid, y,
                                                    n_sc, c1)
    H_final = _build_H_from_detected(final_paths, gains_final, grid, c1, n_sc, dirichlet_r)

    # Stage 1 最终 LMMSE (含 residual_interf 估计)
    y_pilot_only = np.zeros(n_sc, dtype=complex)
    y_pilot_only[pilot_idx] = x_tx[pilot_idx]
    pilot_resid_final = y - H_final @ y_pilot_only
    res_interf_final = max(float(np.mean(np.abs(pilot_resid_final)**2)) - noise_power, 0.0)
    x_soft_final, _ = lmmse_detect_with_confidence(
        y, H_final, data_idx, n_data, noise_power,
        pilot_idx=pilot_idx, pilot_amp=np.sqrt(10.0**3.5),
        residual_interf=res_interf_final, data_amp_sq=data_amp_sq)

    return H_final, final_paths, x_soft_final, len(stable_paths)


# ======================================================================
# Stage 2: 迭代精炼 (Iteration / Refinement)
# ======================================================================

def stage2_iteration(y, x_tx, H_init, stable_paths, atoms, grid,
                     group_atoms, group_masks, cfg):
    """
    CLIP Stage 2: 软判决 + 渐进置信门控 + ID2P 消除 + OMP 重检测。

    参数:
      cfg:
        num_phase2_iter:     int   阶段2迭代次数 (默认 3)
        stage2_keep_min:     float 起始保真比 (默认 0.10)
        stage2_keep_max:     float 最终保真比 (默认 0.80)
        kappa_ind:           float
        max_paths:           int
        enable_residual_inj: bool  是否残差注入 (默认 True)

    返回:
      H_est:       (N,N) 精炼信道估计
      all_paths:   list[int] 更新路径集
      data_soft:   (n_data,) 数据软估计
    """
    n_sc         = cfg["n_sc"]
    n_paths      = cfg["n_paths"]
    dirichlet_r  = cfg["dirichlet_r"]
    pilot_idx    = cfg["pilot_idx"]
    c1           = cfg["_c1"]
    noise_power  = 1.0

    n_iter       = int(cfg.get("num_phase2_iter", 3))
    keep_min     = float(cfg.get("stage2_keep_min", 0.10))
    keep_max     = float(cfg.get("stage2_keep_max", 0.80))
    kappa_ind    = float(cfg.get("kappa_ind", 9.9))
    enable_inj   = cfg.get("enable_residual_inj", True)
    max_doppler  = cfg["max_doppler"]

    data_amp_sq = float(cfg.get("data_amp_sq", 1.0))

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)

    H_est = H_init.copy()
    r_work = y.copy()
    current_paths = list(stable_paths)

    for iteration in range(n_iter):
        # 渐进保真比
        keep_ratio = keep_min + (keep_max - keep_min) * iteration / max(n_iter - 1, 1)

        # --- LMMSE 在原始 y 上做 (减去 pilot), r_work 仅用于 OMP ---
        y_pilot_only_s2 = np.zeros(n_sc, dtype=complex)
        y_pilot_only_s2[pilot_idx] = x_tx[pilot_idx]
        pilot_resid_s2 = y - H_est @ y_pilot_only_s2
        res_interf_s2 = max(float(np.mean(np.abs(pilot_resid_s2)**2)) - noise_power, 0.0)
        x_soft, reliability = lmmse_detect_with_confidence(
            y, H_est, data_idx, n_data, noise_power,
            pilot_idx=pilot_idx, pilot_amp=np.sqrt(10.0**3.5),
            residual_interf=res_interf_s2, data_amp_sq=data_amp_sq)
        rank = np.argsort(reliability)[::-1]
        n_keep = int(n_data * keep_ratio)
        keep_mask = np.zeros(n_data, dtype=bool)
        keep_mask[rank[:n_keep]] = True

        # QPSK 硬判决 (仅对保留符号)
        x_hard = _qpsk_hard(x_soft)
        x_fb = np.zeros(n_data, dtype=complex)
        x_fb[keep_mask] = x_hard[keep_mask]

        # --- ID2P 消除 ---
        id2p = H_est[:, data_idx] @ x_fb
        y_clean = y - id2p

        # --- OMP 重检测 (residualInterf 估计, 截断防爆) ---
        pilot_residual_s2 = r_work - H_est[:, pilot_idx] * np.sqrt(10.0**3.5)
        res_pow_measured = float(np.mean(np.abs(pilot_residual_s2)**2))
        res_pow_s2 = float(np.clip(res_pow_measured, noise_power, 10.0 * noise_power))
        threshold = kappa_ind * res_pow_s2
        r_local = y_clean.copy()
        new_detected = []
        for _ in range(n_paths + 2):
            scores = np.abs(group_atoms @ r_local.conj()) ** 2
            best = int(np.argmax(scores))
            if float(scores[best]) < threshold:
                break
            members = list(group_masks[best])
            A_g = atoms[members, :].T
            gv, _, _, _ = np.linalg.lstsq(A_g, r_local, rcond=None)
            r_local -= A_g @ gv
            new_detected.append(best)

        # 合并新旧路径集
        all_paths = list(set(current_paths + new_detected))
        # 裁剪到 max_paths+2
        if len(all_paths) > n_paths + 2:
            # 按匹配得分排序裁剪
            scores_full = np.abs(atoms[all_paths, :] @ y_clean.conj()) ** 2
            keep_order = np.argsort(scores_full)[::-1][:n_paths + 2]
            all_paths = [all_paths[i] for i in keep_order]

        # --- 残差注入 (额外候选) ---
        if enable_inj:
            r_local2 = y_clean.copy()
            # 移除已检测路径
            for idx in all_paths:
                atom = atoms[idx, :]
                gi = np.dot(atom.conj(), r_local2)
                r_local2 -= gi * atom
            # 在残差中搜索最强额外候选
            extra_scores = np.abs(atoms @ r_local2.conj()) ** 2
            extra_best = int(np.argmax(extra_scores))
            if float(extra_scores[extra_best]) > threshold * 0.5:
                all_paths.append(extra_best)

        current_paths = all_paths

        # --- 重估增益 + 重建 H_est ---
        gains = estimate_gains_pilot_constrained(current_paths, atoms, grid,
                                                  y_clean, n_sc, c1)
        H_est = _build_H_from_detected(current_paths, gains, grid, c1, n_sc, dirichlet_r)

        # --- 更新工作信号 (用全数据硬判决重新消除) ---
        x_hard_all = _qpsk_hard(x_soft)
        id2p_all = H_est[:, data_idx] @ x_hard_all
        r_work = y - id2p_all

        # --- Stage 2 自适应: 下一轮的 res_interf 将从 r_work 重估 ---
        # (在循环顶部 pilot_resid_s2 = y - H_est @ y_pilot 自动反映改善)

    # 最终增益 (从原始 y 重估)
    gains_final = estimate_gains_pilot_constrained(current_paths, atoms, grid,
                                                    y, n_sc, c1)
    H_final = _build_H_from_detected(current_paths, gains_final, grid, c1, n_sc, dirichlet_r)
    # Stage 2 最终 LMMSE
    y_pilot_only_s2f = np.zeros(n_sc, dtype=complex)
    y_pilot_only_s2f[pilot_idx] = x_tx[pilot_idx]
    pilot_resid_s2f = y - H_final @ y_pilot_only_s2f
    res_interf_s2f = max(float(np.mean(np.abs(pilot_resid_s2f)**2)) - noise_power, 0.0)
    x_soft_final, _ = lmmse_detect_with_confidence(
        y, H_final, data_idx, n_data, noise_power,
        pilot_idx=pilot_idx, pilot_amp=np.sqrt(10.0**3.5),
        residual_interf=res_interf_s2f, data_amp_sq=data_amp_sq)

    return H_final, current_paths, x_soft_final


# ======================================================================
# Stage 3: 最终定型 (Polishing)
# ======================================================================

def stage3_polish(y, x_tx, H_in, fixed_paths, atoms, grid, cfg):
    """
    CLIP Stage 3: 固定支撑集 → 联合 LS 增益 → 截尾噪声 → 最终 LMMSE。

    参数:
      cfg:
        trim_alpha: float 截尾比例 (默认 0.15)
        max_paths:  int

    返回:
      ber:  float 最终 BER
      nmse: float 最终 NMSE (需外部传入 H_true)
      H_final: (N,N)
      noise_est: float
    """
    n_sc         = cfg["n_sc"]
    n_paths      = cfg["n_paths"]
    dirichlet_r  = cfg["dirichlet_r"]
    pilot_idx    = cfg["pilot_idx"]
    c1           = cfg["_c1"]
    noise_power  = 1.0
    trim_alpha   = float(cfg.get("trim_alpha", 0.15))

    data_amp_sq = float(cfg.get("data_amp_sq", 1.0))

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)

    # --- 联合 LS 增益重估计 (用 pilot-only 观测) ---
    gains = estimate_gains_pilot_constrained(fixed_paths, atoms, grid, y,
                                              n_sc, c1)
    H_final = _build_H_from_detected(fixed_paths, gains, grid, c1, n_sc, dirichlet_r)

    # --- 截尾噪声估计 ---
    r_final = y - H_final @ x_tx
    residual_powers = np.sort(np.abs(r_final) ** 2)
    n_trim = int(n_sc * trim_alpha)
    n_keep = n_sc - n_trim
    noise_est = float(np.mean(residual_powers[:n_keep]))

    # --- 最终 LMMSE (先减去 pilot 响应, reg 含 residual_interf) ---
    pilot_amp_val = np.sqrt(10.0**3.5)
    y_eff = y - H_final[:, pilot_idx] * pilot_amp_val
    H_data = H_final[:, data_idx]
    # Stage 3 reg: noise_est 即 residual + noise 的截尾估计
    data_amp_sq = float(cfg.get("data_amp_sq", 1.0))
    reg = max(noise_est, noise_power) / max(data_amp_sq, 1.0)
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_final = w_mmse @ y_eff

    return H_final, x_final, noise_est


# ======================================================================
# 完整 CLIP 流水线
# ======================================================================

def run_clip_trial(cfg):
    """
    单次 CLIP 三阶段 trial。

    返回三条路径的 BER/NMSE/Pd:
      S1:     Stage 1 only
      S1+S2:  Stage 1 + Stage 2
      S1+S2+S3: Full CLIP (Stage 1 + 2 + 3)
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
    p_fa         = float(cfg.get("p_fa", 1e-3))

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)

    # --- 信道 + 帧 ---
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    gains_true = bdlr_gains(n_paths, n_clusters, rho, rng)
    H_true = compute_effective_channel(delays, dopplers, gains_true, n_sc, c1, dirichlet_r)

    data_amp_sq = float(cfg.get("data_amp_sq", 1.0))

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)
    data_amp = np.sqrt(10.0 ** (snr_db / 10.0))

    # 导频幅度: 对齐 MATLAB pilot SNR = 35dB → amp = sqrt(10^3.5)
    pilot_amp = np.sqrt(10.0 ** (35.0 / 10.0))  # ≈ 56.2
    x_tx = np.zeros(n_sc, dtype=complex)
    x_tx[pilot_idx] = pilot_amp * (1.0 + 0j)
    bits = rng.integers(0, 2, 2 * n_data)
    sym_data = ((2*bits[:n_data]-1) + 1j*(2*bits[n_data:]-1)) / np.sqrt(2.0)
    x_tx[data_idx] = data_amp * sym_data
    # 噪声功率 = 1.0, 对齐 MATLAB
    noise = (rng.standard_normal(n_sc) + 1j*rng.standard_normal(n_sc)) / np.sqrt(2.0)
    y = H_true @ x_tx + noise

    # --- 原子 ---
    atoms, grid = build_candidate_atoms(n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx)
    group_atoms, group_masks = build_group_atoms(atoms, grid, max_doppler, dirichlet_r)
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    true_set = set(np.where(true_mask)[0])

    M = grid.shape[0]
    kappa_ind = np.log(M / p_fa)
    data_amp_sq_val = data_amp ** 2

    s1_cfg = dict(cfg)
    s1_cfg.update({
        "_c1": c1, "n_sc": n_sc, "n_paths": n_paths,
        "max_delay": max_delay, "max_doppler": max_doppler,
        "dirichlet_r": dirichlet_r, "pilot_idx": pilot_idx,
        "kappa_ind": kappa_ind,
        "data_amp_sq": data_amp_sq_val,
        "progressive_init_scale": float(cfg.get("progressive_init_scale", 3.0)),
        "progressive_final_scale": float(cfg.get("progressive_final_scale", 1.0)),
        "num_phase1_iter": int(cfg.get("num_phase1_iter", 6)),
        "stability_threshold": int(cfg.get("stability_threshold", 2)),
        "stage1_keep_ratio": float(cfg.get("stage1_keep_ratio", 0.30)),
    })

    s2_cfg = dict(cfg)
    s2_cfg.update({
        "_c1": c1, "n_sc": n_sc, "n_paths": n_paths,
        "max_delay": max_delay, "max_doppler": max_doppler,
        "dirichlet_r": dirichlet_r, "pilot_idx": pilot_idx,
        "kappa_ind": kappa_ind,
        "data_amp_sq": data_amp_sq_val,
        "num_phase2_iter": int(cfg.get("num_phase2_iter", 4)),
        "stage2_keep_min": float(cfg.get("stage2_keep_min", 0.10)),
        "stage2_keep_max": float(cfg.get("stage2_keep_max", 0.70)),
        "enable_residual_inj": cfg.get("enable_residual_inj", True),
    })

    s3_cfg = dict(cfg)
    s3_cfg.update({
        "_c1": c1, "n_sc": n_sc, "n_paths": n_paths,
        "max_delay": max_delay, "max_doppler": max_doppler,
        "dirichlet_r": dirichlet_r, "pilot_idx": pilot_idx,
        "data_amp_sq": data_amp_sq_val,
        "trim_alpha": float(cfg.get("trim_alpha", 0.15)),
    })

    # ===== Stage 1 only =====
    H_s1, paths_s1, x_s1, n_stable = stage1_coarse_lock(
        y, x_tx, (n_sc, n_sc), atoms, grid, group_atoms, group_masks, s1_cfg)
    ber_s1 = float(np.mean((np.real(x_s1) > 0).astype(int) !=
                            (np.real(sym_data) > 0).astype(int)))
    nmse_s1 = 10.0 * np.log10(np.sum(np.abs(H_true - H_s1)**2) /
                               max(np.sum(np.abs(H_true)**2), 1e-15))
    pd_s1 = len(set(paths_s1) & true_set) / max(n_paths, 1)

    # ===== Stage 1 + 2 =====
    H_s2, paths_s2, x_s2 = stage2_iteration(
        y, x_tx, H_s1, paths_s1, atoms, grid, group_atoms, group_masks, s2_cfg)
    ber_s2 = float(np.mean((np.real(x_s2) > 0).astype(int) !=
                            (np.real(sym_data) > 0).astype(int)))
    nmse_s2 = 10.0 * np.log10(np.sum(np.abs(H_true - H_s2)**2) /
                               max(np.sum(np.abs(H_true)**2), 1e-15))
    pd_s2 = len(set(paths_s2) & true_set) / max(n_paths, 1)

    # ===== Stage 1 + 2 + 3 =====
    H_s3, x_s3, noise_est = stage3_polish(
        y, x_tx, H_s2, paths_s2, atoms, grid, s3_cfg)
    ber_s3 = float(np.mean((np.real(x_s3) > 0).astype(int) !=
                            (np.real(sym_data) > 0).astype(int)))
    nmse_s3 = 10.0 * np.log10(np.sum(np.abs(H_true - H_s3)**2) /
                               max(np.sum(np.abs(H_true)**2), 1e-15))
    pd_s3 = len(set(paths_s2) & true_set) / max(n_paths, 1)  # Stage 3 不改变路径集

    return {
        "snr": snr_db, "rho": rho,
        "ber_s1": ber_s1, "ber_s2": ber_s2, "ber_s3": ber_s3,
        "nmse_s1": nmse_s1, "nmse_s2": nmse_s2, "nmse_s3": nmse_s3,
        "pd_s1": pd_s1, "pd_s2": pd_s2, "pd_s3": pd_s3,
        "n_s1": len(paths_s1), "n_s2": len(paths_s2), "n_stable": n_stable,
    }


# ======================================================================
# 辅助
# ======================================================================

def _qpsk_hard(x_soft):
    return ((np.real(x_soft) > 0).astype(float)*2-1 +
            1j*(np.imag(x_soft) > 0).astype(float)*2-1j) / np.sqrt(2.0)


# ======================================================================
# 批量 + 绘图
# ======================================================================

def run_grid(cfg, n_trials=30):
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
                trials.append(run_clip_trial(_c))

            row = {"snr_db": snr_val, "rho": rho_val}
            for key in ["ber_s1","ber_s2","ber_s3","nmse_s1","nmse_s2","nmse_s3",
                        "pd_s1","pd_s2","pd_s3"]:
                row[key] = np.mean([r[key] for r in trials])
            results.append(row)
            count += 1
            elapsed = time.time() - t0
            print(f"  [{count}/{total}] SNR={snr_val:2d}dB rho={rho_val:.1f}  "
                  f"BER: {row['ber_s1']:.3f}→{row['ber_s2']:.3f}→{row['ber_s3']:.3f}  "
                  f"PD: {row['pd_s1']:.3f}→{row['pd_s2']:.3f}  "
                  f"({elapsed/count*(total-count):.0f}s left)")

    return results


def plot_strict_clip(results, out_dir):
    snr_vec = sorted(set(r["snr_db"] for r in results))
    rho_vec = sorted(set(r["rho"] for r in results))
    colors = {0.0: "gray", 0.5: "steelblue", 0.9: "crimson"}

    fig, axes = plt.subplots(2, 3, figsize=(22, 14))

    for rho_val in rho_vec:
        rows = [r for r in results if r["rho"] == rho_val]
        snr = np.array([r["snr_db"] for r in rows])
        c = colors[rho_val]

        # BER
        ax = axes[0,0]
        ax.semilogy(snr, [r["ber_s1"] for r in rows], 's--', color=c, linewidth=1, alpha=0.5,
                     label=f"S1 rho={rho_val:.1f}")
        ax.semilogy(snr, [r["ber_s2"] for r in rows], 'D-.', color=c, linewidth=1.5, alpha=0.7,
                     label=f"S2 rho={rho_val:.1f}")
        ax.semilogy(snr, [r["ber_s3"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c,
                     label=f"S3 rho={rho_val:.1f}")

        # Pd
        ax = axes[0,1]
        ax.plot(snr, [r["pd_s1"] for r in rows], 's--', color=c, linewidth=1, alpha=0.5)
        ax.plot(snr, [r["pd_s2"] for r in rows], 'D-.', color=c, linewidth=1.5, alpha=0.7)
        ax.plot(snr, [r["pd_s3"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c)

        # NMSE
        ax = axes[0,2]
        ax.plot(snr, [r["nmse_s1"] for r in rows], 's--', color=c, linewidth=1, alpha=0.5)
        ax.plot(snr, [r["nmse_s2"] for r in rows], 'D-.', color=c, linewidth=1.5, alpha=0.7)
        ax.plot(snr, [r["nmse_s3"] for r in rows], 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c)

        # Pd gain (row 1)
        ax = axes[1,0]
        dpd_s2 = np.array([r["pd_s2"]-r["pd_s1"] for r in rows])
        dpd_s3 = np.array([r["pd_s3"]-r["pd_s1"] for r in rows])
        w = 1.0
        x_shift = (rho_vec.index(rho_val) - 1) * w * 0.3
        ax.bar(snr + x_shift - w*0.15, dpd_s2, width=0.3, color=c, alpha=0.5, label=f"S2 gain rho={rho_val:.1f}")
        ax.bar(snr + x_shift + w*0.15, dpd_s3, width=0.3, color=c, alpha=0.9, label=f"S3 gain rho={rho_val:.1f}")

        # BER gain
        ax = axes[1,1]
        bg_s2 = np.array([r["ber_s1"]/max(r["ber_s2"],1e-10) for r in rows])
        bg_s3 = np.array([r["ber_s1"]/max(r["ber_s3"],1e-10) for r in rows])
        ax.plot(snr, bg_s2, 'D-.', color=c, linewidth=1.5, alpha=0.7, label=f"S2 rho={rho_val:.1f}")
        ax.plot(snr, bg_s3, 'o-', color=c, linewidth=2, markersize=8, markerfacecolor=c, label=f"S3 rho={rho_val:.1f}")

        # NMSE improvement
        ax = axes[1,2]
        dn_s2 = np.array([r["nmse_s2"]-r["nmse_s1"] for r in rows])
        dn_s3 = np.array([r["nmse_s3"]-r["nmse_s1"] for r in rows])
        ax.bar(snr + x_shift - w*0.15, dn_s2, width=0.3, color=c, alpha=0.5)
        ax.bar(snr + x_shift + w*0.15, dn_s3, width=0.3, color=c, alpha=0.9)

    axes[0,0].set_xlabel("Data SNR (dB)"); axes[0,0].set_ylabel("BER"); axes[0,0].set_title("BER: S1 → S2 → S3"); axes[0,0].legend(fontsize=6); axes[0,0].grid(True, alpha=0.3)
    axes[0,1].set_xlabel("Data SNR (dB)"); axes[0,1].set_ylabel("P_d"); axes[0,1].set_title("Detection Probability"); axes[0,1].grid(True, alpha=0.3)
    axes[0,2].set_xlabel("Data SNR (dB)"); axes[0,2].set_ylabel("NMSE (dB)"); axes[0,2].set_title("Channel Estimation NMSE"); axes[0,2].grid(True, alpha=0.3)
    axes[1,0].set_xlabel("Data SNR (dB)"); axes[1,0].set_ylabel("dPd"); axes[1,0].set_title("Pd Gain (over Stage 1)"); axes[1,0].axhline(0,color='k',linewidth=0.5); axes[1,0].legend(fontsize=6); axes[1,0].grid(True,alpha=0.3,axis='y')
    axes[1,1].set_xlabel("Data SNR (dB)"); axes[1,1].set_ylabel("BER ratio"); axes[1,1].set_title("BER_s1/BER ( >1 = better)"); axes[1,1].axhline(1,color='k',linewidth=0.5); axes[1,1].legend(fontsize=6); axes[1,1].grid(True,alpha=0.3)
    axes[1,2].set_xlabel("Data SNR (dB)"); axes[1,2].set_ylabel("dNMSE (dB)"); axes[1,2].set_title("NMSE Improvement"); axes[1,2].axhline(0,color='k',linewidth=0.5); axes[1,2].legend(fontsize=6); axes[1,2].grid(True,alpha=0.3,axis='y')

    fig.suptitle("CLIP Strict: Stage 1 → Stage 2 → Stage 3 Progression",
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0,0,1,0.93])
    path = os.path.join(out_dir, "Fig8_CLIP_Strict.png")
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f"\n  Fig8 saved: {path}")


def print_table(results):
    print(f"\n{'SNR':>4} {'rho':>4} {'BER_s1':>8} {'BER_s2':>8} {'BER_s3':>8} "
          f"{'Pd_s1':>6} {'Pd_s2':>6} {'Pd_s3':>6} {'NMSE_s1':>7} {'NMSE_s2':>7} {'NMSE_s3':>7}")
    print("-" * 95)
    for r in sorted(results, key=lambda x: (x["rho"], x["snr_db"])):
        print(f"{r['snr_db']:4.0f} {r['rho']:4.1f} {r['ber_s1']:8.4f} {r['ber_s2']:8.4f} {r['ber_s3']:8.4f} "
              f"{r['pd_s1']:6.3f} {r['pd_s2']:6.3f} {r['pd_s3']:6.3f} "
              f"{r['nmse_s1']:+7.1f} {r['nmse_s2']:+7.1f} {r['nmse_s3']:+7.1f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--mc", type=int, default=30)
    ap.add_argument("--out-dir", type=str, default=None)
    args = ap.parse_args()

    n_trials = 10 if args.quick else args.mc
    sim_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = args.out_dir or os.path.join(sim_dir, "Results", "Figures")
    os.makedirs(out_dir, exist_ok=True)

    cfg = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "p_fa": 1e-3, "seed": 0, "pilot_idx": 0,
        # Stage 1
        "progressive_init_scale": 3.0,
        "progressive_final_scale": 1.0,
        "num_phase1_iter": 6,
        "stability_threshold": 2,
        "stage1_keep_ratio": 0.30,
        # Stage 2
        "num_phase2_iter": 4,
        "stage2_keep_min": 0.10,
        "stage2_keep_max": 0.70,
        "enable_residual_inj": True,
        # Stage 3
        "trim_alpha": 0.15,
        # Scan
        "snr_vec": [10, 15, 20, 25],
        "rho_vec": [0.0, 0.5, 0.9],
    }

    print(f"=== CLIP Strict: Three-Stage Progression (n={n_trials}) ===\n")
    t0 = time.time()
    results = run_grid(cfg, n_trials)
    print(f"\n  Total: {time.time()-t0:.0f}s")
    print_table(results)
    plot_strict_clip(results, out_dir)
