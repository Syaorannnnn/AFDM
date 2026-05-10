"""
组级 OMP：每轮迭代选最佳 Doppler 邻域组而非最佳单候选。

对比:
  单径 OMP:  best = argmax_m |atom_m^H r|²
  组级 OMP:  best = argmax_G score(r, Phi_G, R_G)
             其中 G = {候选 m 的 Doppler 邻域 ±radius}

与 A4 简化版的不同: 此处用簇内联合残差投影计算组得分，
无需矩阵求逆 (直接用邻域原子叠加权向量)，数值稳定。
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from channel import (
    compute_chirp_params, sample_clustered_geometry,
    bdlr_gains, build_path_matrix, compute_effective_channel,
)
from grid_stats import build_candidate_atoms
from detection_mc import _estimate_true_candidate_mask
from full_frame import build_frame_and_transmit


# ======================================================================
# 组原子构造
# ======================================================================

def build_group_atoms(atoms, grid, max_doppler, dirichlet_r):
    """
    对每个候选格点 m, 构造其 Doppler 邻域组。

    组原子 = 邻域内所有候选原子的复数加权叠加。
    权重 = Dirichlet 核系数 (与距离相关)。

    返回:
      group_atoms: (M, n_sc) 每组一条复合原子向量
      group_masks: list[ndarray] 每组包含的候选索引
    """
    M = grid.shape[0]
    n_sc = atoms.shape[1]
    radius = dirichlet_r

    # 预计算权重: w_k = exp(-|k| / radius) 归一化
    half = 2 * radius + 1
    w_raw = np.exp(-np.abs(np.arange(-radius, radius + 1)) / max(radius, 1))
    w_raw /= w_raw.sum()

    group_atoms = np.zeros((M, n_sc), dtype=complex)
    group_masks = []

    for m in range(M):
        dm = int(grid[m, 0])
        nu_m = int(grid[m, 1])
        members = []
        for di in range(-radius, radius + 1):
            dk = nu_m + di
            if abs(dk) <= max_doppler:
                # 查找网格中匹配 (dm, dk) 的索引
                match = np.where((grid[:, 0] == dm) & (grid[:, 1] == dk))[0]
                if len(match) > 0:
                    members.append((match[0], w_raw[di + radius]))
        if members:
            for idx, w in members:
                group_atoms[m] += w * atoms[idx]
            nrm = np.linalg.norm(group_atoms[m])
            if nrm > 1e-12:
                group_atoms[m] /= nrm
        else:
            group_atoms[m] = atoms[m].copy()
        group_masks.append(np.array([mi for mi, _ in members]))

    return group_atoms, group_masks


# ======================================================================
# 组级 OMP
# ======================================================================

def group_omp_detect(y, atoms, grid, group_atoms, group_masks,
                     max_paths, kappa, residual_power_est=None):
    """
    组级 OMP 检测。

    每轮: 取 group_atoms 与残差的最大内积 → 选取整组 → 联合移除。

    参数:
      y:                  (N,) 观测
      atoms:              (M, N) 单候选原子
      grid:               (M, 2)
      group_atoms:        (M, N) 组原子
      group_masks:        list 每组包含的单候选索引
      max_paths:          int 最大输出组数
      kappa:              float 门限乘数
      residual_power_est: float 残差功率估计（可选, 默认从 y 估计）

    返回:
      detected_groups:  list[dict] 检测到的组 (center_idx, members, gain_vec)
      detected_singles: list[int] 检测到的单候选索引 (去重后)
      n_iter:           int 迭代次数
    """
    M = atoms.shape[0]
    n_sc = len(y)
    r = y.copy()
    detected_groups = []
    detected_singles_set = set()
    detected_gains = {}  # grid_idx → matched-filter gain

    if residual_power_est is None:
        residual_power_est = float(np.mean(np.abs(y) ** 2))

    threshold = kappa * residual_power_est

    for _ in range(max_paths):
        # 组得分
        scores = np.abs(group_atoms @ r.conj()) ** 2
        best_m = int(np.argmax(scores))
        best_score = float(scores[best_m])

        if best_score < threshold:
            break

        # 组内成员
        members = list(group_masks[best_m])

        # 联合增益估计: 用组内所有原子对残差做最小二乘
        A_group = atoms[members, :].T  # (N, K)
        g_vec, _, _, _ = np.linalg.lstsq(A_group, r, rcond=None)

        # 记录增益
        for mi, gv in zip(members, g_vec):
            if mi in detected_gains:
                detected_gains[mi] += gv
            else:
                detected_gains[mi] = gv

        # 移除整组
        r -= A_group @ g_vec

        detected_groups.append({
            "center_idx": best_m,
            "center_delay": int(grid[best_m, 0]),
            "center_doppler": int(grid[best_m, 1]),
            "members": members,
            "gains": g_vec,
            "score": best_score,
        })
        for mi in members:
            detected_singles_set.add(int(mi))

    return detected_groups, list(detected_singles_set), len(detected_groups), detected_gains


def single_omp_detect(y, atoms, max_paths, kappa, residual_power_est=None):
    """
    标准单径 OMP (对照组)。

    返回:
      detected_indices: list[int]
    """
    M = atoms.shape[0]
    r = y.copy()
    detected = []
    detected_gains = {}

    if residual_power_est is None:
        residual_power_est = float(np.mean(np.abs(y) ** 2))
    threshold = kappa * residual_power_est

    for _ in range(max_paths):
        scores = np.abs(atoms @ r.conj()) ** 2
        best_m = int(np.argmax(scores))
        if float(scores[best_m]) < threshold:
            break
        atom_b = atoms[best_m, :]
        g_est = np.dot(atom_b.conj(), r)
        r -= g_est * atom_b
        detected.append(best_m)
        detected_gains[best_m] = g_est

    return detected, detected_gains


# ======================================================================
# 对比实验
# ======================================================================

def run_comparison(cfg):
    """
    单次 trial: 对比单径 OMP vs 组级 OMP 的 Pd + 假径数。
    """
    n_sc         = int(cfg["n_sc"])
    n_paths      = int(cfg["n_paths"])
    max_delay    = int(cfg["max_delay"])
    max_doppler  = int(cfg["max_doppler"])
    doppler_guard = int(cfg["doppler_guard"])
    dirichlet_r  = int(cfg["dirichlet_r"])
    p_fa         = float(cfg["p_fa"])
    pilot_idx    = int(cfg.get("pilot_idx", 0))
    seed         = int(cfg.get("seed", 42))
    use_calibrated_kappa = cfg.get("use_calibrated_kappa", True)

    rng = np.random.default_rng(seed)
    y, x_tx, sym, H, delays, dopplers, gains, N0, c1 = \
        build_frame_and_transmit(cfg, rng)

    atoms, grid = build_candidate_atoms(
        n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx)
    group_atoms, group_masks = build_group_atoms(atoms, grid, max_doppler, dirichlet_r)

    M = grid.shape[0]
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    true_set = set(np.where(true_mask)[0])

    # 残差功率: 使用噪声底，不是全帧能量
    # noise_power = 1.0 (build_frame_and_transmit 中固定)
    res_pow = 1.0  # σ² = N0

    # κ 选择
    if use_calibrated_kappa:
        kappa_single = 24.6 if dirichlet_r == 3 else 31.0
        kappa_group  = kappa_single * 0.85  # 组原子 SNR 增益: 门限可略降
    else:
        kappa_single = np.log(M / p_fa)
        kappa_group  = np.log(M / p_fa)

    # 单径 OMP
    detected_s = single_omp_detect(y, atoms, n_paths + 2, kappa_single, res_pow)
    true_s  = len(set(detected_s) & true_set)
    false_s = len(set(detected_s) - true_set)

    # 组级 OMP
    groups, detected_g, n_iter = group_omp_detect(
        y, atoms, grid, group_atoms, group_masks,
        n_paths + 2, kappa_group, res_pow)
    true_g  = len(set(detected_g) & true_set)
    false_g = len(set(detected_g) - true_set)

    return {
        "pd_single": true_s / max(n_paths, 1),
        "pd_group":  true_g / max(n_paths, 1),
        "false_single": false_s,
        "false_group": false_g,
        "n_detected_single": len(detected_s),
        "n_detected_group": len(detected_g),
        "n_groups": n_iter,
        "pd_diff": (true_g - true_s) / max(n_paths, 1),
    }


def run_batch(cfg, n_trials=20):
    """批量对比。"""
    rows = []
    print(f"{'trial':>5} {'Pd_s':>6} {'Pd_g':>6} {'dPd':>6} "
          f"{'F_s':>5} {'F_g':>5} {'n_s':>4} {'n_g':>4} {'grp':>4}")
    print("-" * 55)
    for t in range(n_trials):
        _c = dict(cfg)
        _c["seed"] = cfg.get("seed", 42) + t
        r = run_comparison(_c)
        rows.append(r)
        print(f"{t:5d} {r['pd_single']:6.3f} {r['pd_group']:6.3f} "
              f"{r['pd_diff']:+6.3f} {r['false_single']:5d} {r['false_group']:5d} "
              f"{r['n_detected_single']:4d} {r['n_detected_group']:4d} "
              f"{r['n_groups']:4d}")

    # 汇总
    n = len(rows)
    print(f"\n=== 汇总 (n={n}) ===")
    for key in ["pd_single", "pd_group", "false_single", "false_group", "pd_diff"]:
        vals = [r[key] for r in rows]
        print(f"  {key:20s}: mean={np.mean(vals):.4f}  std={np.std(vals):.4f}")
    print(f"  Pd 提升: {np.mean([r['pd_diff'] for r in rows]):+.4f} "
          f"({np.mean([r['pd_diff'] for r in rows])*100:+.1f} pp)")
    print(f"  假径减少: {np.mean([r['false_single'] for r in rows]):.1f} → "
          f"{np.mean([r['false_group'] for r in rows]):.1f}")

    return rows


if __name__ == "__main__":
    base = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "rho": 0.5, "p_fa": 1e-3, "seed": 42, "pilot_idx": 0,
        "use_calibrated_kappa": True,
    }

    # --- SNR 扫描 (固定 κ=24.6) ---
    print("=== SNR 扫描: 单径 vs 组级 OMP ===\n")
    for snr in [0, 5, 10, 15]:
        _c = dict(base)
        _c["snr_db"] = snr
        _c["seed"] = 100
        rows = run_batch(_c, n_trials=20)
        pd_s = np.mean([r["pd_single"] for r in rows])
        pd_g = np.mean([r["pd_group"] for r in rows])
        f_s  = np.mean([r["false_single"] for r in rows])
        f_g  = np.mean([r["false_group"] for r in rows])
        print(f"  SNR={snr:2d}dB: Pd_s={pd_s:.3f} Pd_g={pd_g:.3f} "
              f"(dPd={pd_g-pd_s:+.3f})  F_s={f_s:.1f} F_g={f_g:.1f}\n")

    # --- κ 扫描 (固定 SNR=10dB) ---
    print("=== κ 扫描: SNR=10dB ===\n")
    print(f"{'kappa':>8} {'Pd_s':>7} {'Pd_g':>7} {'dPd':>7} {'F_s':>6} {'F_g':>6}")
    print("-" * 48)
    for kappa in [10, 15, 20, 25, 30, 40, 60]:
        _c = dict(base)
        _c["snr_db"] = 10
        _c["use_calibrated_kappa"] = False
        _c["kappa_override"] = kappa
        _c["seed"] = 200
        # 临时覆盖
        import full_frame

        # 手动运行带自定义 kappa 的检测
        rows_list = []
        for t in range(15):
            _c2 = dict(_c)
            _c2["seed"] = _c["seed"] + t
            rng = np.random.default_rng(_c2["seed"])
            y, x_tx, sym, H, delays, dopplers, gains, N0, c1 = \
                full_frame.build_frame_and_transmit(_c2, rng)
            n_sc = _c2["n_sc"]
            atoms, grid = build_candidate_atoms(
                n_sc, c1, _c2["dirichlet_r"], _c2["max_delay"],
                _c2["max_doppler"], _c2["pilot_idx"])
            group_atoms, group_masks = build_group_atoms(
                atoms, grid, _c2["max_doppler"], _c2["dirichlet_r"])
            true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
            true_set = set(np.where(true_mask)[0])

            res_pow = 1.0
            n_paths = _c2["n_paths"]
            k_s = kappa
            k_g = kappa * 0.85

            det_s = single_omp_detect(y, atoms, n_paths + 2, k_s, res_pow)
            groups, det_g, n_iter = group_omp_detect(
                y, atoms, grid, group_atoms, group_masks,
                n_paths + 2, k_g, res_pow)

            rows_list.append({
                "pd_single": len(set(det_s) & true_set) / n_paths,
                "pd_group": len(set(det_g) & true_set) / n_paths,
                "false_single": len(set(det_s) - true_set),
                "false_group": len(set(det_g) - true_set),
            })
        pd_s = np.mean([r["pd_single"] for r in rows_list])
        pd_g = np.mean([r["pd_group"] for r in rows_list])
        f_s  = np.mean([r["false_single"] for r in rows_list])
        f_g  = np.mean([r["false_group"] for r in rows_list])
        print(f"{kappa:8.0f} {pd_s:7.3f} {pd_g:7.3f} {pd_g-pd_s:+7.3f} "
              f"{f_s:6.1f} {f_g:6.1f}")
