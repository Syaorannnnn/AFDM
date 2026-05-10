"""
组证据统计量 ΔL_G：完整实现（含组内相关矩阵 R_c）。

来自四篇精读总结 §观测层:

    ΔL_G = y^H C0^{-1} Phi_G (R_G^{-1} + Phi_G^H C0^{-1} Phi_G)^{-1} Phi_G^H C0^{-1} y
           - log|I + R_G Phi_G^H C0^{-1} Phi_G|

单径 ID2P 是 |G|=1 且 R_G=sigma² 的特例。

本模块对比简化版（簇中心原子）和完整版（R_c 矩阵）的分离性能。
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from channel import (
    compute_chirp_params, sample_clustered_geometry,
    bdlr_gains, build_path_matrix,
)
from grid_stats import build_candidate_atoms, build_candidate_grid
from detection_mc import _estimate_true_candidate_mask
from cfar import compute_effective_candidates


def compute_full_group_evidence(
    delays, dopplers, cluster_labels, n_sc, c1, dirichlet_r,
    pilot_idx, n_clusters, rho, gain_var, rng,
    max_delay, max_doppler
):
    """
    完整版 ΔL_G：对每个簇构建 Phi_G 和 R_G，计算边际似然增益。

    参数:
      与 simplified 版本相同。

    返回:
      group_evidence:  (n_clusters, M_cand)  每个簇在每个候选格点的 ΔL_G
      single_evidence: (n_paths, M_cand)     单径匹配功率（归一化）
      group_simple:    (n_clusters, M_cand)  简化版（簇中心原子）
    """
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    atoms, grid = build_candidate_atoms(n_sc, c1, dirichlet_r,
                                        max_delay, max_doppler, pilot_idx)
    n_paths = len(delays)
    pc = n_paths // n_clusters
    M_cand = grid.shape[0]

    # --- 真实 pilot 响应 ---
    path_mats = [
        build_path_matrix(int(delays[i]), dopplers[i], n_sc, loc_step, dirichlet_r)
        for i in range(n_paths)
    ]
    gains = bdlr_gains(n_paths, n_clusters, rho, rng)
    pilot_row = np.zeros(n_sc, dtype=complex)
    for i, phi in enumerate(path_mats):
        pilot_row += gains[i] * phi[pilot_idx, :]

    # --- 残差协方差 C0 = sigma² I ---
    # 用 pilot_row 非导频位置的行范数估计
    non_pilot_mask = np.ones(n_sc, dtype=bool)
    non_pilot_mask[pilot_idx] = False
    sigma2 = float(np.mean(np.abs(pilot_row[non_pilot_mask]) ** 2)) + 1e-12
    C0_inv = np.eye(n_sc) / sigma2

    # --- 单径证据 ---
    single_evidence = np.abs(atoms @ pilot_row.conj()) ** 2 / sigma2

    # --- 组证据（简化版：簇中心原子）---
    group_simple = np.zeros((n_clusters, M_cand))
    for c in range(n_clusters):
        lo = c * pc
        hi = lo + pc
        c_delays = delays[lo:hi]
        c_dopp = dopplers[lo:hi]
        base_d = int(np.median(c_delays))
        base_nu = float(np.median(c_dopp))
        phi_c = build_path_matrix(base_d, base_nu, n_sc, loc_step, dirichlet_r)
        c_atom = phi_c[pilot_idx, :]
        norm_val = np.linalg.norm(c_atom)
        if norm_val > 1e-12:
            c_atom = c_atom / norm_val
        group_simple[c, :] = np.abs(atoms @ c_atom.conj()) ** 2 / sigma2

    # --- 完整 ΔL_G: 候选邻域组 ---
    # 对每个候选格点 m=(delay_d, doppler_k)，构造以 m 为中心的 Doppler 邻域组。
    # 邻域半径 = dirichlet_r，模拟分数 Doppler 展宽导致的相邻候选相关性。
    # Phi_G_m 的每一列是邻域内一个格点的归一化原子。
    # R_G 编码邻域内原子之间的 Dirichlet 相关结构。

    group_radius = dirichlet_r
    group_size = 2 * group_radius + 1  # 组大小
    group_evidence_full = np.zeros((1, M_cand))  # 所有候选共享同一组结构

    # R_G: Dirichlet 核相关矩阵 —— 相邻整数 Doppler 候选的分数泄漏相关
    # 简化: 用指数衰减 |i-j|/radius
    R_G = np.zeros((group_size, group_size))
    for i in range(group_size):
        for j in range(group_size):
            R_G[i, j] = np.exp(-2.0 * abs(i - j) / max(group_radius, 1))
    # 添加单位矩阵确保满秩
    R_G = 0.9 * R_G + 0.1 * np.eye(group_size)

    for m in range(M_cand):
        delay_m = int(grid[m, 0])
        doppler_m = int(grid[m, 1])

        # 邻域组: 固定 delay, 相邻 Doppler
        Phi_G_m = np.zeros((n_sc, group_size), dtype=complex)
        valid_cols = 0
        col_map = []
        for di in range(-group_radius, group_radius + 1):
            dk = doppler_m + di
            # 检查是否在候选网格内
            if abs(dk) <= max_doppler:
                tmp_phi = build_path_matrix(
                    delay_m, float(dk), n_sc, loc_step, dirichlet_r
                )
                col_vec = tmp_phi[pilot_idx, :]
                nrm = np.linalg.norm(col_vec)
                if nrm > 1e-12:
                    col_vec = col_vec / nrm
                Phi_G_m[:, valid_cols] = col_vec
                col_map.append(valid_cols)
                valid_cols += 1

        if valid_cols < 2:
            group_evidence_full[0, m] = single_evidence[m]
            continue

        # 裁剪 Phi_G 和 R_G
        Phi_G_use = Phi_G_m[:, :valid_cols]
        R_G_use = R_G[np.ix_(col_map, col_map)]

        # Gram
        Phi_b = Phi_G_use.conj().T @ C0_inv  # (Pc, N)
        G_m = Phi_b @ Phi_G_use               # (Pc, Pc)

        R_inv = np.linalg.inv(R_G_use)
        M_mat = R_inv + G_m

        try:
            M_inv = np.linalg.inv(M_mat)
        except np.linalg.LinAlgError:
            M_inv = np.linalg.pinv(M_mat)

        # 得分项
        yb = pilot_row.conj() @ C0_inv
        v = yb @ Phi_G_use
        score = float(np.real(v @ M_inv @ v.conj()))

        # 惩罚项
        _, logdet = np.linalg.slogdet(np.eye(valid_cols) + R_G_use @ G_m)
        penalty = float(logdet)

        group_evidence_full[0, m] = max(score - penalty, 0.0)

    return group_evidence_full, single_evidence, group_simple, atoms, grid, gains, pilot_row


def compare_evidence_methods(cfg):
    """
    对比单径、简化组、完整 ΔL_G 三种方法在真实候选上的分离性能。

    返回:
      dict 含均值、最大化候选索引、误匹配率。
    """
    n_sc = int(cfg.get("n_sc", 64))
    n_paths = int(cfg.get("n_paths", 6))
    n_clusters = int(cfg.get("n_clusters", 2))
    max_delay = int(cfg.get("max_delay", 3))
    max_doppler = int(cfg.get("max_doppler", 2))
    doppler_guard = int(cfg.get("doppler_guard", 3))
    dirichlet_r = int(cfg.get("dirichlet_r", 3))
    rho = float(cfg.get("rho", 0.5))
    seed = int(cfg.get("seed", 42))
    pilot_idx = int(cfg.get("pilot_idx", 0))

    rng = np.random.default_rng(seed)
    c1, _ = compute_chirp_params(n_sc, max_doppler, doppler_guard)
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    cluster_labels = np.repeat(np.arange(n_clusters), n_paths // n_clusters)
    gain_var = 1.0 / n_paths

    g_full, s_ev, g_simple, atoms, grid, gains, pilot_row = \
        compute_full_group_evidence(
            delays, dopplers, cluster_labels, n_sc, c1, dirichlet_r,
            pilot_idx, n_clusters, rho, gain_var, rng,
            max_delay, max_doppler)

    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    false_mask = ~true_mask

    # 每种方法在 true/false 候选上的得分
    s_true  = s_ev[true_mask]
    s_false = s_ev[false_mask]
    gs_true  = g_simple.max(axis=0)[true_mask]
    gs_false = g_simple.max(axis=0)[false_mask]
    gf_true  = g_full.max(axis=0)[true_mask]
    gf_false = g_full.max(axis=0)[false_mask]

    # 最强候选匹配
    best_single = int(np.argmax(s_ev))
    best_group_s = int(np.argmax(g_simple.max(axis=0)))
    best_group_f = int(np.argmax(g_full.max(axis=0)))

    true_indices = np.where(true_mask)[0]
    hit_single = int(best_single in true_indices)
    hit_group_s = int(best_group_s in true_indices)
    hit_group_f = int(best_group_f in true_indices)

    return {
        "rho": rho,
        "dir_r": dirichlet_r,
        # 均值
        "s_mean_true": float(np.mean(s_true)),
        "s_mean_false": float(np.mean(s_false)),
        "gs_mean_true": float(np.mean(gs_true)),
        "gs_mean_false": float(np.mean(gs_false)),
        "gf_mean_true": float(np.mean(gf_true)),
        "gf_mean_false": float(np.mean(gf_false)),
        # 最大/最小分离比
        "sep_single": float(np.mean(s_true) / max(np.mean(s_false), 1e-15)),
        "sep_group_simple": float(np.mean(gs_true) / max(np.mean(gs_false), 1e-15)),
        "sep_group_full": float(np.mean(gf_true) / max(np.mean(gf_false), 1e-15)),
        # top-1 命中
        "hit_single": hit_single,
        "hit_group_simple": hit_group_s,
        "hit_group_full": hit_group_f,
        # 有效候选
        "n_true": int(true_mask.sum()),
        "n_false": int(false_mask.sum()),
        # 完整证据的对比度
        "gf_contrast": float(
            np.max(g_full[:, true_mask]) / max(np.max(g_full[:, false_mask]), 1e-15)
        ),
    }


def run_evidence_scan(base_cfg):
    """扫描 rho × dir_r, 输出对比表。"""
    print(f"{'rho':>5} {'dir_r':>5} {'sep_s':>8} {'sep_gs':>8} {'sep_gf':>8} "
          f"{'hit_s':>6} {'hit_gs':>7} {'hit_gf':>7} {'gf_contrast':>11}")
    print("-" * 70)
    rows = []
    for dr in [1, 3]:
        for rho_v in [0.0, 0.5, 0.9]:
            _c = dict(base_cfg)
            _c["dirichlet_r"] = dr
            _c["rho"] = rho_v
            r = compare_evidence_methods(_c)
            rows.append(r)
            print(f"{r['rho']:5.1f} {r['dir_r']:5d} "
                  f"{r['sep_single']:8.1f} {r['sep_group_simple']:8.1f} {r['sep_group_full']:8.1f} "
                  f"{r['hit_single']:6d} {r['hit_group_simple']:7d} {r['hit_group_full']:7d} "
                  f"{r['gf_contrast']:11.2f}")
    return rows


if __name__ == '__main__':
    base = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "rho": 0.5, "seed": 42, "pilot_idx": 0,
    }
    print("=== A4: Full ΔL_G Group Evidence Comparison ===\n")
    run_evidence_scan(base)
