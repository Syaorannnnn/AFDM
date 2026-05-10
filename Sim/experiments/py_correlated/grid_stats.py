"""
延迟-多普勒候选网格上的 pilot 匹配统计量。

该模块用于区分两类有效候选数：
  1. ID2P 数据子载波维度的 M_eff（见 id2p.py）
  2. OMP/CFAR 搜索网格维度的 M_eff（见本模块）

本模块只使用 Python 仿真，不调用 MATLAB 主干脚本。
"""

import numpy as np

from channel import (
    bdlr_gains,
    build_path_matrix,
    compute_chirp_params,
    sample_clustered_geometry,
)


def build_candidate_grid(max_delay, max_doppler):
    """
    构造延迟-整数 Doppler 候选网格。

    M = (max_delay + 1) × (2 × max_doppler + 1)

    参数:
      max_delay:   int 最大延迟索引
      max_doppler: int 最大 Doppler 索引

    返回:
      grid: (M, 2) int，每行 [delay, doppler_int]
    """
    rows = []
    for delay in range(max_delay + 1):
        for doppler_int in range(-max_doppler, max_doppler + 1):
            rows.append([delay, doppler_int])
    return np.asarray(rows, dtype=int)


def build_candidate_atoms(n_sc, c1, dirichlet_r, max_delay, max_doppler,
                          pilot_idx=0):
    """
    构造 pilot 行上的归一化候选原子。

    每个候选网格点 (delay, doppler_int) 对应一个长度为 N 的复向量，
    即对应 Dirichlet 核展宽后的 pilot 行响应（归一化到单位范数）。

    参数:
      n_sc:         int 子载波数
      c1:           float chirp 参数
      dirichlet_r:  int Dirichlet 展宽半径
      max_delay:    int
      max_doppler:  int
      pilot_idx:    int 导频位置

    返回:
      atoms: (M, N) complex128 每行是单位范数的候选原子
      grid:  (M, 2) int 对应的候选网格
    """
    c1_eff, loc_step = compute_chirp_params(n_sc, max_doppler, 0)
    # 始终使用传入的 c1 计算 loc_step
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    grid = build_candidate_grid(max_delay, max_doppler)
    M = grid.shape[0]
    atoms = np.zeros((M, n_sc), dtype=complex)
    for idx, (delay, doppler_int) in enumerate(grid):
        phi = build_path_matrix(
            int(delay), float(doppler_int), n_sc, loc_step, dirichlet_r
        )
        row = phi[pilot_idx, :]
        norm_val = np.linalg.norm(row)
        if norm_val > 1e-12:
            row = row / norm_val
        atoms[idx, :] = row
    return atoms, grid


def estimate_grid_stat_covariance(cfg):
    """
    Monte Carlo 估计候选网格匹配统计量的协方差。

    流程:
      1. 采样簇化几何
      2. 构建归一化候选原子矩阵
      3. 对每次 MC: 生成 BDLR 增益 → 计算 pilot 响应 → 取原子内积 → 功率
      4. 从 MC 功率样本计算协方差

    参数:
      cfg: dict，必须包含:
        n_sc, n_paths, n_clusters, max_delay, max_doppler,
        doppler_guard, dirichlet_r, rho, mc_samples, seed

    返回:
      dict:
        cov:           (M, M) 协方差矩阵
        power_samples: (mc_samples, M) 原始功率样本
        grid:          (M, 2) 候选网格
        delays:        (n_paths,) 采样延迟
        dopplers:      (n_paths,) 采样 Doppler
    """
    n_sc         = int(cfg.get("n_sc", 64))
    n_paths      = int(cfg.get("n_paths", 6))
    n_clusters   = int(cfg.get("n_clusters", 2))
    max_delay    = int(cfg.get("max_delay", 3))
    max_doppler  = int(cfg.get("max_doppler", 2))
    doppler_guard = int(cfg.get("doppler_guard", 3))
    dirichlet_r  = int(cfg.get("dirichlet_r", 3))
    rho          = float(cfg.get("rho", 0.5))
    mc_samples   = int(cfg.get("mc_samples", 200))
    seed         = int(cfg.get("seed", 42))
    pilot_idx    = int(cfg.get("pilot_idx", 0))

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng
    )
    atoms, grid = build_candidate_atoms(
        n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx
    )

    path_mats = [
        build_path_matrix(int(delays[idx]), dopplers[idx], n_sc, loc_step,
                          dirichlet_r)
        for idx in range(n_paths)
    ]
    M_cand = grid.shape[0]
    power_samples = np.zeros((mc_samples, M_cand))
    for sample_idx in range(mc_samples):
        gains = bdlr_gains(n_paths, n_clusters, rho, rng)
        pilot_row = np.zeros(n_sc, dtype=complex)
        for path_idx, phi in enumerate(path_mats):
            pilot_row += gains[path_idx] * phi[pilot_idx, :]
        matched = atoms.conj() @ pilot_row
        power_samples[sample_idx, :] = np.abs(matched) ** 2

    cov = np.cov(power_samples, rowvar=False)
    return {
        "cov": cov,
        "power_samples": power_samples,
        "grid": grid,
        "delays": delays,
        "dopplers": dopplers,
    }
