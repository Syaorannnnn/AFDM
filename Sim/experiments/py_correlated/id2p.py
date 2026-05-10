"""
ID2P 频谱分解：将数据子载波在导频位置的泄漏分解为 i.i.d. 基线 B(j)
和 ρ-调制项 Delta(j)。

核心公式（来自四篇精读总结 §建模层）:

    对第 c 个簇, 数据子载波 j:
      A_c(j) = sum_{i∈c} Phi_i(pilot, data_j)           (相干叠加)
      B_c(j) = sum_{i∈c} |Phi_i(pilot, data_j)|²        (非相干叠加)
      Delta_c(j) = |A_c(j)|² - B_c(j)                     (ρ-调制项)

    Delta(j) < 0 → 该子载波为 ID2P 的"天然零陷"(null)
    Delta(j) > 0 → 该子载波为 ID2P 的"热点"(hotspot)

    增益相关时: E[|H(pilot, j)|²] = gain_var × (B_total + ρ × Delta_total)

参考:
  id2p_correlation_probe.py 中的 compute_ab_delta() 和 mc_id2p_vector()
"""

import numpy as np


def compute_id2p_spectrum(Phi_list, cluster_labels, n_data, pilot_idx, data_idx):
    """
    逐簇分解 ID2P 频谱。

    参数:
      Phi_list:       length P list of (N,N) 每径的 Dirichlet 展宽矩阵
      cluster_labels: (P,) int 每径归属的簇编号 0..C-1
      n_data:         int 数据子载波数
      pilot_idx:      int 导频在 DAFT 域的索引
      data_idx:       (n_data,) int 数据子载波索引

    返回:
      B_by_cluster: (C, n_data) float 每簇每子载波的 B_c(j)
      D_by_cluster: (C, n_data) float 每簇每子载波的 Delta_c(j)
      B_total:      (n_data,) float    sum_c B_c(j)
      D_total:      (n_data,) float    sum_c Delta_c(j)
    """
    n_clusters = cluster_labels.max() + 1
    A_by_cluster = np.zeros((n_clusters, n_data), dtype=complex)
    B_by_cluster = np.zeros((n_clusters, n_data))

    for i, Phi_i in enumerate(Phi_list):
        c = cluster_labels[i]
        row = Phi_i[pilot_idx, :][data_idx]  # 导频位置到各数据子载波的耦合
        A_by_cluster[c] += row
        B_by_cluster[c] += np.abs(row) ** 2

    D_by_cluster = np.abs(A_by_cluster) ** 2 - B_by_cluster
    B_total = B_by_cluster.sum(axis=0)
    D_total = D_by_cluster.sum(axis=0)
    return B_by_cluster, D_by_cluster, B_total, D_total


def compute_effective_id2p(B_total, D_total, rho, gain_var):
    """
    给定相关度 rho 下的期望 ID2P 功率。

    E[|H(0, j)|²] = gain_var × (B_total(j) + rho × D_total(j))

    验证: rho=0 → gain_var × B_total (i.i.d. 基线)
          rho=1 → gain_var × (B_total + D_total) (全相关)
    """
    return gain_var * (B_total + rho * D_total)


def compute_null_hotspot_map(D_total):
    """
    分类每个数据子载波为 null 或 hotspot。

    返回:
      is_null:    (n_data,) bool
      is_hotspot: (n_data,) bool
      n_null:     int
      n_hotspot:  int
    """
    is_null = D_total < 0
    is_hotspot = D_total > 0
    return is_null, is_hotspot, int(is_null.sum()), int(is_hotspot.sum())


def compute_spatial_cov(
    delays, dopplers, cluster_labels, n_sc, c1, radius,
    pilot_idx, data_idx, n_mc, rho, n_clusters, rng=None
):
    """
    Monte Carlo 估计 ID2P 检测统计量的空间协方差矩阵。

    对每对数据子载波 (j, j'), 估计 Cov(|H(pilot, j)|², |H(pilot, j')|²)。

    参数:
      delays, dopplers: 路径几何 (P,)
      cluster_labels:   (P,) int
      n_sc:             int 子载波数
      c1:               float chirp 参数
      radius:           int Dirichlet 展宽半径
      pilot_idx:        int
      data_idx:         (n_data,)
      n_mc:             int Monte Carlo 次数
      rho:              float BDLR ρ
      n_clusters:       int
      rng:              numpy Generator

    返回:
      cov: (n_data, n_data) float64 协方差矩阵
    """
    try:
        from .channel import build_path_matrix, bdlr_gains
    except ImportError:
        from channel import build_path_matrix, bdlr_gains  # 交互运行后备

    if rng is None:
        rng = np.random.default_rng()

    n_data = len(data_idx)
    n_paths = len(delays)

    # 预先构建所有路径的固定-几何 Phi 矩阵
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    Phi_fixed = [
        build_path_matrix(int(delays[i]), dopplers[i], n_sc, loc_step, radius)
        for i in range(n_paths)
    ]

    power_samples = np.zeros((n_mc, n_data))
    for t in range(n_mc):
        gains = bdlr_gains(n_paths, n_clusters, rho, rng)
        row = np.zeros(n_data, dtype=complex)
        for i in range(n_paths):
            row += gains[i] * Phi_fixed[i][pilot_idx, :][data_idx]
        power_samples[t, :] = np.abs(row) ** 2

    # 无偏协方差
    cov = np.cov(power_samples, rowvar=False)
    return cov


def decompose_by_cluster(B_cluster, D_cluster, cluster_names=None):
    """
    将逐簇 B 和 Delta 汇总为可读字典。

    返回:
      summary: dict, 含 'cluster_0_B', 'cluster_0_D', ...
    """
    n_clusters = B_cluster.shape[0]
    if cluster_names is None:
        cluster_names = [f'cluster_{c}' for c in range(n_clusters)]
    out = {}
    for c in range(n_clusters):
        out[f'{cluster_names[c]}_B'] = B_cluster[c]
        out[f'{cluster_names[c]}_D'] = D_cluster[c]
    return out
