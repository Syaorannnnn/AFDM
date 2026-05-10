"""
BDLR 相关信道模型 + AFDM Dirichlet 核有效信道构造。

BDLR 模型将 P 条路径划分为 C 个簇（cluster），每簇 Pc 条子径。
簇内增益共享一个公共复高斯分量，叠加独立扰动：

    h_{c,i} = common_c + independent_{c,i}
    common_c        ~ CN(0, gain_var × ρ)
    independent_{c,i} ~ CN(0, gain_var × (1-ρ))

ρ=0 → 完全 i.i.d.（退化为 GiFreeChannelSampler 行为）
ρ=1 → 同一簇内所有子径增益完全相同。

参考文献:
  Medeiros et al., "CA-CFAR Detection for SAR Systems Over Correlated
  Gamma-Distributed Clutter", IEEE GRSL, 2024.  (相关 Gamma 模型)
  Sheikh & Budhiraja, "Correlated Block-Sparse Channel Estimation ...",
  IEEE TWC, 2025.  (BDLR block-sparse 先验)
"""

import numpy as np


# ======================================================================
# Chirp 参数
# ======================================================================

def compute_chirp_params(n_sc, max_doppler, doppler_guard):
    """
    AFDM chirp 参数 c1, loc_step。

    与 MATLAB GiFreeConfig 等价:
      c1 = (2×(k_max + guard) + 1) / (2×N)
      loc_step = (2×N×c1) mod N

    参数:
      n_sc:         子载波数 N
      max_doppler:  最大 Doppler 索引 k_max
      doppler_guard: Doppler 保护间隔

    返回:
      c1:       chirp 参数 1
      loc_step: DAFT 域每延迟步进对应的位置偏移
    """
    c1 = (2.0 * (max_doppler + doppler_guard) + 1.0) / (2.0 * n_sc)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    return c1, loc_step


# ======================================================================
# 路径几何采样
# ======================================================================

def sample_path_geometry(n_paths, max_delay, max_doppler, rng=None):
    """
    采样延迟-多普勒几何位置。

    延迟从中均匀无放回抽取；多普勒按 cos(θ) 分布。
    避免 (delay, round(doppler)) 冲突。

    参数:
      n_paths:     路径数 P
      max_delay:   最大延迟索引 l_max
      max_doppler: 最大 Doppler 索引 k_max
      rng:         numpy Generator（可选）

    返回:
      delays:   (n_paths,) int 延迟索引
      dopplers: (n_paths,) float 分数 Doppler 值
    """
    if rng is None:
        rng = np.random.default_rng()

    all_delays = np.arange(max_delay + 1)
    n_avail = len(all_delays)
    repeat = int(np.ceil(n_paths / n_avail))
    pool = np.tile(all_delays, repeat)
    rng.shuffle(pool)
    delays = pool[:n_paths].copy()

    dopplers = np.zeros(n_paths)
    for i in range(n_paths):
        for _ in range(100):
            phase = rng.uniform(-np.pi, np.pi)
            dop = max_doppler * np.cos(phase)
            conflict = any(
                delays[i] == delays[j] and round(dop) == round(dopplers[j])
                for j in range(i)
            )
            if not conflict:
                dopplers[i] = dop
                break
        else:
            dopplers[i] = dop  # 100 次未解决则接受
    return delays, dopplers


# ======================================================================
# BDLR 增益
# ======================================================================

def bdlr_gains(n_paths, n_clusters, rho, rng=None):
    """
    BDLR 簇内相关复增益。

    每簇: common ~ CN(0, var×ρ), independent ~ CN(0, var×(1-ρ))。

    参数:
      n_paths:    路径总数 P
      n_clusters: 簇数 C（P 必须能被 C 整除）
      rho:        簇内相关系数 ∈ [0, 1]
      rng:        numpy Generator（可选）

    返回:
      gains: (n_paths,) complex128 复路径增益
    """
    if n_paths % n_clusters != 0:
        raise ValueError(
            f"n_paths={n_paths} 不能被 n_clusters={n_clusters} 整除"
        )
    if not 0.0 <= rho <= 1.0:
        raise ValueError(f"rho={rho} 超出 [0, 1]")

    if rng is None:
        rng = np.random.default_rng()

    pc = n_paths // n_clusters
    gain_var = 1.0 / n_paths  # 归一化总功率 = 1

    h = np.zeros(n_paths, dtype=complex)
    for c in range(n_clusters):
        lo = c * pc
        hi = lo + pc

        common = (
            np.sqrt(gain_var * rho / 2.0)
            * (rng.standard_normal() + 1j * rng.standard_normal())
        )
        independent = (
            np.sqrt(gain_var * (1.0 - rho) / 2.0)
            * (rng.standard_normal(pc) + 1j * rng.standard_normal(pc))
        )
        h[lo:hi] = common + independent
    return h


# ======================================================================
# 簇化几何采样（路径在 DAFT 域有可观测重叠）
# ======================================================================

def sample_clustered_geometry(n_paths, n_clusters, max_delay, max_doppler,
                               rng=None):
    """
    采样簇化延迟-多普勒几何，使同簇内路径在 DAFT 域有重叠。

    延迟: 每簇分配一个基准延迟，簇内子径在 ±1 范围内微调。
    多普勒: 每簇分配一个基准 Doppler，子径在 ±0.3×max_doppler 内变化。

    参数:
      n_paths:     总径数（必须能被 n_clusters 整除）
      n_clusters:  簇数
      max_delay:   最大延迟索引
      max_doppler: 最大 Doppler 索引
      rng:         numpy Generator

    返回:
      delays, dopplers: (n_paths,) arrays
    """
    if rng is None:
        rng = np.random.default_rng()

    pc = n_paths // n_clusters

    # 为每个 cluster 分配互不相同的基准延迟
    base_delays = rng.choice(max_delay + 1, size=n_clusters, replace=False)
    # 为每个 cluster 分配基准 Doppler（随机相位 cos 映射）
    base_phases = rng.uniform(-np.pi, np.pi, n_clusters)
    base_dopplers = max_doppler * np.cos(base_phases)

    delays = np.zeros(n_paths, dtype=int)
    dopplers = np.zeros(n_paths)

    for c in range(n_clusters):
        lo = c * pc
        hi = lo + pc
        # 簇内延迟在基准附近 ±1（裁剪到 [0, max_delay]）
        d_offsets = rng.integers(-1, 2, size=pc)
        delays[lo:hi] = np.clip(base_delays[c] + d_offsets, 0, max_delay)
        # 簇内 Doppler 在基准附近
        dopplers[lo:hi] = base_dopplers[c] + rng.uniform(
            -0.3, 0.3, size=pc
        ) * max_doppler

    return delays, dopplers


# ======================================================================
# Dirichlet 核 + 路径矩阵 + 有效信道
# ======================================================================

def _dirichlet_coeffs(frac, radius, n_sc):
    """
    Dirichlet 核系数。

    对分数 Doppler 偏移 frac ∈ [-0.5, 0.5]，计算展宽半径 radius 内
    各整数偏移 k ∈ [-radius, +radius] 的权重。

    与 MATLAB FractionalChannelBuilder.computeDirichletCoeff 等价。

    返回:
      coeffs: (2*radius+1,) complex128 复数权重（已归一化）
    """
    k_vec = np.arange(-radius, radius + 1)
    alpha = frac - k_vec  # 各偏移对应的分数部分

    num = np.sin(np.pi * alpha)
    den = n_sc * np.sin(np.pi * alpha / n_sc)

    # 避免除零：den≈0 时用极限值
    safe = np.abs(den) > 1e-15
    d = np.empty(len(k_vec), dtype=complex)
    d[safe] = (num[safe] / den[safe]) * np.exp(
        1j * np.pi * alpha[safe] * (n_sc - 1) / n_sc
    )
    d[~safe] = np.exp(1j * np.pi * alpha[~safe] * (n_sc - 1) / n_sc)
    return d


def build_path_matrix(delay, doppler, n_sc, loc_step, radius):
    """
    单径 Dirichlet 核展宽矩阵 Phi_i ∈ C^{N×N}。

    Phi_i[dst, src] 非零当 src + delay ≡ 某 DAFT 域位置，
    dst = (loc_index + j + dk) mod N，权重由 Dirichlet 核给出。

    参数:
      delay:     int 延迟索引
      doppler:   float 分数 Doppler 值
      n_sc:      int 子载波数
      loc_step:  int (2N×c1) mod N
      radius:    int Dirichlet 展宽半径（0 = 无展宽）

    返回:
      Phi: (n_sc, n_sc) complex128 稀疏矩阵
    """
    alpha_int = int(round(doppler))
    frac = doppler - alpha_int
    loc_index = (alpha_int + loc_step * delay) % n_sc

    if radius <= 0:
        # 整数 Doppler：单点映射
        d_coeffs = np.array([1.0 + 0j])
        k_range = np.array([0])
    else:
        k_range = np.arange(-radius, radius + 1)
        d_coeffs = _dirichlet_coeffs(frac, radius, n_sc)

    Phi = np.zeros((n_sc, n_sc), dtype=complex)
    for j in range(n_sc):
        src = (j + delay) % n_sc
        base = (loc_index + j) % n_sc
        for ki, dk in enumerate(k_range):
            dst = (base + dk) % n_sc
            Phi[dst, src] += d_coeffs[ki]
    return Phi


def compute_effective_channel(
    delays, dopplers, gains, n_sc, c1, radius
):
    """
    叠加所有路径构造 DAFT 域有效信道。

    H_eff = sum_i h_i × Phi_i,  Phi_i 由 build_path_matrix 给出。

    参数:
      delays:   (P,) int 延迟索引
      dopplers: (P,) float 分数 Doppler
      gains:    (P,) complex128 径增益
      n_sc:     int 子载波数
      c1:       float chirp 参数 1
      radius:   int Dirichlet 展宽半径

    返回:
      H_eff: (n_sc, n_sc) complex128 有效信道矩阵
    """
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    n_paths = len(delays)
    H = np.zeros((n_sc, n_sc), dtype=complex)
    for i in range(n_paths):
        Phi_i = build_path_matrix(
            int(delays[i]), dopplers[i], n_sc, loc_step, radius
        )
        H += gains[i] * Phi_i
    return H
