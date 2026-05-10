"""
CFAR 有效自由度：从候选点协方差矩阵计算 M_eff 修正独立假设门限。

核心公式（来自四篇精读总结 §CFAR 层）:

    M_eff = (Σλ_i)² / Σλ_i²

    lambda_i 是候选窗检测统计量协方差矩阵的特征值。
    M_eff ∈ [1, M]，完全独立时 M_eff = M，完全相关时 M_eff = 1。

    独立门限: γ_ind = P_res × A_p² × log(M / P_fa)
    修正门限: γ_corr = P_res × A_p² × log(M_eff / P_fa)
    门限修正比: r = log(M_eff / P_fa) / log(M / P_fa) ∈ (0, 1]

    当 r < 1: 独立假设过高估计了门限→更容易漏检。
    当 r ≈ 1: 相关性可忽略→当前门限有效。

参考:
  Raghavan, "CFAR Detection in Clutter With a Kronecker Covariance
  Structure", IEEE TAES, 2017.  (特征值分解, CFAR 门限保持)
  Medeiros et al., "CA-CFAR Detection for SAR Systems Over Correlated
  Gamma-Distributed Clutter", IEEE GRSL, 2024.  (M_eff 修正)
"""

import numpy as np


def compute_effective_candidates(cov_matrix):
    """
    从协方差矩阵计算有效候选数。

    M_eff = (Σ λ_i)² / Σ λ_i²

    参数:
      cov_matrix: (M, M) 候选窗检测统计量的协方差矩阵

    返回:
      m_eff: float 有效候选数 ∈ [1, M]
    """
    cov_matrix = np.asarray(cov_matrix, dtype=float)
    if cov_matrix.ndim != 2 or cov_matrix.shape[0] != cov_matrix.shape[1]:
        raise ValueError("cov_matrix must be a square 2-D matrix")
    if not np.all(np.isfinite(cov_matrix)):
        raise ValueError("cov_matrix contains non-finite values")

    M = cov_matrix.shape[0]
    # 强制对称化以避免数值不对称导致 eigvalsh 偏差
    eigvals = np.linalg.eigvalsh(0.5 * (cov_matrix + cov_matrix.T))
    # 截断负值为零
    eigvals = np.maximum(eigvals, 0.0)
    sum_eig = eigvals.sum()
    if sum_eig < 1e-15:
        return 1.0  # 全零协方差 → 全相关
    m_eff = sum_eig ** 2 / np.sum(eigvals ** 2)
    return float(np.clip(m_eff, 1.0, M))


def compute_independent_threshold(n_candidates, p_fa, residual_power, pilot_amp_sq):
    """
    独立假设下的 CFAR 门限（与 MATLAB GiFreeEstimator 等价）。

    γ_ind = P_res × A_p² × ln(M / P_fa)

    参数:
      n_candidates:  int 候选总数 M = (l_max+1) × (2k_max+1)
      p_fa:          float 目标虚警率（如 1e-3）
      residual_power: float 残差功率 P_res
      pilot_amp_sq:  float 导频幅度平方 A_p²

    返回:
      threshold: float CFAR 门限
    """
    return residual_power * pilot_amp_sq * np.log(n_candidates / p_fa)


def compute_corrected_threshold(
    cov_matrix, p_fa, residual_power, pilot_amp_sq
):
    """
    M_eff 修正后的 CFAR 门限。

    γ_corr = P_res × A_p² × ln(M_eff / P_fa)

    参数:
      cov_matrix:    (M,M) 协方差矩阵
      p_fa:          float
      residual_power: float
      pilot_amp_sq:  float

    返回:
      threshold: float 修正后的门限
    """
    m_eff = compute_effective_candidates(cov_matrix)
    return residual_power * pilot_amp_sq * np.log(m_eff / p_fa)


def compute_threshold_ratio(cov_matrix, p_fa, n_candidates=None):
    """
    门限修正比。

    r = ln(M_eff / P_fa) / ln(M / P_fa)

    r < 1  → 独立假设门限偏高（更容易漏检）
    r ≈ 1 → 相关性可忽略

    参数:
      cov_matrix:   (M,M) 协方差矩阵
      p_fa:         float
      n_candidates: int 若不提供则取 cov_matrix.shape[0]

    返回:
      ratio: float ∈ (0, 1]
      m_eff: float 有效候选数
      m_raw: int   原始候选数
    """
    M = n_candidates if n_candidates is not None else cov_matrix.shape[0]
    m_eff = compute_effective_candidates(cov_matrix)

    log_eff = np.log(max(m_eff, 1e-15) / p_fa)
    log_raw = np.log(M / p_fa)

    ratio = float(log_eff / log_raw)
    return ratio, m_eff, M


def compute_m_eff_from_eigvals(eigvals):
    """
    直接从特征值计算 M_eff（避免重复特征分解）。

    参数:
      eigvals: (M,) 已排序或未排序的特征值

    返回:
      m_eff: float
    """
    eigvals = np.maximum(np.asarray(eigvals), 0.0)
    s = eigvals.sum()
    if s < 1e-15:
        return 1.0
    return float(s ** 2 / np.sum(eigvals ** 2))


def bootstrap_m_eff(power_samples, n_boot=200, seed=0):
    """
    对功率样本进行 bootstrap，估计 M_eff 的经验区间。

    参数:
      power_samples: (n_mc, n_candidates) 非负功率样本
      n_boot:        bootstrap 重采样次数
      seed:          随机种子

    返回:
      summary: dict，包含 mean、p05、p50、p95
    """
    power_samples = np.asarray(power_samples, dtype=float)
    if power_samples.ndim != 2:
        raise ValueError("power_samples must be a 2-D matrix")
    if power_samples.shape[0] < 3:
        raise ValueError("at least three Monte Carlo samples are required")
    if not np.all(np.isfinite(power_samples)):
        raise ValueError("power_samples contains non-finite values")

    rng = np.random.default_rng(seed)
    n_mc = power_samples.shape[0]
    estimates = np.zeros(n_boot)
    for idx in range(n_boot):
        sample_idx = rng.integers(0, n_mc, size=n_mc)
        cov = np.cov(power_samples[sample_idx, :], rowvar=False)
        estimates[idx] = compute_effective_candidates(cov)

    return {
        "mean": float(np.mean(estimates)),
        "p05": float(np.percentile(estimates, 5)),
        "p50": float(np.percentile(estimates, 50)),
        "p95": float(np.percentile(estimates, 95)),
    }
