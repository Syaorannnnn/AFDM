"""
经验 Pfa/Pd 验证。

比较三类 CFAR 门限在候选网格维度的虚警率与检测概率：
  γ_ind  = P_res × ln(M / P_fa)        独立假设门限
  γ_corr = P_res × ln(M_eff / P_fa)    M_eff 修正门限
  γ_cal  = P_res × κ_target             经验标定门限（反解使 Pfa 达标）

P_res 从假候选（非真实路径位置）的统计量均值估计。
本模块只使用 Python 仿真，不调用 MATLAB 主干接收机。
"""

import numpy as np

from cfar import compute_effective_candidates
from grid_stats import estimate_grid_stat_covariance


def _estimate_true_candidate_mask(grid, delays, dopplers):
    """
    根据采样几何标记真实候选。

    分数 Doppler 取整后，如果某候选格点与采样路径的 (delay, round(doppler))
    匹配，则标记为 true path。
    """
    true_pairs = {
        (int(delay), int(round(doppler)))
        for delay, doppler in zip(delays, dopplers)
    }
    mask = np.zeros(grid.shape[0], dtype=bool)
    for idx, (delay, doppler_int) in enumerate(grid):
        mask[idx] = (int(delay), int(doppler_int)) in true_pairs
    return mask


def run_detection_threshold_mc(cfg):
    """
    计算独立门限与 M_eff 修正门限下的经验 Pfa/Pd。
    """
    p_fa = float(cfg.get("p_fa", 1e-3))
    stats = estimate_grid_stat_covariance(cfg)
    power_samples = stats["power_samples"]
    grid = stats["grid"]
    true_mask = _estimate_true_candidate_mask(
        grid, stats["delays"], stats["dopplers"]
    )
    false_mask = ~true_mask

    cov = stats["cov"]
    m_eff = compute_effective_candidates(cov)
    m_raw = grid.shape[0]

    false_power = np.mean(power_samples[:, false_mask])
    residual_power = float(false_power)

    gamma_ind  = residual_power * np.log(m_raw / p_fa)
    gamma_corr = residual_power * np.log(m_eff / p_fa)

    false_stats = power_samples[:, false_mask]
    true_stats  = power_samples[:, true_mask]

    empirical_pfa_ind = float(np.mean(false_stats > gamma_ind))
    empirical_pfa_corr = float(np.mean(false_stats > gamma_corr))
    empirical_pd_ind = float(np.mean(true_stats > gamma_ind))
    empirical_pd_corr = float(np.mean(true_stats > gamma_corr))

    return {
        "m_raw_grid": float(m_raw),
        "m_eff_grid": float(m_eff),
        "threshold_ratio": float(np.log(m_eff / p_fa) / np.log(m_raw / p_fa)),
        "gamma_independent": float(gamma_ind),
        "gamma_corrected": float(gamma_corr),
        "empirical_pfa_independent": empirical_pfa_ind,
        "empirical_pfa_corrected": empirical_pfa_corr,
        "empirical_pd_independent": empirical_pd_ind,
        "empirical_pd_corrected": empirical_pd_corr,
    }


# ======================================================================
# 经验门限标定
# ======================================================================

def calibrate_threshold_multiplier(cfg):
    """
    反解使经验 Pfa 达到目标值的 CFAR 门限乘数 κ。

    流程:
      1. 对 H0 统计量（假候选功率 / residual_power）做经验 CDF
      2. 在 CDF 的 1-P_fa 分位点读取 κ

    κ 的物理含义: γ_cal = P_res × κ,  使 Pfa_empirical = target_pfa

    与公式型门限的关系:
      γ_ind  = P_res × ln(M/P_fa)      → κ_ind  = ln(M/P_fa)
      γ_corr = P_res × ln(M_eff/P_fa)  → κ_corr = ln(M_eff/P_fa)
      γ_cal  = P_res × κ_cal           → κ_cal = 经验反解


    返回:
      dict:
        kappa_cal:          float 标定后的 CFAR 乘数
        kappa_ind:          float ln(M/P_fa)
        kappa_corr:         float ln(M_eff/P_fa)
        kappa_ratio:        float κ_cal / κ_ind
        empirical_pfa_cal:  float 标定门限的经验 Pfa（应 ≈ target）
        empirical_pd_cal:   float 标定门限的经验 Pd
        empirical_pfa_ind:  float
        empirical_pd_ind:   float
    """
    p_fa = float(cfg.get("p_fa", 1e-3))
    stats = estimate_grid_stat_covariance(cfg)
    power_samples = stats["power_samples"]
    grid = stats["grid"]
    true_mask = _estimate_true_candidate_mask(grid, stats["delays"], stats["dopplers"])
    false_mask = ~true_mask

    m_raw = grid.shape[0]
    cov = stats["cov"]
    m_eff = compute_effective_candidates(cov)

    # 残差功率
    false_power = float(np.mean(power_samples[:, false_mask]))
    residual_power = false_power

    # 归一化 H0 统计量: Z_j = power / P_res
    false_norm = power_samples[:, false_mask] / residual_power
    z_all = false_norm.ravel()

    # 经验 CCDF: 排序后从右取分位
    z_sorted = np.sort(z_all)
    # 目标: P(Z > κ) = p_fa
    target_idx = int((1.0 - p_fa) * len(z_sorted))
    target_idx = min(target_idx, len(z_sorted) - 1)
    kappa_cal = float(z_sorted[target_idx])

    # 公式型 κ
    kappa_ind  = np.log(m_raw / p_fa)
    kappa_corr = np.log(m_eff / p_fa)

    # 标定门限下的 Pfa / Pd
    gamma_cal = residual_power * kappa_cal
    false_stats = power_samples[:, false_mask]
    true_stats  = power_samples[:, true_mask]
    gamma_ind  = residual_power * kappa_ind
    gamma_corr = residual_power * kappa_corr

    return {
        "kappa_cal": kappa_cal,
        "kappa_ind": kappa_ind,
        "kappa_corr": kappa_corr,
        "kappa_ratio": float(kappa_cal / kappa_ind),
        "m_raw_grid": float(m_raw),
        "m_eff_grid": float(m_eff),
        "empirical_pfa_cal": float(np.mean(false_stats > gamma_cal)),
        "empirical_pd_cal":  float(np.mean(true_stats > gamma_cal)),
        "empirical_pfa_ind": float(np.mean(false_stats > gamma_ind)),
        "empirical_pd_ind":  float(np.mean(true_stats > gamma_ind)),
        "empirical_pfa_corr": float(np.mean(false_stats > gamma_corr)),
        "empirical_pd_corr":  float(np.mean(true_stats > gamma_corr)),
    }


def run_calibration_matrix(base_cfg):
    """
    扫描 ρ × dir_r 矩阵，全部用经验标定。

    扫描维度:
      rho: [0.0, 0.5, 0.9]
      dir_r: [0, 1, 3]

    返回:
      rows: list[dict]，每行含 rho, dirichlet_r, kappa_cal, kappa_ind,
                        kappa_ratio, pfa/cal/ind/corr, pd_cal
    """
    rows = []
    for rho in [0.0, 0.5, 0.9]:
        for radius in [0, 1, 3]:
            _cfg = dict(base_cfg)
            _cfg["rho"] = rho
            _cfg["dirichlet_r"] = radius
            result = calibrate_threshold_multiplier(_cfg)
            result["rho"] = rho
            result["dirichlet_r"] = radius
            rows.append(result)
    return rows
