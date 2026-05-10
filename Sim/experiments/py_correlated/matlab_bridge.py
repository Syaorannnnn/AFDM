#!/usr/bin/env python3
"""
MATLAB pyrunfile() 桥接入口。

用法 (MATLAB):
  cfg = py.dict(pyargs('n_sc', 64, 'n_paths', 6, 'rho', 0.5));
  result = pyrunfile('Sim/experiments/py_correlated/matlab_bridge.py', ...
                     'result', config=cfg);
  disp(result{'threshold_ratio'});

输入 config 键 (全部可选, 有默认值):
  n_sc:           int   子载波数 (默认 64)
  n_paths:        int   总径数 (默认 6)
  n_clusters:     int   簇数 (默认 2)
  rho:            float BDLR ρ (默认 0.5)
  max_delay:      int   最大延迟索引 (默认 3)
  max_doppler:    int   最大 Doppler 索引 (默认 2)
  doppler_guard:  int   Doppler 保护间隔 (默认 3)
  dirichlet_r:    int   Dirichlet 展宽半径 (默认 3)
  mc_samples:     int   Monte Carlo 次数 (默认 500, 标定时建议 >= 5000)
  seed:           int   随机种子 (默认 42)
  p_fa:           float CFAR P_fa (默认 1e-3)
  do_calibrate:   float 设为 1.0 时执行 κ 标定（需更多 MC）

返回 result (py.dict):
  'm_eff':             子载波维度有效候选数
  'm_raw':             原始候选数 M = N-1
  'threshold_ratio':   门限修正比 ∈ (0, 1]
  'n_null':            null 子载波数 (Delta<0)
  'n_hotspot':         hotspot 子载波数 (Delta>0)
  'delta_min':         Delta_total 最小值
  'delta_max':         Delta_total 最大值
  'm_eff_grid':        候选网格维度有效候选数
  'kappa_ind':         公式 κ (ln(M/P_fa))
  'kappa_corr':        M_eff 修正 κ (ln(M_eff/P_fa))
  [若 do_calibrate=1.0]:
    'kappa_cal':       经验标定 κ（使 Pfa 达标）
    'kappa_ratio':     κ_cal / κ_ind（当前 CLIP s 因子的对标值）
    'pd_at_calibrated': 标定门限下的 Pd
  'ok':                1.0 表示正常运行
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from channel import (
    compute_chirp_params, sample_clustered_geometry,
    build_path_matrix,
)
from id2p import (
    compute_id2p_spectrum, compute_effective_id2p,
    compute_null_hotspot_map, compute_spatial_cov,
)
from cfar import compute_threshold_ratio, compute_effective_candidates
from grid_stats import estimate_grid_stat_covariance
from detection_mc import _estimate_true_candidate_mask


def matlab_bridge(config):
    """
    供 MATLAB 调用的单入口函数。
    """
    n_sc         = int(config.get('n_sc', 64))
    n_paths      = int(config.get('n_paths', 6))
    n_clusters   = int(config.get('n_clusters', 2))
    rho          = float(config.get('rho', 0.5))
    max_delay    = int(config.get('max_delay', 3))
    max_doppler  = int(config.get('max_doppler', 2))
    doppler_guard = int(config.get('doppler_guard', 3))
    dirichlet_r  = int(config.get('dirichlet_r', 3))
    mc_samples   = int(config.get('mc_samples', 500))
    seed         = int(config.get('seed', 42))
    p_fa         = float(config.get('p_fa', 1e-3))
    do_calibrate = bool(int(float(config.get('do_calibrate', 0.0))))

    try:
        n_data = n_sc - 1
        pilot_idx = 0
        data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
        gain_var = 1.0 / n_paths
        rng = np.random.default_rng(seed)
        c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)

        # 固定簇化几何
        delays, dopplers = sample_clustered_geometry(
            n_paths, n_clusters, max_delay, max_doppler, rng
        )
        cluster_labels = np.repeat(np.arange(n_clusters), n_paths // n_clusters)

        Phi_list = [
            build_path_matrix(int(delays[i]), dopplers[i], n_sc, loc_step,
                              dirichlet_r)
            for i in range(n_paths)
        ]

        _, Dc, B_total, D_total = compute_id2p_spectrum(
            Phi_list, cluster_labels, n_data, pilot_idx, data_idx
        )
        _, _, n_null, n_hotspot = compute_null_hotspot_map(D_total)

        cov = compute_spatial_cov(
            delays, dopplers, cluster_labels, n_sc, c1,
            dirichlet_r, pilot_idx, data_idx,
            mc_samples, rho, n_clusters, rng
        )
        ratio, m_eff, m_raw = compute_threshold_ratio(cov, p_fa, n_data)

        # --- 候选网格维度 ---
        grid_cfg = {
            "n_sc": n_sc, "n_paths": n_paths, "n_clusters": n_clusters,
            "max_delay": max_delay, "max_doppler": max_doppler,
            "doppler_guard": doppler_guard, "dirichlet_r": dirichlet_r,
            "rho": rho, "mc_samples": mc_samples, "seed": seed,
            "pilot_idx": pilot_idx,
        }
        grid_stats = estimate_grid_stat_covariance(grid_cfg)
        m_eff_grid = compute_effective_candidates(grid_stats["cov"])
        m_raw_grid = float(grid_stats["grid"].shape[0])

        kappa_ind  = np.log(m_raw_grid / p_fa)
        kappa_corr = np.log(m_eff_grid / p_fa)

        result = {
            "m_eff": float(m_eff),
            "m_raw": float(m_raw),
            "threshold_ratio": float(ratio),
            "n_null": float(n_null),
            "n_hotspot": float(n_hotspot),
            "delta_min": float(D_total.min()),
            "delta_max": float(D_total.max()),
            "m_eff_grid": float(m_eff_grid),
            "m_raw_grid": float(m_raw_grid),
            "kappa_ind": float(kappa_ind),
            "kappa_corr": float(kappa_corr),
            "ok": 1.0,
        }

        # --- 经验 κ 标定 ---
        if do_calibrate:
            true_mask = _estimate_true_candidate_mask(
                grid_stats["grid"], grid_stats["delays"], grid_stats["dopplers"]
            )
            false_mask = ~true_mask
            power_samples = grid_stats["power_samples"]

            false_power = float(np.mean(power_samples[:, false_mask]))
            residual_power = max(false_power, 1e-15)

            false_norm = (power_samples[:, false_mask] / residual_power).ravel()
            z_sorted = np.sort(false_norm)
            target_idx = int((1.0 - p_fa) * len(z_sorted))
            target_idx = min(target_idx, len(z_sorted) - 1)
            kappa_cal = float(z_sorted[target_idx])

            gamma_cal = residual_power * kappa_cal
            result["kappa_cal"] = kappa_cal
            result["kappa_ratio"] = float(kappa_cal / kappa_ind)
            result["pd_at_calibrated"] = float(
                np.mean(power_samples[:, true_mask] > gamma_cal)
            )

        return result

    except Exception as e:
        return {
            "m_eff": 0.0, "m_raw": 0.0, "threshold_ratio": 0.0,
            "n_null": 0.0, "n_hotspot": 0.0,
            "delta_min": 0.0, "delta_max": 0.0,
            "m_eff_grid": 0.0, "m_raw_grid": 0.0,
            "kappa_ind": 0.0, "kappa_corr": 0.0,
            "ok": 0.0, "error": str(e),
        }


if __name__ == "__main__":
    print("=== matlab_bridge 快速测试 ===")
    for do_cal in [0.0, 1.0]:
        test_cfg = {
            "n_sc": 64, "n_paths": 6, "n_clusters": 2,
            "rho": 0.5, "mc_samples": 2000,
            "do_calibrate": do_cal,
        }
        r = matlab_bridge(test_cfg)
        print(f"\n  do_calibrate={int(do_cal)}:")
        for k in sorted(r):
            print(f"    {k}: {r[k]}")
