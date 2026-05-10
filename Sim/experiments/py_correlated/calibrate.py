#!/usr/bin/env python3
"""
CFAR 经验门限标定：反解使 Pfa 达标所需的 κ 乘数。

用法:
  python calibrate.py
  python calibrate.py --mc 20000 --seed 1

输出:
  calibrate_matrix.csv — κ_cal vs κ_ind/vs κ_corr 全矩阵
  Fig5e_ThresholdCalibration.png — κ vs ρ 曲线
"""

import argparse
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _plot_setup import setup_plots
setup_plots()
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from detection_mc import run_calibration_matrix, calibrate_threshold_multiplier
from cfar import compute_effective_candidates, bootstrap_m_eff
from grid_stats import estimate_grid_stat_covariance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=str, default=None)
    args = ap.parse_args()

    # 输出目录
    sim_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = args.out_dir or os.path.join(sim_dir, "Results", "Figures")
    os.makedirs(out_dir, exist_ok=True)

    base = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "rho": 0.5, "mc_samples": args.mc,
        "seed": args.seed, "p_fa": 1e-3, "pilot_idx": 0,
    }

    t0 = time.time()
    print(f"=== CFAR Threshold Calibration Matrix (mc={args.mc}) ===")

    rows = run_calibration_matrix(base)

    # --- 表格 ---
    header = ("rho,dir_r,M_raw,M_eff,"
              "kappa_ind,kappa_corr,kappa_cal,kappa_ratio,"
              "pfa_ind,pfa_corr,pfa_cal,pd_ind,pd_corr,pd_cal")
    print(f"\n{'rho':>5} {'r':>3} {'M_raw':>5} {'M_eff':>6} "
          f"{'κ_ind':>7} {'κ_corr':>6} {'κ_cal':>7} {'ratio':>6} "
          f"{'Pfa_ind':>9} {'Pfa_corr':>9} {'Pfa_cal':>9} "
          f"{'Pd_ind':>7} {'Pd_corr':>7} {'Pd_cal':>7}")
    print("-" * 100)

    csv_path = os.path.join(out_dir, "calibrate_matrix.csv")
    with open(csv_path, "w") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(
                f"{r['rho']},{r['dirichlet_r']},{r['m_raw_grid']},{r['m_eff_grid']:.3f},"
                f"{r['kappa_ind']:.3f},{r['kappa_corr']:.3f},{r['kappa_cal']:.3f},{r['kappa_ratio']:.3f},"
                f"{r['empirical_pfa_ind']:.5f},{r['empirical_pfa_corr']:.5f},{r['empirical_pfa_cal']:.5f},"
                f"{r['empirical_pd_ind']:.4f},{r['empirical_pd_corr']:.4f},{r['empirical_pd_cal']:.4f}\n"
            )
            print(
                f"{r['rho']:5.1f} {r['dirichlet_r']:3d} "
                f"{r['m_raw_grid']:5.0f} {r['m_eff_grid']:6.2f} "
                f"{r['kappa_ind']:7.2f} {r['kappa_corr']:6.2f} {r['kappa_cal']:7.2f} {r['kappa_ratio']:6.2f} "
                f"{r['empirical_pfa_ind']:9.5f} {r['empirical_pfa_corr']:9.5f} {r['empirical_pfa_cal']:9.5f} "
                f"{r['empirical_pd_ind']:7.4f} {r['empirical_pd_corr']:7.4f} {r['empirical_pd_cal']:7.4f}"
            )

    t1 = time.time()
    print(f"\nCSV saved: {csv_path}")
    print(f"Elapsed: {t1 - t0:.1f}s")

    # --- 图: κ vs ρ ---
    _plot_kappa_comparison(rows, out_dir)


def _plot_kappa_comparison(rows, out_dir):
    """Fig5e: κ_cal vs κ_ind vs κ_corr 对比。"""
    # 只取 dir_r=3
    r3 = [r for r in rows if r["dirichlet_r"] == 3]
    rhos = np.array([r["rho"] for r in r3])

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # 左: κ 三条曲线
    ax = axes[0]
    k_ind  = np.array([r["kappa_ind"] for r in r3])
    k_corr = np.array([r["kappa_corr"] for r in r3])
    k_cal  = np.array([r["kappa_cal"] for r in r3])
    ax.plot(rhos, k_ind, 's-', color='gray', linewidth=2, markersize=8, label=r'$\kappa_{\rm ind}=\ln(M/P_{fa})$')
    ax.plot(rhos, k_corr, 'D--', color='steelblue', linewidth=2, markersize=8, label=r'$\kappa_{\rm corr}=\ln(M_{\rm eff}/P_{fa})$')
    ax.plot(rhos, k_cal, 'o-', color='crimson', linewidth=2, markersize=10, markerfacecolor='crimson', label=r'$\kappa_{\rm cal}$ (经验标定)')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$\kappa$ (CFAR multiplier)')
    ax.set_title('CFAR multiplier comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 中: κ ratio
    ax = axes[1]
    for r_val, marker, label in [(0, 'steelblue', 'dir_r=0'), (1, 'green', 'dir_r=1'), (3, 'crimson', 'dir_r=3')]:
        subset = [r for r in rows if r["dirichlet_r"] == r_val]
        rho_v = [r["rho"] for r in subset]
        ratio_v = [r["kappa_ratio"] for r in subset]
        ax.plot(rho_v, ratio_v, 'o-', color=marker, linewidth=2, markersize=8, label=label)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$\kappa_{\rm cal} / \kappa_{\rm ind}$')
    ax.set_title(r'Calibration ratio: how much $\kappa$ must increase')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 右: Pd 三组对比 (dir_r=3)
    ax = axes[2]
    x = np.arange(3)
    width = 0.25
    pd_ind  = np.array([r["empirical_pd_ind"] for r in r3])
    pd_corr = np.array([r["empirical_pd_corr"] for r in r3])
    pd_cal  = np.array([r["empirical_pd_cal"] for r in r3])
    ax.bar(x - width, pd_ind,  width, color='gray',      label=r'$P_d$ 独立')
    ax.bar(x,          pd_corr, width, color='steelblue', label=r'$P_d$ $M_{\rm eff}$')
    ax.bar(x + width,  pd_cal,  width, color='crimson',  label=r'$P_d$ 标定')
    ax.set_xticks(x)
    ax.set_xticklabels([r'$\rho=0.0$', r'$\rho=0.5$', r'$\rho=0.9$'])
    ax.set_ylabel(r'$P_d$')
    ax.set_title('Detection probability comparison (dir_r=3)')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle(
        'CFAR Threshold Calibration under Fractional Doppler (dir_r=3)',
        fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5e_ThresholdCalibration.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5e saved: {path}')


if __name__ == '__main__':
    main()
