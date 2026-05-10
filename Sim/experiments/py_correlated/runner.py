#!/usr/bin/env python3
"""
ID2P 相关径推广 —— 实验编排与可视化入口。

用法:
  # 完整 rho 扫描 + 功率分配测试 (默认参数)
  python runner.py --all

  # 小规模快速诊断
  python runner.py --quick --n-sc 64 --n-paths 6 --rho 0.5

  # 仅 CFAR 修正分析
  python runner.py --cfar-only --rho 0.8

  # MATLAB 调用 (静默模式, 仅输出 CSV)
  python runner.py --quick --matlab-out

依赖: numpy, scipy, matplotlib
"""

import argparse
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
# 确保同目录模块可导入
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _plot_setup import setup_plots
setup_plots()
import matplotlib.pyplot as plt

# 将包目录加入路径, 支持直接脚本运行
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from channel import (
    compute_chirp_params, sample_clustered_geometry, bdlr_gains,
    build_path_matrix, compute_effective_channel,
)
from id2p import (
    compute_id2p_spectrum, compute_effective_id2p,
    compute_null_hotspot_map, compute_spatial_cov,
)
from cfar import (
    compute_effective_candidates, compute_threshold_ratio,
    compute_independent_threshold, compute_corrected_threshold,
)
from detection_mc import run_detection_threshold_mc
from report import write_markdown_report


# ======================================================================
# 默认参数
# ======================================================================

DEFAULTS = {
    'n_sc':          64,
    'max_delay':     3,
    'max_doppler':   2,
    'doppler_guard': 3,
    'dirichlet_r':   3,
    'n_paths':       6,
    'n_clusters':    2,
    'rho_low':       0.05,
    'rho_high':      0.85,
    'mc_gains':      500,
    'pilot_idx':     0,
    'seed':          42,
    'p_fa':          1e-3,
    'gain_var':      None,   # 默认 = 1/n_paths
    'snr_db_vec':    [0, 5, 10, 15, 20, 25],
}


# ======================================================================
# rho 扫描核心
# ======================================================================

def run_rho_sweep(cfg):
    """
    rho 扫描：计算 M_eff、门限修正比、null/hotspot 分布。

    返回:
      dict 含 'rho_vec', 'm_eff_vec', 'ratio_vec', 'n_null_vec',
      'n_hotspot_vec', 'delta_range_vec', 'cov_low', 'cov_high',
      'B_total', 'D_total', ...
    """
    n_sc       = cfg['n_sc']
    n_paths    = cfg['n_paths']
    n_clusters = cfg['n_clusters']
    rho_low    = cfg['rho_low']
    rho_high   = cfg['rho_high']
    mc_gains   = cfg['mc_gains']
    pilot_idx  = cfg['pilot_idx']
    seed       = cfg['seed']

    gain_var = cfg['gain_var'] or (1.0 / n_paths)
    pc = n_paths // n_clusters
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, cfg['max_doppler'], cfg['doppler_guard'])

    # --- 固定簇化几何 ---
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, cfg['max_delay'], cfg['max_doppler'], rng
    )
    cluster_labels = np.repeat(np.arange(n_clusters), pc)

    Phi_list = [
        build_path_matrix(int(delays[i]), dopplers[i], n_sc, loc_step,
                          cfg['dirichlet_r'])
        for i in range(n_paths)
    ]

    Bc, Dc, B_total, D_total = compute_id2p_spectrum(
        Phi_list, cluster_labels, n_data, pilot_idx, data_idx
    )

    # --- rho 扫描: 协方差 + M_eff ---
    rho_vec = np.array([0.0, 0.1, 0.3, 0.5, 0.7, 0.9])
    m_eff_vec = np.zeros(len(rho_vec))
    n_null_vec = np.zeros(len(rho_vec), dtype=int)
    n_hotspot_vec = np.zeros(len(rho_vec), dtype=int)
    delta_range_vec = np.zeros(len(rho_vec))
    cov_low = None   # rho_vec[0] 的协方差
    cov_high = None  # rho_vec[-1] 的协方差

    for ri, rho in enumerate(rho_vec):
        cov = compute_spatial_cov(
            delays, dopplers, cluster_labels, n_sc, c1,
            cfg['dirichlet_r'], pilot_idx, data_idx,
            mc_gains, rho, n_clusters, rng
        )
        m_eff_vec[ri] = compute_effective_candidates(cov)
        if ri == 0:
            cov_low = cov
        if ri == len(rho_vec) - 1:
            cov_high = cov

        # rho 下的期望 ID2P
        eff = compute_effective_id2p(B_total, D_total, rho, gain_var)
        is_null, is_hot, nn, nh = compute_null_hotspot_map(D_total)
        n_null_vec[ri] = nn
        n_hotspot_vec[ri] = nh
        delta_range_vec[ri] = D_total.max() - D_total.min()

    return {
        'rho_vec': rho_vec,
        'm_eff_vec': m_eff_vec,
        'n_null_vec': n_null_vec,
        'n_hotspot_vec': n_hotspot_vec,
        'delta_range_vec': delta_range_vec,
        'cov_low': cov_low,
        'cov_high': cov_high,
        'B_total': B_total,
        'D_total': D_total,
        'B_cluster': Bc,
        'D_cluster': Dc,
        'delays': delays,
        'dopplers': dopplers,
        'cluster_labels': cluster_labels,
        'n_data': n_data,
        'data_idx': data_idx,
        'cfg': cfg,
    }


# ======================================================================
# 功率分配测试
# ======================================================================

def run_power_allocation_test(cfg, B_total, D_total):
    """
    对比 null-weighted vs uniform 功率分配的有效互信息。

    简化模型: MI = Σ_j log₂(1 + w_j / N0_eff)
    其中 N0_eff = N0 + P_ID2P, P_ID2P = Σ_j w_j × σ_j²(ρ)

    返回:
      dict 含 'snr_db_vec', 'mi_uniform', 'mi_null', 'w_null', ...
    """
    rho = cfg['rho_high']
    gain_var = cfg['gain_var'] or (1.0 / cfg['n_paths'])
    snr_db_vec = cfg['snr_db_vec']
    n_data = len(B_total)

    sigma2_rho = gain_var * (B_total + rho * D_total)
    w_uniform = np.ones(n_data)
    w_null = 1.0 / (sigma2_rho + 1e-4 * sigma2_rho.mean())

    mi_uniform = np.zeros(len(snr_db_vec))
    mi_null = np.zeros(len(snr_db_vec))
    for si, snr_db in enumerate(snr_db_vec):
        snr_lin = 10.0 ** (snr_db / 10.0)
        mi_uniform[si] = _compute_mi(w_uniform, sigma2_rho, snr_lin)
        mi_null[si] = _compute_mi(w_null, sigma2_rho, snr_lin)

    return {
        'snr_db_vec': np.array(snr_db_vec),
        'mi_uniform': mi_uniform,
        'mi_null': mi_null,
        'w_null': w_null,
        'sigma2_rho': sigma2_rho,
    }


def _compute_mi(weights, sigma2_vec, total_snr_lin):
    """有效互信息: MI = Σ_j log₂(1 + w_j / N0_eff)。"""
    w_norm = weights / weights.sum() * len(weights)
    n0_eff = 1.0 / total_snr_lin + np.sum(w_norm * sigma2_vec)
    snr_eff = w_norm / n0_eff
    return float(np.sum(np.log2(1.0 + snr_eff)))


# ======================================================================
# 绘图
# ======================================================================

def plot_all(sweep_results, pa_results, out_dir):
    """生成全套 Fig5 系列图。"""
    os.makedirs(out_dir, exist_ok=True)

    _plot_delta_decomp(sweep_results, out_dir)
    _plot_spatial_cov(sweep_results, out_dir)
    _plot_power_allocation(pa_results, sweep_results['cfg'], out_dir)
    _plot_cfar_correction(sweep_results, out_dir)
    plt.close('all')


def _plot_delta_decomp(res, out_dir):
    """Fig5a: Delta 分解热力图 + 直方图 + 排序剖面。"""
    Bc = res['B_cluster']
    Dc = res['D_cluster']
    Bt = res['B_total']
    Dt = res['D_total']
    n_data = res['n_data']
    n_clusters = Bc.shape[0]
    delta_range = max(abs(Dt.min()), abs(Dt.max()), 1e-10)

    fig, axes = plt.subplots(2, 3, figsize=(20, 10))

    for ci in range(min(n_clusters, 2)):
        ax = axes[0, ci]
        im = ax.imshow(Dc[ci].reshape(1, -1), aspect='auto', cmap='RdBu_r',
                       vmin=-delta_range, vmax=delta_range,
                       extent=[0, n_data, 0, 1])
        ax.set_title(fr'Cluster {ci}: $\Delta_c(j)$', fontsize=12)
        ax.set_yticks([])

    ax = axes[0, 2]
    im = ax.imshow(Dt.reshape(1, -1), aspect='auto', cmap='RdBu_r',
                   vmin=-delta_range, vmax=delta_range,
                   extent=[0, n_data, 0, 1])
    ax.set_title(r'Total $\sum_c \Delta_c(j)$', fontsize=12)
    ax.set_yticks([])
    fig.colorbar(im, ax=axes[0, :], orientation='horizontal', pad=0.15,
                 shrink=0.6, label=r'$\Delta_c(j)$')

    # hist
    ax = axes[1, 0]
    ax.hist(Dt, bins=30, color='steelblue', edgecolor='k', alpha=0.8)
    ax.axvline(0, color='red', linestyle='--', linewidth=1.2)
    n_null = int((Dt < 0).sum())
    n_hot = int((Dt > 0).sum())
    ax.set_xlabel(r'$\sum_c \Delta_c(j)$')
    ax.set_ylabel('Count')
    ax.set_title('Null/hotspot distribution')
    ax.text(0.98, 0.95, f'null: {n_null}  hotspot: {n_hot}',
            transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    # sorted Delta
    ax = axes[1, 1]
    sorted_idx = np.argsort(Dt)
    ax.plot(Dt[sorted_idx], '.-', color='steelblue', markersize=3)
    ax.axhline(0, color='red', linestyle='--', linewidth=1.2)
    ax.set_xlabel('Subcarrier index (sorted by Δ)')
    ax.set_ylabel(r'$\sum_c \Delta_c(j)$')
    ax.set_title('Sorted Delta profile')

    # B + D overlay
    ax = axes[1, 2]
    ax2 = ax.twinx()
    ax.plot(Bt, '.-', color='gray', markersize=3, alpha=0.6,
            label=r'$B_{total}$ (baseline)')
    ax2.plot(Dt, '.-', color='steelblue', markersize=3, alpha=0.8,
             label=r'$D_{total}$ (modulation)')
    ax.set_xlabel('Data subcarrier index j')
    ax.set_ylabel(r'$B_{total}(j)$', color='gray')
    ax2.set_ylabel(r'$D_{total}(j)$', color='steelblue')
    ax.set_title('Baseline vs modulation term')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    cfg = res['cfg']
    fig.suptitle(
        f'Delta-Decomposition (N={cfg["n_sc"]}, Dirichlet r={cfg["dirichlet_r"]}, '
        f'fractional Doppler)',
        fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5a_DeltaDecomp.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5a saved: {path}')


def _plot_spatial_cov(res, out_dir):
    """Fig5b: 空间协方差矩阵 (low-ρ vs high-ρ)。"""
    cov_low = res['cov_low']
    cov_high = res['cov_high']
    cfg = res['cfg']

    if cov_low is None or cov_high is None:
        print('  Fig5b SKIP: no cov data')
        return

    vlim = max(
        np.percentile(np.abs(cov_high), 98),
        np.percentile(np.abs(cov_low), 98),
    )
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    for ax, cov, ri, title in [
        (axes[0], cov_low, 0, 'low'),
        (axes[1], cov_high, -1, 'high'),
    ]:
        im = ax.imshow(cov, cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                       aspect='equal', interpolation='nearest')
        ax.set_title(
            rf'Spatial covariance $\rho={res["rho_vec"][ri]:.1f}$ ({title})',
            fontsize=13)
        ax.set_xlabel("Data subcarrier j'")
        ax.set_ylabel('Data subcarrier j')
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle(r'$\mathrm{Cov}(|H(0,j)|^2, |H(0,j^\prime)|^2)$ — Cluster Imprint',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5b_SpatialCov.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5b saved: {path}')


def _plot_power_allocation(pa, cfg, out_dir):
    """Fig5c: null-weighted 功率分配 vs uniform。"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax = axes[0]
    ax.plot(pa['w_null'], '.-', color='steelblue', markersize=3)
    ax.axhline(1, color='gray', linestyle='--', linewidth=0.8, label='uniform')
    ax.set_xlabel('Data subcarrier index j')
    ax.set_ylabel(r'Power weight $w_j$')
    ax.set_title('Null-weighted power allocation')
    ax.legend()

    ax = axes[1]
    ax.plot(pa['snr_db_vec'], pa['mi_uniform'], 's-', color='gray',
            linewidth=2, markersize=8, label='Uniform power')
    ax.plot(pa['snr_db_vec'], pa['mi_null'], 'o-', color='steelblue',
            linewidth=2, markersize=8, label='Null-weighted')
    ax.set_xlabel('Total data SNR (dB)')
    ax.set_ylabel('Effective MI (bits/channel use)')
    ax.set_title('Power allocation: Uniform vs Null-weighted')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f'Power Allocation Test (rho={cfg["rho_high"]}, N={cfg["n_sc"]})',
        fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5c_PowerAlloc.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5c saved: {path}')


def _plot_cfar_correction(res, out_dir):
    """Fig5d: M_eff / M 修正比 vs ρ。"""
    rho_vec = res['rho_vec']
    m_eff_vec = res['m_eff_vec']
    M = res['n_data']

    p_fa = res['cfg']['p_fa']
    ratio_vec = np.array([
        np.log(max(me, 1e-15) / p_fa) / np.log(M / p_fa)
        for me in m_eff_vec
    ])

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax = axes[0]
    ax.plot(rho_vec, m_eff_vec, 'o-', color='steelblue', linewidth=2,
            markersize=8, markerfacecolor='steelblue')
    ax.axhline(M, color='gray', linestyle='--', alpha=0.5,
               label=f'M={M} (独立上限)')
    ax.axhline(1, color='red', linestyle=':', alpha=0.5,
               label='M=1 (全相关下界)')
    ax.set_xlabel(r'$\rho$ (intra-cluster correlation)')
    ax.set_ylabel(r'$M_{\rm eff}$')
    ax.set_title(r'Effective candidates $M_{\rm eff}$ vs $\rho$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(rho_vec, ratio_vec, 's-', color='crimson', linewidth=2,
            markersize=8, markerfacecolor='crimson')
    ax.axhline(1, color='gray', linestyle='--', alpha=0.5, label='ratio=1')
    ax.set_xlabel(r'$\rho$ (intra-cluster correlation)')
    ax.set_ylabel(r'Threshold ratio $\gamma_{\rm corr} / \gamma_{\rm ind}$')
    ax.set_title('CFAR threshold correction ratio')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f'CFAR Correction Analysis (M={M}, P_fa={p_fa})',
        fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5d_CfarCorrection.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5d saved: {path}')


# ======================================================================
# CSV 导出
# ======================================================================

def save_matlab_csv(sweep_results, pa_results, out_dir):
    """导出数值报告 CSV 供 MATLAB 后处理。"""
    os.makedirs(out_dir, exist_ok=True)

    # 主表: rho 扫描
    path_main = os.path.join(out_dir, 'rho_sweep.csv')
    header = 'rho,m_eff,threshold_ratio,n_null,n_hotspot,delta_range\n'
    rows = []
    rho = sweep_results['rho_vec']
    me = sweep_results['m_eff_vec']
    M = sweep_results['n_data']
    p_fa = sweep_results['cfg']['p_fa']
    for i in range(len(rho)):
        ratio = np.log(max(me[i], 1e-15)/p_fa) / np.log(M/p_fa)
        rows.append(
            f'{rho[i]:.2f},{me[i]:.2f},{ratio:.4f},'
            f'{sweep_results["n_null_vec"][i]},'
            f'{sweep_results["n_hotspot_vec"][i]},'
            f'{sweep_results["delta_range_vec"][i]:.6f}'
        )
    with open(path_main, 'w') as f:
        f.write(header + '\n'.join(rows) + '\n')
    print(f'  CSV saved: {path_main}')

    # 功率分配表
    path_pa = os.path.join(out_dir, 'power_allocation.csv')
    with open(path_pa, 'w') as f:
        f.write('snr_db,mi_uniform,mi_null,mi_gain\n')
        for si, snr in enumerate(pa_results['snr_db_vec']):
            gain = pa_results['mi_null'][si] - pa_results['mi_uniform'][si]
            f.write(
                f'{snr},{pa_results["mi_uniform"][si]:.4f},'
                f'{pa_results["mi_null"][si]:.4f},{gain:+.4f}\n'
            )
    print(f'  CSV saved: {path_pa}')


# ======================================================================
# 命令行入口
# ======================================================================

def build_config(args):
    cfg = dict(DEFAULTS)
    for key in DEFAULTS:
        val = getattr(args, key, None)
        if val is not None:
            cfg[key] = val
    cfg['gain_var'] = 1.0 / cfg['n_paths']
    return cfg


def main():
    ap = argparse.ArgumentParser(
        description='ID2P 相关径推广 — 实验运行器'
    )
    ap.add_argument('--all', action='store_true',
                    help='运行完整扫描')
    ap.add_argument('--quick', action='store_true',
                    help='小规模快速诊断 (N=64, mc=100)')
    ap.add_argument('--cfar-only', action='store_true',
                    help='仅 CFAR 修正分析')
    ap.add_argument('--matlab-out', action='store_true',
                    help='静默模式, 仅输出 CSV')
    ap.add_argument('--n-sc', type=int, default=None)
    ap.add_argument('--n-paths', type=int, default=None)
    ap.add_argument('--n-clusters', type=int, default=None)
    ap.add_argument('--rho', type=float, default=None,
                    help='rho_high 值')
    ap.add_argument('--out-dir', type=str, default=None)
    ap.add_argument('--seed', type=int, default=None)
    ap.add_argument("--grid-cfar", action="store_true",
                    help="运行延迟-Doppler 候选网格上的 M_eff/Pfa/Pd 验证")
    ap.add_argument("--mc-grid", type=int, default=None,
                    help="候选网格 Monte Carlo 次数")
    ap.add_argument("--scan-grid", action="store_true",
                    help="执行 Python 级 rho/dirichlet/geometry 扫描")

    args = ap.parse_args()
    cfg = build_config(args)

    if args.rho is not None:
        cfg['rho_high'] = args.rho
    if args.n_clusters is not None:
        cfg['n_clusters'] = args.n_clusters
    if args.quick:
        cfg['mc_gains'] = 100
    if args.seed is not None:
        cfg['seed'] = args.seed
    if args.mc_grid is not None:
        cfg["mc_samples_grid"] = int(args.mc_grid)

    # 输出目录: py_correlated/ -> experiments/ -> CODES/Sim/
    sim_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = args.out_dir or os.path.join(sim_dir, 'Results', 'Figures')
    os.makedirs(out_dir, exist_ok=True)

    verbose = not args.matlab_out
    if verbose:
        print(f'=== ID2P Correlation Probe ===')
        print(f'  N={cfg["n_sc"]}, P={cfg["n_paths"]}, '
              f'C={cfg["n_clusters"]}, mc={cfg["mc_gains"]}')
        print(f'  rho_high={cfg["rho_high"]}')
        t0 = time.time()

    # --- 运行 ---
    sweep = run_rho_sweep(cfg)
    pa = run_power_allocation_test(cfg, sweep['B_total'], sweep['D_total'])

    if verbose:
        t1 = time.time()
        print(f'  Sweep done ({t1-t0:.1f}s)')

    # --- 输出 ---
    if not args.cfar_only:
        plot_all(sweep, pa, out_dir)

    save_matlab_csv(sweep, pa, out_dir)

    # --- 候选网格 CFAR 验证 ---
    grid_result = None
    grid_rows = []
    if args.grid_cfar:
        grid_cfg = dict(cfg)
        grid_cfg["rho"] = cfg["rho_high"]
        grid_cfg["mc_samples"] = int(cfg.get("mc_samples_grid", cfg["mc_gains"]))
        grid_result = run_detection_threshold_mc(grid_cfg)
        if verbose:
            print("\n=== Grid-CFAR Monte Carlo ===")
            for key in sorted(grid_result):
                print(f"  {key}: {grid_result[key]:.6g}")
        row = dict(grid_result)
        row["rho"] = float(cfg["rho_high"])
        grid_rows.append(row)

    # --- 参数扫描矩阵 ---
    if args.scan_grid:
        scan_rows = run_python_validation_matrix(cfg)
        scan_path = save_validation_matrix_csv(scan_rows, out_dir)
        if verbose:
            print(f"\n  Scan CSV saved: {scan_path}")
        for srow in scan_rows:
            grid_rows.append(srow)

    # --- Markdown 报告 ---
    id2p_rows = []
    for i, rho_v in enumerate(sweep["rho_vec"]):
        ratio = np.log(max(sweep["m_eff_vec"][i], 1e-15) / cfg["p_fa"]) / np.log(
            sweep["n_data"] / cfg["p_fa"]
        )
        id2p_rows.append({
            "rho": float(rho_v),
            "m_eff": float(sweep["m_eff_vec"][i]),
            "threshold_ratio": float(ratio),
        })

    report_path = write_markdown_report(
        os.path.join(out_dir, "python_cfar_meff_report.md"),
        cfg, id2p_rows, grid_rows,
    )
    if verbose:
        print(f"  Report saved: {report_path}")

    # --- 摘要 ---
    rho_v = sweep['rho_vec']
    me_v  = sweep['m_eff_vec']
    M = sweep['n_data']
    p_fa = cfg['p_fa']

    if verbose:
        print(f'\n=== Summary ===')
        print(f'{"rho":>6} {"M_eff":>8} {"ratio":>8} {"n_null":>7} {"n_hot":>7}')
        for i in range(len(rho_v)):
            ratio = np.log(max(me_v[i], 1e-15)/p_fa) / np.log(M/p_fa)
            print(f'{rho_v[i]:6.2f} {me_v[i]:8.1f} {ratio:8.4f} '
                  f'{sweep["n_null_vec"][i]:7d} {sweep["n_hotspot_vec"][i]:7d}')

        print(f'\nPower allocation @ rho={cfg["rho_high"]:.2f}:')
        for si, snr in enumerate(pa['snr_db_vec']):
            gain = pa['mi_null'][si] - pa['mi_uniform'][si]
            print(f'  SNR={snr:2d} dB: MI uniform={pa["mi_uniform"][si]:.3f}  '
                  f'null-weighted={pa["mi_null"][si]:.3f}  gain={gain:+.3f}')


# ======================================================================
# 参数扫描矩阵
# ======================================================================

def run_python_validation_matrix(base_cfg):
    """
    执行 Python 级验证矩阵扫描。

    扫描维度:
      rho: [0.0, 0.5, 0.9]
      dirichlet_r: [0, 1, 3]

    返回:
      rows: list[dict]
    """
    rows = []
    for rho in [0.0, 0.5, 0.9]:
        for radius in [0, 1, 3]:
            _cfg = dict(base_cfg)
            _cfg["rho"] = rho
            _cfg["dirichlet_r"] = radius
            _cfg["mc_samples"] = int(base_cfg.get("mc_samples_grid", 300))
            result = run_detection_threshold_mc(_cfg)
            result["rho"] = rho
            result["dirichlet_r"] = radius
            rows.append(result)
    return rows


def save_validation_matrix_csv(rows, out_dir):
    """导出扫描矩阵 CSV。"""
    path = os.path.join(out_dir, "python_cfar_validation_matrix.csv")
    header = [
        "rho", "dirichlet_r", "m_raw_grid", "m_eff_grid",
        "threshold_ratio", "empirical_pfa_independent",
        "empirical_pfa_corrected", "empirical_pd_independent",
        "empirical_pd_corrected",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(",".join(header) + "\n")
        for row in rows:
            vals = [str(row.get(key, "")) for key in header]
            handle.write(",".join(vals) + "\n")
    return path


if __name__ == '__main__':
    main()
