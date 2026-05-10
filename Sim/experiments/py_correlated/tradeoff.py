#!/usr/bin/env python3
"""
κ-Pfa-Pd 权衡曲线 + 组检测 ΔL_G 原型。

A2: 扫描 κ 乘数，生成 ROC 风格 Pfa-Pd 权衡曲线。
A3: 将候选网格点按簇分组，用组证据统计量替代单候选判定。

用法:
  python tradeoff.py --mc 10000
  python tradeoff.py --mc 10000 --rho 0.5
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _plot_setup import setup_plots
setup_plots()
import matplotlib.pyplot as plt

from channel import (
    compute_chirp_params, sample_clustered_geometry,
    bdlr_gains, build_path_matrix,
)
from grid_stats import build_candidate_atoms
from cfar import compute_effective_candidates
from detection_mc import _estimate_true_candidate_mask


# ======================================================================
# A2: κ 扫描权衡曲线
# ======================================================================

def run_kappa_tradeoff(cfg):
    """
    扫描 κ 乘数，计算每个 κ 下的 Pfa 和 Pd。

    返回:
      kappa_vec: (n_scan,) 扫描的 κ 值
      pfa_vec:   (n_scan,) 每个 κ 的经验 Pfa
      pd_vec:    (n_scan,) 每个 κ 的经验 Pd
      kappa_ind: float ln(M/P_fa)
    """
    from grid_stats import estimate_grid_stat_covariance

    p_fa_target = float(cfg.get("p_fa", 1e-3))
    mc_samples  = int(cfg.get("mc_samples", 10000))

    stats = estimate_grid_stat_covariance(cfg)
    power_samples = stats["power_samples"]
    grid = stats["grid"]
    true_mask = _estimate_true_candidate_mask(
        grid, stats["delays"], stats["dopplers"])
    false_mask = ~true_mask

    false_power = float(np.mean(power_samples[:, false_mask]))
    residual_power = false_power
    M = float(grid.shape[0])
    kappa_ind = np.log(M / p_fa_target)

    false_stats = power_samples[:, false_mask] / residual_power
    true_stats  = power_samples[:, true_mask] / residual_power

    z_false = false_stats.ravel()
    z_true  = true_stats.ravel()

    # 在 κ_ind 的 0.3x ~ 5x 范围内取 80 个点
    kappa_vec = np.linspace(kappa_ind * 0.3, kappa_ind * 5.0, 80)
    pfa_vec = np.array([float(np.mean(z_false > k)) for k in kappa_vec])
    pd_vec  = np.array([float(np.mean(z_true  > k)) for k in kappa_vec])

    return kappa_vec, pfa_vec, pd_vec, float(kappa_ind), float(M)


def plot_tradeoff(results, out_dir):
    """Fig5f: κ-Pfa-Pd 权衡曲线。"""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # 场景标签
    labels = {
        (0.0, 3): 'rho=0.0, dir_r=3',
        (0.5, 3): 'rho=0.5, dir_r=3',
        (0.9, 3): 'rho=0.9, dir_r=3',
        (0.0, 1): 'rho=0.0, dir_r=1',
        (0.9, 1): 'rho=0.9, dir_r=1',
    }

    colors = plt.cm.tab10(np.linspace(0, 1, len(results)))

    # 左: Pfa vs κ
    ax = axes[0]
    for (key, (kvec, pfa, pd, kind, M)), c in zip(results.items(), colors):
        ax.semilogy(kvec, pfa, color=c, linewidth=1.5,
                     label=labels.get(key, str(key)))
    ax.axhline(1e-3, color='gray', linestyle='--', alpha=0.5,
               label=r'design $P_{fa}=10^{-3}$')
    ax.axhline(1e-2, color='gray', linestyle=':', alpha=0.3)
    ax.set_xlabel(r'$\kappa$')
    ax.set_ylabel(r'$P_{fa}$')
    ax.set_title(r'$P_{fa}$ vs $\kappa$')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 中: Pd vs κ
    ax = axes[1]
    for (key, (kvec, pfa, pd, kind, M)), c in zip(results.items(), colors):
        ax.plot(kvec, pd, color=c, linewidth=1.5,
                label=labels.get(key, str(key)))
    ax.set_xlabel(r'$\kappa$')
    ax.set_ylabel(r'$P_d$')
    ax.set_title(r'$P_d$ vs $\kappa$')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 右: ROC 风格 Pfa-Pd
    ax = axes[2]
    for (key, (kvec, pfa, pd, kind, M)), c in zip(results.items(), colors):
        ax.semilogx(pfa, pd, color=c, linewidth=1.5,
                     label=labels.get(key, str(key)))
    ax.axvline(1e-3, color='gray', linestyle='--', alpha=0.5,
               label=r'$P_{fa}=10^{-3}$')
    ax.set_xlabel(r'$P_{fa}$')
    ax.set_ylabel(r'$P_d$')
    ax.set_title(r'Detection trade-off (ROC style)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(r'CFAR $\kappa$-Pfa-Pd Trade-off under Fractional Doppler',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, 'Fig5f_KappaTradeoff.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print(f'  Fig5f saved: {path}')


# ======================================================================
# A3: 组检测 ΔL_G 原型
# ======================================================================

def compute_group_evidence(delays, dopplers, cluster_labels, n_sc, c1,
                           dirichlet_r, pilot_idx, data_idx, n_clusters,
                           rho, gain_var, rng, max_delay, max_doppler):
    """
    组证据统计量 ΔL_G：按簇分组后计算联合边际似然增益。

    对每个簇 c（包含 Pc 条路径），构造组观测矩阵 Phi_c ∈ C^{N × N×Pc}，
    计算:
      S_c = R_c ⊗ (Phi_c^H C0^{-1} Phi_c)
      ΔL_c = y_c^H (Phi_c / C0) y_c 的修正形式

    简化版: 直接用 pilot 行响应，比较组内相干叠加 vs 非相干叠加的增益比。

    返回:
      group_scores:    (n_clusters, n_candidates) 组级得分
      single_scores:   (n_paths, n_candidates) 单径得分
    """
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc

    # 候选网格原子 - 使用配置中的 max_delay/max_doppler，确保与后续 grid 一致
    atoms, _ = build_candidate_atoms(
        n_sc, c1, dirichlet_r, max_delay, max_doppler,
        pilot_idx
    )

    # 真实路径的 pilot 响应
    n_paths = len(delays)
    pc = n_paths // n_clusters
    path_mats = [
        build_path_matrix(int(delays[i]), dopplers[i], n_sc, loc_step,
                          dirichlet_r)
        for i in range(n_paths)
    ]
    gains = bdlr_gains(n_paths, n_clusters, rho, rng)
    pilot_row = np.zeros(n_sc, dtype=complex)
    for path_idx, phi in enumerate(path_mats):
        pilot_row += gains[path_idx] * phi[pilot_idx, :]

    # 单径匹配得分 (归一化原子内积的模方)
    single_scores = np.abs(atoms @ pilot_row.conj()) ** 2

    # 组匹配: 同簇内所有路径的原子行做相干叠加
    group_scores = np.zeros((n_clusters, atoms.shape[0]))
    for c in range(n_clusters):
        lo = c * pc
        hi = lo + pc
        # 簇内路径的联合原子: 取每径 atom 行的平均方向
        # 简化: 用簇内所有径的 (delay, doppler) 中心作为组代表
        c_delays = delays[lo:hi]
        c_dopp = dopplers[lo:hi]
        base_d = int(np.median(c_delays))
        base_nu = float(np.median(c_dopp))
        # 组原子: 以簇中心构建
        phi_group = build_path_matrix(base_d, base_nu, n_sc, loc_step,
                                      dirichlet_r)
        group_atom = phi_group[pilot_idx, :]
        norm_val = np.linalg.norm(group_atom)
        if norm_val > 1e-12:
            group_atom = group_atom / norm_val
        group_scores[c, :] = np.abs(atoms @ group_atom.conj()) ** 2

    return group_scores, single_scores, gains, pilot_row, atoms


def compare_group_vs_single(cfg):
    """
    对比组检测 vs 单径检测的匹配得分提升。

    返回:
      dict 含 single_score_true (真实候选得分), group_score_true,
           single_score_false, group_score_false 的对比
    """
    n_sc         = int(cfg.get("n_sc", 64))
    n_paths      = int(cfg.get("n_paths", 6))
    n_clusters   = int(cfg.get("n_clusters", 2))
    max_delay    = int(cfg.get("max_delay", 3))
    max_doppler  = int(cfg.get("max_doppler", 2))
    doppler_guard = int(cfg.get("doppler_guard", 3))
    dirichlet_r  = int(cfg.get("dirichlet_r", 3))
    rho          = float(cfg.get("rho", 0.5))
    seed         = int(cfg.get("seed", 42))
    pilot_idx    = int(cfg.get("pilot_idx", 0))

    rng = np.random.default_rng(seed)
    c1, _ = compute_chirp_params(n_sc, max_doppler, doppler_guard)
    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng
    )
    cluster_labels = np.repeat(np.arange(n_clusters), n_paths // n_clusters)

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    gain_var = 1.0 / n_paths

    group_scores, single_scores, gains, pilot_row, atoms = \
        compute_group_evidence(
            delays, dopplers, cluster_labels, n_sc, c1,
            dirichlet_r, pilot_idx, data_idx, n_clusters, rho, gain_var, rng,
            max_delay, max_doppler
        )

    _, grid = build_candidate_atoms(
        n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx
    )
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    false_mask = ~true_mask

    # 单径得分
    ss_true  = single_scores[true_mask]
    ss_false = single_scores[false_mask]

    # 组得分（取所有簇的最大值）
    gs = group_scores.max(axis=0)
    gs_true  = gs[true_mask]
    gs_false = gs[false_mask]

    return {
        "single_true_mean": float(np.mean(ss_true)),
        "single_false_mean": float(np.mean(ss_false)),
        "group_true_mean": float(np.mean(gs_true)),
        "group_false_mean": float(np.mean(gs_false)),
        "single_true_max": float(np.max(ss_true)),
        "group_true_max": float(np.max(gs_true)),
        "separation_single": float(np.mean(ss_true) / max(np.mean(ss_false), 1e-15)),
        "separation_group": float(np.mean(gs_true) / max(np.mean(gs_false), 1e-15)),
        "n_true": int(true_mask.sum()),
        "n_false": int(false_mask.sum()),
    }


# ======================================================================
# main
# ======================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", type=int, default=10000)
    ap.add_argument("--rho", type=float, default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=str, default=None)
    ap.add_argument("--group", action="store_true",
                    help="运行组检测对比 (A3)")
    args = ap.parse_args()

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

    # A2: κ 权衡曲线
    print(f"=== A2: kappa-Pfa-Pd trade-off (mc={args.mc}) ===")
    # 扫描 dir_r=1,3 和 rho=0.0,0.9
    configs = {}
    for dr in [1, 3]:
        for rho_v in [0.0, 0.9]:
            _c = dict(base)
            _c["dirichlet_r"] = dr
            _c["rho"] = rho_v
            key = (rho_v, dr)
            configs[key] = _c

    if args.rho is not None:
        _c = dict(base)
        _c["rho"] = args.rho
        configs = {(args.rho, 3): _c}

    results = {}
    for key, cfg in configs.items():
        print(f"  rho={key[0]:.1f}, dir_r={key[1]}")
        kvec, pfa, pd_vec, kind, M = run_kappa_tradeoff(cfg)
        results[key] = (kvec, pfa, pd_vec, kind, M)
        # 输出关键点
        idx_1e3 = np.argmin(np.abs(pfa - 1e-3))
        idx_1e2 = np.argmin(np.abs(pfa - 1e-2))
        print(f"    at Pfa=1e-3: kappa={kvec[idx_1e3]:.1f}, Pd={pd_vec[idx_1e3]:.4f}")
        print(f"    at Pfa=1e-2: kappa={kvec[idx_1e2]:.1f}, Pd={pd_vec[idx_1e2]:.4f}")

    plot_tradeoff(results, out_dir)

    # A3: 组检测对比
    if args.group:
        print(f"\n=== A3: Group detection vs single-path detection ===")
        for dr in [1, 3]:
            for rho_v in [0.0, 0.5, 0.9]:
                _c = dict(base)
                _c["dirichlet_r"] = dr
                _c["rho"] = rho_v
                _c["mc_samples"] = min(args.mc, 5000)
                cmp = compare_group_vs_single(_c)
                print(f"  rho={rho_v:.1f} dir_r={dr}: "
                      f"single_sep={cmp['separation_single']:.2f} "
                      f"group_sep={cmp['separation_group']:.2f} "
                      f"| ss_true_mean={cmp['single_true_mean']:.4f} "
                      f"gs_true_mean={cmp['group_true_mean']:.4f}")


if __name__ == '__main__':
    main()
