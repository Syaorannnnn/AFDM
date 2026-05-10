"""
完整 DAFT 域收发 + 组证据 ΔL_G 对比。

两条路径:
  A. 完美 CSI 基线: H_est=H_eff, 残差=纯噪声, C0=σ²I
  B. OMP 估计:      H_est=OMP(H_eff), 残差含 ID2P, C0=σ²I + H_err*P_data*H_err^H

每条路径对比单候选 vs 邻域组的证据统计量分离比。
本模块不调用 MATLAB 主干。
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from channel import (
    compute_chirp_params, sample_clustered_geometry,
    bdlr_gains, build_path_matrix, compute_effective_channel,
)
from grid_stats import build_candidate_atoms
from detection_mc import _estimate_true_candidate_mask
from cfar import compute_effective_candidates


# ======================================================================
# 公共: 帧构造 + 单径/组证据
# ======================================================================

def build_frame_and_transmit(cfg, rng):
    """
    构造 DAFT 帧 + 真实信道 + 传输。

    返回:
      y, x_tx, sym_data, H_eff, delays, dopplers, gains, noise_power, c1
    """
    n_sc         = int(cfg["n_sc"])
    n_paths      = int(cfg["n_paths"])
    n_clusters   = int(cfg["n_clusters"])
    rho          = float(cfg["rho"])
    max_delay    = int(cfg["max_delay"])
    max_doppler  = int(cfg["max_doppler"])
    doppler_guard = int(cfg["doppler_guard"])
    dirichlet_r  = int(cfg["dirichlet_r"])
    snr_db       = float(cfg["snr_db"])
    pilot_idx    = int(cfg.get("pilot_idx", 0))

    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)

    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    gains = bdlr_gains(n_paths, n_clusters, rho, rng)
    H_eff = compute_effective_channel(delays, dopplers, gains, n_sc, c1, dirichlet_r)

    # 帧: 导频(单元功率) + QPSK 数据(SNR 功率)
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)
    data_amp = np.sqrt(10.0 ** (snr_db / 10.0))

    x_tx = np.zeros(n_sc, dtype=complex)
    x_tx[pilot_idx] = 1.0 + 0j
    bits = rng.integers(0, 2, 2 * n_data)
    sym_data = ((2 * bits[:n_data] - 1) + 1j * (2 * bits[n_data:] - 1)) / np.sqrt(2.0)
    x_tx[data_idx] = data_amp * sym_data

    noise_power = 1.0
    noise = np.sqrt(noise_power / 2.0) * (
        rng.standard_normal(n_sc) + 1j * rng.standard_normal(n_sc))
    y = H_eff @ x_tx + noise

    return y, x_tx, sym_data, H_eff, delays, dopplers, gains, noise_power, c1


def compute_evidence_pair(r_residual, C0, atoms, grid, max_doppler,
                          dirichlet_r, n_sc, loc_step, pilot_idx=0):
    """
    单次残差 + C0 上的单径/组证据。

    返回:
      g_ev, s_ev: (M,) 组/单证据向量
    """
    M = grid.shape[0]
    rad = dirichlet_r
    gs = 2 * rad + 1

    # R_G
    R_G = np.zeros((gs, gs))
    for i in range(gs):
        for j in range(gs):
            R_G[i, j] = np.exp(-2.0 * abs(i - j) / max(rad, 1))
    R_G = 0.9 * R_G + 0.1 * np.eye(gs)

    try:
        C0_inv = np.linalg.inv(C0)
    except np.linalg.LinAlgError:
        reg = 1e-6 * np.trace(C0) / n_sc
        C0_inv = np.linalg.inv(C0 + reg * np.eye(n_sc))

    yb = r_residual.conj() @ C0_inv  # (N,)
    s_ev = np.abs(atoms @ r_residual.conj()) ** 2
    g_ev = np.zeros(M)

    for m in range(M):
        dm = int(grid[m, 0])
        nu_m = int(grid[m, 1])

        Phi = np.zeros((n_sc, gs), dtype=complex)
        vc = 0; cm = []
        for di in range(-rad, rad + 1):
            dk = nu_m + di
            if abs(dk) <= max_doppler:
                tmp = build_path_matrix(dm, float(dk), n_sc, loc_step, dirichlet_r)
                col = tmp[:, pilot_idx]  # 导频列：输入导频→输出全 DAFT 域
                nrm = np.linalg.norm(col)
                if nrm > 1e-12: col /= nrm
                Phi[:, vc] = col
                cm.append(vc)
                vc += 1

        if vc < 2:
            g_ev[m] = s_ev[m]
            continue

        Phi_u = Phi[:, :vc]
        R_u   = R_G[np.ix_(cm, cm)]
        Pb    = Phi_u.conj().T @ C0_inv
        G_m   = Pb @ Phi_u
        M_mat = np.linalg.inv(R_u) + G_m

        try:
            Minv = np.linalg.inv(M_mat)
        except np.linalg.LinAlgError:
            Minv = np.linalg.pinv(M_mat)

        v = yb @ Phi_u
        score   = float(np.real(v @ Minv @ v.conj()))
        _, ld   = np.linalg.slogdet(np.eye(vc) + R_u @ G_m)
        penalty = float(ld)
        g_ev[m] = max(score - penalty, 0.0)

    return g_ev, s_ev


def compute_separation(g_ev, s_ev, grid, delays, dopplers):
    """计算 true/false 分离比。"""
    mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    gt, gf = g_ev[mask], g_ev[~mask]
    st, sf = s_ev[mask], s_ev[~mask]
    eps = 1e-15
    return {
        "sep_group":  float(np.mean(gt) / max(np.mean(gf), eps)),
        "sep_single": float(np.mean(st) / max(np.mean(sf), eps)),
        "g_true_mean": float(np.mean(gt)),
        "s_true_mean": float(np.mean(st)),
        "g_max": float(np.max(g_ev)),
        "s_max": float(np.max(s_ev)),
        "g_peak_at_true": float(gt.max() / max(g_ev.max(), eps)),
    }


# ======================================================================
# 路径 A: 已知路径, 保留 ID2P (r = H_eff@x_data + noise)
# ======================================================================

def run_known_paths_id2p(cfg, rng):
    """
    H_eff 已知, 只移除导频分量, 保留数据干扰。
    残差 r = y - H_eff @ x_pilot = H_eff @ x_data + noise
    这就是 ID2P 的物理定义。
    """
    y, x_tx, sym, H, d, nu, g, N0, c1 = build_frame_and_transmit(cfg, rng)

    n_sc = cfg["n_sc"]
    pilot_idx = cfg.get("pilot_idx", 0)
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])

    # 只移除导频
    x_pilot_only = np.zeros(n_sc, dtype=complex)
    x_pilot_only[pilot_idx] = x_tx[pilot_idx]
    r = y - H @ x_pilot_only  # = H @ x_data + noise

    # C0 = σ²I + H_data @ H_data^H (数据干扰的结构化协方差)
    H_data = H[:, data_idx]
    data_amp = np.sqrt(10.0 ** (cfg["snr_db"] / 10.0))
    C0 = N0 * np.eye(n_sc) + (data_amp ** 2) * np.real(H_data @ H_data.conj().T)

    atoms, grid = build_candidate_atoms(
        n_sc, c1, cfg["dirichlet_r"], cfg["max_delay"], cfg["max_doppler"], pilot_idx)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc

    g_ev, s_ev = compute_evidence_pair(
        r, C0, atoms, grid, cfg["max_doppler"], cfg["dirichlet_r"],
        n_sc, loc_step, pilot_idx)

    sep = compute_separation(g_ev, s_ev, grid, d, nu)
    sep["label"] = "known paths + ID2P"
    return sep


# ======================================================================
# 路径 A': 已知路径 + ID2P, 标量 C0 (对照组)
# ======================================================================

def run_id2p_scalar_c0(cfg, rng):
    """
    与 run_known_paths_id2p 相同的残差, 但用标量 C0 = σ²I。
    用于对比结构化 C0 是否带来组证据改善。
    """
    y, x_tx, sym, H, d, nu, g, N0, c1 = build_frame_and_transmit(cfg, rng)

    n_sc = cfg["n_sc"]
    pilot_idx = cfg.get("pilot_idx", 0)
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])

    x_pilot_only = np.zeros(n_sc, dtype=complex)
    x_pilot_only[pilot_idx] = x_tx[pilot_idx]
    r = y - H @ x_pilot_only

    C0 = N0 * np.eye(n_sc)  # 标量版本

    atoms, grid = build_candidate_atoms(
        n_sc, c1, cfg["dirichlet_r"], cfg["max_delay"], cfg["max_doppler"], pilot_idx)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc

    g_ev, s_ev = compute_evidence_pair(
        r, C0, atoms, grid, cfg["max_doppler"], cfg["dirichlet_r"],
        n_sc, loc_step, pilot_idx)

    sep = compute_separation(g_ev, s_ev, grid, d, nu)
    sep["label"] = "known paths + ID2P (scalar C0)"
    return sep


# ======================================================================
# 路径 B: OMP 估计, 结构化 C0
# ======================================================================

def run_omp_estimate(cfg, rng):
    """OMP 粗估计 → H_est → r + C0_structured。"""
    y, x_tx, sym, H, d, nu, g, N0, c1 = build_frame_and_transmit(cfg, rng)

    n_sc = cfg["n_sc"]
    n_paths = cfg["n_paths"]
    pilot_idx = cfg.get("pilot_idx", 0)
    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)
    data_amp = np.sqrt(10.0 ** (cfg["snr_db"] / 10.0))

    atoms, grid = build_candidate_atoms(
        n_sc, c1, cfg["dirichlet_r"], cfg["max_delay"], cfg["max_doppler"], pilot_idx)
    loc_step = int(round(2.0 * n_sc * c1)) % n_sc
    M = grid.shape[0]

    # OMP
    r_omp = y.copy()
    H_est = np.zeros((n_sc, n_sc), dtype=complex)
    kappa = np.log(M / cfg["p_fa"])
    for _ in range(n_paths):
        scores = np.abs(atoms @ r_omp.conj()) ** 2
        best = int(np.argmax(scores))
        atom_b = atoms[best, :]
        g_est = np.dot(atom_b.conj(), r_omp)
        r_omp -= g_est * atom_b
        db = int(grid[best, 0])
        nb = float(grid[best, 1])
        H_est += g_est * build_path_matrix(db, nb, n_sc, loc_step, cfg["dirichlet_r"])

    true_mask = _estimate_true_candidate_mask(grid, d, nu)
    pd = float(true_mask.sum()) / n_paths

    # 残差 + C0
    r_res = y - H_est @ x_tx
    H_err = H - H_est
    C0_data = (data_amp ** 2) * np.real(H_err @ H_err.conj().T)
    C0 = N0 * np.eye(n_sc) + C0_data

    g_ev, s_ev = compute_evidence_pair(
        r_res, C0, atoms, grid, cfg["max_doppler"], cfg["dirichlet_r"],
        n_sc, loc_step, pilot_idx)

    sep = compute_separation(g_ev, s_ev, grid, d, nu)
    sep["label"] = "OMP + structured C0"
    sep["pd"] = pd
    sep["residual_power"] = float(np.mean(np.abs(r_res) ** 2))
    return sep


# ======================================================================
# 批量对比
# ======================================================================

def compare_paths(cfg, n_trials=5):
    """
    运行 n_trials 次，对比路径 A vs B 的组/单证据分离比。
    """
    base = dict(cfg)
    results = {"id2p_structured": [], "id2p_scalar": [], "omp": []}

    for t in range(n_trials):
        _c = dict(base)
        _c["seed"] = base.get("seed", 42) + t

        # A: 已知路径 + ID2P, 结构化 C0
        rng = np.random.default_rng(_c["seed"])
        pa = run_known_paths_id2p(_c, rng)
        results["id2p_structured"].append(pa)

        # A': 已知路径 + ID2P, 标量 C0 (对照组)
        rng2 = np.random.default_rng(_c["seed"])
        pc = run_id2p_scalar_c0(_c, rng2)
        results["id2p_scalar"].append(pc)

        # B: OMP + 结构化 C0
        rng3 = np.random.default_rng(_c["seed"])
        pb = run_omp_estimate(_c, rng3)
        results["omp"].append(pb)

    # 汇总
    for key in ["id2p_structured", "id2p_scalar", "omp"]:
        rows = results[key]
        n = len(rows)
        sep_g = np.mean([r["sep_group"] for r in rows])
        sep_s = np.mean([r["sep_single"] for r in rows])
        label = rows[0]["label"]
        extra = ""
        if key == "omp":
            pd_m = np.mean([r["pd"] for r in rows])
            rp_m = np.mean([r["residual_power"] for r in rows])
            extra = f"  Pd={pd_m:.2f}  res_pow={rp_m:.1f}"
        print(f"  {label:20s}: sep_group={sep_g:7.1f}  sep_single={sep_s:7.1f}  "
              f"gain={sep_g/max(sep_s,1e-15):.1f}x{extra}")

    return results


if __name__ == "__main__":
    base = {
        "n_sc": 64, "n_paths": 6, "n_clusters": 2,
        "max_delay": 3, "max_doppler": 2,
        "doppler_guard": 3, "dirichlet_r": 3,
        "rho": 0.5, "snr_db": 15.0, "p_fa": 1e-3,
        "seed": 42, "pilot_idx": 0,
    }
    print("=== 全帧 ΔL_G 对比: 完美CSI vs OMP+结构化C0 ===\n")
    compare_paths(base, n_trials=5)
