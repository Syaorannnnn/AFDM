"""
CLIP systematic debug — 逐层打印中间量，定位 NMSE/BER 高居不下的根因。

检查顺序:
  D1. 信道构造正确性 (H_true vs build_path_matrix)
  D2. 增益估计误差 (g_est vs g_true, 在纯导频信号上)
  D3. 增益估计误差 (g_est vs g_true, 在完整帧上)
  D4. H_est 质量 (NMSE 分解: 真径误差 vs 假径污染)
  D5. LMMSE 输入质量 (H_est 列空间与 H_true 列空间对齐度)
  D6. ID2P 消除效果 (消除前后残差功率)
  D7. 数据检测 SNR (等化后每符号 SNR)
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
from group_omp import build_group_atoms


# ======================================================================
# 固定 seed 的单次 trial 环境
# ======================================================================

def build_env(snr_db=15.0, rho=0.5, seed=42):
    n_sc, n_paths, n_clusters = 64, 6, 2
    max_delay, max_doppler, doppler_guard, dirichlet_r = 3, 2, 3, 3
    pilot_idx = 0
    pilot_amp = np.sqrt(10.0 ** 3.5)   # 35 dB
    data_amp  = np.sqrt(10.0 ** (snr_db / 10.0))

    rng = np.random.default_rng(seed)
    c1, loc_step = compute_chirp_params(n_sc, max_doppler, doppler_guard)

    delays, dopplers = sample_clustered_geometry(
        n_paths, n_clusters, max_delay, max_doppler, rng)
    gains_true = bdlr_gains(n_paths, n_clusters, rho, rng)
    H_true = compute_effective_channel(delays, dopplers, gains_true, n_sc, c1, dirichlet_r)

    data_idx = np.array([j for j in range(n_sc) if j != pilot_idx])
    n_data = len(data_idx)

    bits = rng.integers(0, 2, 2 * n_data)
    sym_data = ((2*bits[:n_data]-1) + 1j*(2*bits[n_data:]-1)) / np.sqrt(2.0)

    x_tx = np.zeros(n_sc, dtype=complex)
    x_tx[pilot_idx] = pilot_amp
    x_tx[data_idx]  = data_amp * sym_data

    noise = (rng.standard_normal(n_sc) + 1j*rng.standard_normal(n_sc)) / np.sqrt(2.0)
    y = H_true @ x_tx + noise

    atoms, grid = build_candidate_atoms(n_sc, c1, dirichlet_r, max_delay, max_doppler, pilot_idx)
    group_atoms, group_masks = build_group_atoms(atoms, grid, max_doppler, dirichlet_r)
    true_mask = _estimate_true_candidate_mask(grid, delays, dopplers)
    true_set  = set(np.where(true_mask)[0])

    return dict(
        n_sc=n_sc, n_paths=n_paths, pilot_idx=pilot_idx,
        pilot_amp=pilot_amp, data_amp=data_amp,
        c1=c1, loc_step=loc_step, dirichlet_r=dirichlet_r,
        delays=delays, dopplers=dopplers, gains_true=gains_true,
        H_true=H_true, x_tx=x_tx, sym_data=sym_data,
        data_idx=data_idx, n_data=n_data, bits=bits,
        y=y, atoms=atoms, grid=grid,
        group_atoms=group_atoms, group_masks=group_masks,
        true_set=true_set,
    )


def sep(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print('='*60)


# ======================================================================
# D1: 信道构造正确性
# ======================================================================

def d1_channel_construction(env):
    sep("D1: 信道构造正确性")
    H_true = env["H_true"]
    gains_true = env["gains_true"]
    delays = env["delays"]
    dopplers = env["dopplers"]
    c1 = env["c1"]
    loc_step = env["loc_step"]
    n_sc = env["n_sc"]
    dirichlet_r = env["dirichlet_r"]

    # 手动重建 H_true
    H_rebuild = np.zeros((n_sc, n_sc), dtype=complex)
    for i in range(env["n_paths"]):
        H_rebuild += gains_true[i] * build_path_matrix(
            delays[i], dopplers[i], n_sc, loc_step, dirichlet_r)

    err = np.linalg.norm(H_true - H_rebuild) / np.linalg.norm(H_true)
    print(f"  H_true 重建相对误差: {err:.2e}  (应 < 1e-10)")

    # 信道能量分布
    print(f"  ||H_true||_F²  = {np.sum(np.abs(H_true)**2):.4f}")
    print(f"  ||H_true[:,0]||² (pilot列) = {np.sum(np.abs(H_true[:,0])**2):.4f}")
    print(f"  真实增益 |g_i|: {np.abs(gains_true)}")
    print(f"  真实延迟/Doppler: {list(zip(delays, dopplers))}")


# ======================================================================
# D2: 纯导频信号上的增益估计
# ======================================================================

def d2_gain_estimation_pilot_only(env):
    sep("D2: 纯导频信号上的增益估计 (无数据干扰)")
    H_true = env["H_true"]
    pilot_amp = env["pilot_amp"]
    n_sc = env["n_sc"]
    atoms = env["atoms"]
    grid = env["grid"]
    true_set = env["true_set"]
    gains_true = env["gains_true"]
    c1 = env["c1"]
    loc_step = env["loc_step"]
    dirichlet_r = env["dirichlet_r"]

    # 纯导频发射
    x_pilot = np.zeros(n_sc, dtype=complex)
    x_pilot[0] = pilot_amp
    y_pilot = H_true @ x_pilot  # 无噪声, 无数据

    # 在真实候选上估计增益
    true_list = sorted(true_set)
    r = y_pilot.copy()
    pilot_amp_est = np.sqrt(10.0 ** 3.5)
    gains_est = {}
    for idx in true_list:
        atom = atoms[idx, :]
        gi = np.dot(atom.conj(), r) / pilot_amp_est
        r -= gi * pilot_amp_est * atom
        gains_est[idx] = gi

    print(f"  真实候选索引: {true_list}")
    print(f"  真实增益 |g_true|: {np.abs(gains_true)}")
    print(f"  估计增益 |g_est|:  {[abs(gains_est.get(i, 0)) for i in true_list]}")

    # 重建 H_est
    H_est = np.zeros((n_sc, n_sc), dtype=complex)
    for idx in true_list:
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        H_est += gains_est[idx] * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)

    nmse = 10*np.log10(np.sum(np.abs(H_true - H_est)**2) /
                       max(np.sum(np.abs(H_true)**2), 1e-15))
    print(f"  NMSE (纯导频, 真实支撑): {nmse:+.2f} dB  (应 < -20 dB)")
    return gains_est, H_est


# ======================================================================
# D3: 完整帧上的增益估计 (含数据干扰)
# ======================================================================

def d3_gain_estimation_full_frame(env):
    sep("D3: 完整帧上的增益估计 (含数据干扰)")
    H_true = env["H_true"]
    y = env["y"]
    atoms = env["atoms"]
    grid = env["grid"]
    true_set = env["true_set"]
    gains_true = env["gains_true"]
    c1 = env["c1"]
    loc_step = env["loc_step"]
    dirichlet_r = env["dirichlet_r"]
    n_sc = env["n_sc"]

    true_list = sorted(true_set)
    r = y.copy()
    pilot_amp_est = np.sqrt(10.0 ** 3.5)
    gains_est = {}
    for idx in true_list:
        atom = atoms[idx, :]
        gi = np.dot(atom.conj(), r) / pilot_amp_est
        r -= gi * pilot_amp_est * atom
        gains_est[idx] = gi

    print(f"  真实增益 |g_true|: {np.abs(gains_true)}")
    print(f"  估计增益 |g_est|:  {[abs(gains_est.get(i, 0)) for i in true_list]}")

    # 增益误差
    for i, idx in enumerate(true_list):
        g_t = gains_true[i]
        g_e = gains_est.get(idx, 0)
        err = abs(g_e - g_t) / max(abs(g_t), 1e-15)
        print(f"    path {i}: |g_true|={abs(g_t):.4f}  |g_est|={abs(g_e):.4f}  "
              f"phase_err={np.angle(g_e/g_t)*180/np.pi:.1f}°  rel_err={err:.2%}")

    H_est = np.zeros((n_sc, n_sc), dtype=complex)
    for idx in true_list:
        d_est = int(grid[idx, 0])
        nu_est = float(grid[idx, 1])
        H_est += gains_est[idx] * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)

    nmse = 10*np.log10(np.sum(np.abs(H_true - H_est)**2) /
                       max(np.sum(np.abs(H_true)**2), 1e-15))
    print(f"  NMSE (完整帧, 真实支撑): {nmse:+.2f} dB")
    return gains_est, H_est


# ======================================================================
# D4: NMSE 分解 — 真径误差 vs 假径污染
# ======================================================================

def d4_nmse_decomposition(env, detected_paths):
    sep("D4: NMSE 分解 (真径误差 vs 假径污染)")
    H_true = env["H_true"]
    y = env["y"]
    atoms = env["atoms"]
    grid = env["grid"]
    true_set = env["true_set"]
    c1 = env["c1"]
    loc_step = env["loc_step"]
    dirichlet_r = env["dirichlet_r"]
    n_sc = env["n_sc"]

    true_in_det  = [i for i in detected_paths if i in true_set]
    false_in_det = [i for i in detected_paths if i not in true_set]
    print(f"  检测路径总数: {len(detected_paths)}")
    print(f"  真径: {len(true_in_det)} / {len(true_set)}  假径: {len(false_in_det)}")

    pilot_amp_est = np.sqrt(10.0 ** 3.5)

    def build_H_subset(indices):
        r = y.copy()
        gains = {}
        for idx in indices:
            atom = atoms[idx, :]
            gi = np.dot(atom.conj(), r) / pilot_amp_est
            r -= gi * pilot_amp_est * atom
            gains[idx] = gi
        H = np.zeros((n_sc, n_sc), dtype=complex)
        for idx in indices:
            d_est = int(grid[idx, 0])
            nu_est = float(grid[idx, 1])
            H += gains[idx] * build_path_matrix(d_est, nu_est, n_sc, loc_step, dirichlet_r)
        return H

    H_true_part  = build_H_subset(true_in_det)
    H_false_part = build_H_subset(false_in_det)
    H_full       = build_H_subset(detected_paths)

    norm_true = np.sum(np.abs(H_true)**2)
    nmse_true_only = 10*np.log10(np.sum(np.abs(H_true - H_true_part)**2) / norm_true)
    nmse_false_cont = 10*np.log10(np.sum(np.abs(H_false_part)**2) / norm_true)
    nmse_full = 10*np.log10(np.sum(np.abs(H_true - H_full)**2) / norm_true)

    print(f"  NMSE (仅真径部分):   {nmse_true_only:+.2f} dB  ← 真径估计误差")
    print(f"  假径能量 / ||H_true||²: {nmse_false_cont:+.2f} dB  ← 假径污染")
    print(f"  NMSE (全检测集):     {nmse_full:+.2f} dB")
    return H_full


# ======================================================================
# D5: LMMSE 输入质量
# ======================================================================

def d5_lmmse_input_quality(env, H_est):
    sep("D5: LMMSE 输入质量")
    H_true = env["H_true"]
    data_idx = env["data_idx"]
    n_data = env["n_data"]
    n_sc = env["n_sc"]

    # 列空间对齐度: 每列的余弦相似度
    cos_sims = []
    for j in range(n_sc):
        h_t = H_true[:, j]
        h_e = H_est[:, j]
        nt, ne = np.linalg.norm(h_t), np.linalg.norm(h_e)
        if nt > 1e-12 and ne > 1e-12:
            cos_sims.append(abs(np.dot(h_t.conj(), h_e)) / (nt * ne))

    print(f"  列余弦相似度: mean={np.mean(cos_sims):.4f}  "
          f"min={np.min(cos_sims):.4f}  (1.0 = 完美对齐)")

    # 数据列的条件数
    H_data = H_est[:, data_idx]
    sv = np.linalg.svd(H_data, compute_uv=False)
    print(f"  H_est 数据列奇异值: max={sv[0]:.4f}  min={sv[-1]:.6f}  "
          f"cond={sv[0]/max(sv[-1],1e-15):.1f}")

    # 有效 SNR: ||H_true @ x_data||² / noise
    y = env["y"]
    x_tx = env["x_tx"]
    pilot_idx = env["pilot_idx"]
    x_data_only = x_tx.copy()
    x_data_only[pilot_idx] = 0
    id2p_power = np.mean(np.abs(H_true @ x_data_only)**2)
    noise_power = 1.0
    print(f"  ID2P 功率 / 噪声功率: {id2p_power:.2f}  "
          f"({10*np.log10(id2p_power):.1f} dB 高于噪声底)")


# ======================================================================
# D6: ID2P 消除效果
# ======================================================================

def d6_id2p_cancellation(env, H_est):
    sep("D6: ID2P 消除效果")
    y = env["y"]
    x_tx = env["x_tx"]
    H_true = env["H_true"]
    data_idx = env["data_idx"]
    n_data = env["n_data"]
    n_sc = env["n_sc"]
    pilot_idx = env["pilot_idx"]
    noise_power = 1.0

    # 理想消除 (用 H_true + 真实数据)
    x_data_true = x_tx.copy()
    x_data_true[pilot_idx] = 0
    y_ideal_clean = y - H_true @ x_data_true
    print(f"  理想消除后残差功率: {np.mean(np.abs(y_ideal_clean)**2):.4f}  "
          f"(应 ≈ noise_power={noise_power})")

    # 用 H_est + 硬判决消除
    H_data = H_est[:, data_idx]
    reg = noise_power
    HhH = H_data.conj().T @ H_data
    w_mmse = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
    x_soft = w_mmse @ y
    x_hard = ((np.real(x_soft) > 0).astype(float)*2-1 +
              1j*(np.imag(x_soft) > 0).astype(float)*2-1j) / np.sqrt(2.0)

    # 硬判决正确率
    sym_data = env["sym_data"]
    data_amp = env["data_amp"]
    x_hard_true = ((np.real(sym_data) > 0).astype(float)*2-1 +
                   1j*(np.imag(sym_data) > 0).astype(float)*2-1j) / np.sqrt(2.0)
    correct_rate = float(np.mean(x_hard == x_hard_true))
    print(f"  LMMSE 硬判决正确率: {correct_rate:.3f}  (0.5=随机, 1.0=完美)")

    x_fb = np.zeros(n_sc, dtype=complex)
    x_fb[data_idx] = x_hard * data_amp
    id2p_est = H_est @ x_fb
    y_clean = y - id2p_est
    print(f"  H_est 消除后残差功率: {np.mean(np.abs(y_clean)**2):.4f}")

    # 残差中真实导频响应的保留程度
    x_pilot_only = np.zeros(n_sc, dtype=complex)
    x_pilot_only[pilot_idx] = env["pilot_amp"]
    pilot_response = H_true @ x_pilot_only
    cos_sim = abs(np.dot(pilot_response.conj(), y_clean)) / (
        np.linalg.norm(pilot_response) * np.linalg.norm(y_clean) + 1e-15)
    print(f"  残差与真实导频响应余弦相似度: {cos_sim:.4f}  (1.0 = 消除完美)")


# ======================================================================
# D7: 等化后每符号 SNR
# ======================================================================

def d7_post_equalization_snr(env, H_est):
    sep("D7: 等化后每符号 SNR")
    y = env["y"]
    H_true = env["H_true"]
    data_idx = env["data_idx"]
    n_data = env["n_data"]
    sym_data = env["sym_data"]
    data_amp = env["data_amp"]
    noise_power = 1.0

    def eval_ber_snr(H, label):
        H_data = H[:, data_idx]
        reg = noise_power
        HhH = H_data.conj().T @ H_data
        w = np.linalg.solve(HhH + reg * np.eye(n_data), H_data.conj().T)
        x_est = w @ y
        # 信号分量: w @ H_true @ x_data_true
        x_data_true = np.zeros(len(y), dtype=complex)
        x_data_true[data_idx] = data_amp * sym_data
        sig = w @ (H_true @ x_data_true)
        noise_out = x_est - sig
        snr_per_sym = np.mean(np.abs(sig)**2) / max(np.mean(np.abs(noise_out)**2), 1e-15)
        ber = float(np.mean((np.real(x_est) > 0).astype(int) !=
                            (np.real(sym_data) > 0).astype(int)))
        print(f"  [{label}] 等化后 SNR: {10*np.log10(snr_per_sym):.1f} dB  BER: {ber:.4f}")

    eval_ber_snr(H_true, "H_true (理想)")
    eval_ber_snr(H_est,  "H_est  (估计)")


# ======================================================================
# 主流程
# ======================================================================

if __name__ == "__main__":
    print("\n" + "="*60)
    print("  CLIP Systematic Debug")
    print("="*60)

    for snr_db in [15.0, 25.0]:
        print(f"\n\n{'#'*60}")
        print(f"  SNR = {snr_db} dB,  rho = 0.5")
        print(f"{'#'*60}")

        env = build_env(snr_db=snr_db, rho=0.5, seed=42)

        d1_channel_construction(env)
        _, H_pilot = d2_gain_estimation_pilot_only(env)
        _, H_full  = d3_gain_estimation_full_frame(env)

        # 模拟 Pd=0.33 的检测结果: 2 真径 + 4 假径
        true_list  = sorted(env["true_set"])
        M = env["grid"].shape[0]
        rng2 = np.random.default_rng(99)
        false_pool = [i for i in range(M) if i not in env["true_set"]]
        false_sample = list(rng2.choice(false_pool, size=4, replace=False))
        detected_2t4f = true_list[:2] + false_sample
        detected_all_true = true_list  # 理想情况

        print("\n  --- 场景 A: 2 真径 + 4 假径 (Pd=0.33) ---")
        H_2t4f = d4_nmse_decomposition(env, detected_2t4f)
        d5_lmmse_input_quality(env, H_2t4f)
        d6_id2p_cancellation(env, H_2t4f)
        d7_post_equalization_snr(env, H_2t4f)

        print("\n  --- 场景 B: 全真径 (Pd=1.0, 理想上界) ---")
        H_ideal = d4_nmse_decomposition(env, detected_all_true)
        d5_lmmse_input_quality(env, H_ideal)
        d6_id2p_cancellation(env, H_ideal)
        d7_post_equalization_snr(env, H_ideal)

        print("\n  --- 场景 C: 纯导频增益估计 (D2 结果) ---")
        d7_post_equalization_snr(env, H_pilot)
