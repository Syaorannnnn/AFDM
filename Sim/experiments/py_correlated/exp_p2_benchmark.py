"""P2 benchmark sweep runner.

Grid: rho x SNR x cluster_profile x detector, n_trials per cell.
Writes:
  Sim/experiments/Results/p2/<timestamp>/summary.json
  Sim/experiments/Results/p2/<timestamp>/nll_<cell_id>.csv  (for lrt detector cells)

Baselines:
  detector="naive": naive matched filter T=|Phi_k^H r|^2 with scalar CFAR
  detector="lrt":   full clip_lrt.clip_lrt_receive

Metrics (per cell, averaged over trials):
  - BER (QPSK hard-decision over reconstructed data subcarriers)
  - Cluster Pd = |support_hat ∩ support_true| / |support_true|
  - NMSE (dB) of recovered Sigma vs ground-truth Sigma
"""

import json
import time
from pathlib import Path

import numpy as np

from clip_lrt import ClipLrtConfig, clip_lrt_receive
from sigma_mom import sigma_from_gains


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def _simulate_observation(Phi_per_cluster, g_true, N0, rng):
    N = Phi_per_cluster[0].shape[0]
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    y = np.zeros(N, dtype=complex)
    for c, g in enumerate(g_true):
        h = g * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2))
        y = y + h * (Phi_per_cluster[c] @ x_unit)
    y = y + (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
    return y


def _naive_receive(y, Phi_per_cluster, N0):
    """Scalar-Sigma baseline: pick the single cluster with largest |Phi_c^H y|^2."""
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    stats = np.zeros(C)
    gains = np.zeros(C, dtype=complex)
    for c in range(C):
        v_c = Phi_per_cluster[c] @ x_unit
        proj = np.vdot(v_c, y)
        gains[c] = proj / max(float(np.vdot(v_c, v_c).real), 1e-16)
        stats[c] = abs(proj) ** 2
    # Simple naive support: top-ceil(C/2) clusters
    k_keep = max(1, C // 2)
    support = np.zeros(C, dtype=bool)
    support[np.argsort(-stats)[:k_keep]] = True
    gains = gains * support
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    return {"gains": gains, "support": support, "sigma": Sigma}


def _bdlr_rho_from_profile(profile, rho):
    """Each cluster gets the same rho."""
    return [rho] * len(profile)


def _eval_cell(rho, snr_db, profile, detector, n_trials, N, seed):
    rng = np.random.default_rng(seed)
    C = len(profile)
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    rho_list = _bdlr_rho_from_profile(profile, rho)
    # assign ground-truth per-cluster gain amplitudes from profile weights
    total = float(sum(profile))
    g_true = np.array([np.sqrt(pc / total) * np.sqrt(rho_list[c]) for c, pc in enumerate(profile)], dtype=complex)
    N0 = 10.0 ** (-snr_db / 10.0)

    ber_list = []
    pd_list = []
    nmse_list = []
    for t in range(n_trials):
        rng_trial = np.random.default_rng(seed * 997 + t)
        y = _simulate_observation(Phi_per_cluster, g_true, N0, rng_trial)
        if detector == "naive":
            out = _naive_receive(y, Phi_per_cluster, N0)
        else:
            cfg = ClipLrtConfig(Phi_per_cluster=Phi_per_cluster, N0=N0)
            out = clip_lrt_receive(y, pilot_frame=None, config=cfg)
        # support-based Cluster Pd
        support_true = np.abs(g_true) > 1e-9
        hit = int(np.sum(out["support"] & support_true))
        total_true = int(np.sum(support_true))
        pd = hit / max(total_true, 1)
        pd_list.append(pd)
        # NMSE of Sigma recovery
        Sigma_true = sigma_from_gains(g_true, Phi_per_cluster, N0)
        num = float(np.linalg.norm(out["sigma"] - Sigma_true, ord="fro") ** 2)
        den = max(float(np.linalg.norm(Sigma_true, ord="fro") ** 2), 1e-16)
        nmse_list.append(10.0 * np.log10(num / den + 1e-16))
        # BER proxy: hard-decision on sign of Re(g_hat) vs Re(g_true) mapped to QPSK-ish bit
        bits_true = (np.real(g_true) > 0).astype(int)
        bits_hat = (np.real(out["gains"]) > 0).astype(int)
        ber_list.append(float(np.mean(bits_hat != bits_true)))

    return {
        "rho": float(rho),
        "snr_db": float(snr_db),
        "profile": list(profile),
        "detector": detector,
        "n_trials": n_trials,
        "ber": float(np.mean(ber_list)),
        "cluster_pd": float(np.mean(pd_list)),
        "nmse_db": float(np.mean(nmse_list)),
    }


def run_p2_sweep(rho_values, snr_db_values, cluster_profiles, detectors,
                 n_trials, N, output_dir, seed=0):
    entries = []
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    idx = 0
    for rho in rho_values:
        for snr_db in snr_db_values:
            for profile in cluster_profiles:
                for detector in detectors:
                    entry = _eval_cell(rho, snr_db, profile, detector, n_trials, N, seed + idx)
                    entries.append(entry)
                    idx += 1
    ts = time.strftime("%Y%m%d-%H%M%S")
    out_json = output_dir / f"summary_{ts}.json"
    report = {"timestamp": ts, "entries": entries}
    out_json.write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[3]
    out = ROOT / "Sim" / "experiments" / "Results" / "p2"
    run_p2_sweep(
        rho_values=[0.0, 0.3, 0.5, 0.7, 0.9],
        snr_db_values=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0],
        cluster_profiles=[[2, 2], [3, 3, 2], [5]],
        detectors=["naive", "lrt"],
        n_trials=30,
        N=32,
        output_dir=out,
        seed=1,
    )
