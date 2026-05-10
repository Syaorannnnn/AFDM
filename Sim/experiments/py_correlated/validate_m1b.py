"""End-to-end M1-B numerical validation sweep.

Usage
-----
    python Sim/experiments/py_correlated/validate_m1b.py

Produces
--------
    Sim/experiments/Results/m1b_validation_<timestamp>.json
    Sim/experiments/Results/m1b_kappa_vs_rho_<timestamp>.png
"""

import json
import time
from pathlib import Path

import numpy as np

from sigma_closed_form import closed_form_sigma, mc_sigma, bdlr_decompose
from kappa_bound import kappa_weyl_upper_bound


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def _cluster_R(cluster_sizes, rho_val, gain_var=1.0):
    P = int(np.sum(cluster_sizes))
    R = np.zeros((P, P), dtype=complex)
    offset = 0
    for sz in cluster_sizes:
        block = gain_var * (
            (1 - rho_val) * np.eye(sz) + rho_val * np.ones((sz, sz))
        )
        R[offset:offset + sz, offset:offset + sz] = block
        offset += sz
    return R


def run_sweep(rho_values, snr_db_values, cluster_sizes, N, n_trials,
              output_dir, seed=0):
    rng = np.random.default_rng(seed)
    P = int(np.sum(cluster_sizes))
    Phi = _make_random_unitary_phis(P, N, rng)
    labels = np.repeat(np.arange(len(cluster_sizes)), cluster_sizes)
    P_per = np.array([sz / P for sz in cluster_sizes], dtype=float)
    entries = []
    for rho_val in rho_values:
        for snr_db in snr_db_values:
            sigma_d2 = 1.0
            N0 = sigma_d2 * 10 ** (-snr_db / 10.0)
            R = _cluster_R(cluster_sizes, rho_val)
            Sigma_cf = closed_form_sigma(Phi, R, sigma_d2, N0)
            Sigma_mc = mc_sigma(Phi, R, sigma_d2, N0, n_trials, rng=rng)
            rel = (
                np.linalg.norm(Sigma_mc - Sigma_cf, "fro")
                / np.linalg.norm(Sigma_cf, "fro")
            )
            coherent, scalar, N0_term = bdlr_decompose(
                Phi, labels, P_per, rho_val * np.ones(len(cluster_sizes)),
                sigma_d2, N0,
            )
            kappa_bound = kappa_weyl_upper_bound(
                coherent, scalar_weight=scalar[0, 0].real, N0=N0,
            )
            kappa_emp = float(np.linalg.cond(Sigma_cf))
            entries.append({
                "rho": float(rho_val),
                "snr_db": float(snr_db),
                "sigma_rel_err": float(rel),
                "kappa_empirical": kappa_emp,
                "kappa_bound": float(kappa_bound),
            })
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "timestamp": timestamp,
        "config": {
            "rho_values": list(rho_values),
            "snr_db_values": list(snr_db_values),
            "cluster_sizes": list(cluster_sizes),
            "N": N, "n_trials": n_trials, "seed": seed,
        },
        "entries": entries,
    }
    out_json = output_dir / f"m1b_validation_{timestamp}.json"
    out_json.write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[3]
    out = ROOT / "Sim" / "experiments" / "Results"
    run_sweep(
        rho_values=[0.0, 0.3, 0.5, 0.7, 0.9],
        snr_db_values=[0.0, 10.0, 20.0],
        cluster_sizes=[3, 3, 2],
        N=32, n_trials=2000, output_dir=out, seed=1,
    )
