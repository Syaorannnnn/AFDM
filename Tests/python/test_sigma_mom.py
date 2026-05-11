"""Method-of-moments cluster-gain estimator tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_mom import estimate_cluster_gains_mom, sigma_from_gains


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_sigma_from_gains_zero_gains_gives_noise_floor():
    rng = np.random.default_rng(300)
    N, C = 12, 3
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    gains = np.zeros(C, dtype=complex)
    N0 = 0.05
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    assert np.allclose(Sigma, N0 * np.eye(N), atol=1e-12)


def test_sigma_from_gains_reconstruction_matches_formula():
    rng = np.random.default_rng(301)
    N, C = 16, 2
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    gains = np.array([0.4 + 0.1j, 0.2 - 0.3j])
    N0 = 0.02
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    expected = N0 * np.eye(N, dtype=complex)
    for c in range(C):
        expected += (abs(gains[c]) ** 2) * (Phi_per_cluster[c] @ Phi_per_cluster[c].conj().T)
    assert np.allclose(Sigma, expected, atol=1e-12)


def test_mom_estimator_recovers_gains_within_MC_bias():
    """Oracle-support MoM: average |g_hat_c|^2 over MC should approach |g_true_c|^2."""
    rng = np.random.default_rng(302)
    N, C = 24, 2
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.7 + 0.0j, 0.5 + 0.0j])
    N0 = 0.01
    n_trials = 400
    accum_abs2 = np.zeros(C)
    for _ in range(n_trials):
        # draw coherent cluster gains and data
        h_per = np.array([g_true[c] * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2)) for c in range(C)])
        z = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
        r = sum(h_per[c] * (Phi_per_cluster[c] @ np.ones(N, dtype=complex) / np.sqrt(N)) for c in range(C)) + z
        support = np.ones(C, dtype=bool)
        g_hat = estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support)
        accum_abs2 += np.abs(g_hat) ** 2
    emp = accum_abs2 / n_trials
    true_abs2 = np.abs(g_true) ** 2
    rel = np.max(np.abs(emp - true_abs2) / (true_abs2 + 1e-12))
    assert rel < 0.35, f"MoM MC bias too large: rel={rel:.3f}"
