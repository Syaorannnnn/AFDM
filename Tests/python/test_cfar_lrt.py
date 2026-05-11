"""BDLR-aware CFAR threshold tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from cfar_lrt import cfar_threshold_lrt, cfar_threshold_lrt_batch


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_scalar_sigma_matches_p1_cfar_formula():
    """When Sigma = alpha * I, eta_new(k) = alpha * ln(M/P_FA)."""
    rng = np.random.default_rng(200)
    N, M_cand, P_fa = 20, 7, 1e-3
    Phi_k = _rand_unitary(N, rng)
    alpha = 0.6
    Sigma = alpha * np.eye(N, dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    eta = cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0)
    # eta_new(k) = 1 * (||Phi_k||^2 / alpha) * ln(M/P_FA)
    expected = (np.linalg.norm(Phi_k, "fro") ** 2 / alpha) * np.log(M_cand / P_fa)
    assert np.isclose(eta, expected, rtol=1e-10)


def test_batch_matches_singles():
    rng = np.random.default_rng(201)
    N, M_cand, P_fa = 16, 5, 1e-3
    Phi_stack = np.stack([_rand_unitary(N, rng) for _ in range(M_cand)], axis=0)
    Sigma = 0.2 * np.eye(N, dtype=complex) + 0.05 * np.ones((N, N), dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    batch = cfar_threshold_lrt_batch(Phi_stack, Sigma_inv, M_cand, P_fa)
    singles = np.array([
        cfar_threshold_lrt(Phi_stack[m], Sigma_inv, M_cand, P_fa) for m in range(M_cand)
    ])
    assert batch.shape == (M_cand,)
    assert np.allclose(batch, singles, rtol=1e-12)
