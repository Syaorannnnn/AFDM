"""Sanity test that the new module imports and exposes the expected API."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_module_importable():
    import sigma_closed_form
    assert hasattr(sigma_closed_form, "closed_form_sigma")
    assert hasattr(sigma_closed_form, "mc_sigma")


from sigma_closed_form import closed_form_sigma, mc_sigma
import numpy as np


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def test_diagonal_R_yields_scalar_sigma():
    rng = np.random.default_rng(0)
    N, P = 16, 4
    Phi = _make_random_unitary_phis(P, N, rng)
    R = np.diag([0.25, 0.25, 0.25, 0.25]).astype(complex)
    sigma_d2, N0 = 1.0, 0.01
    Sigma = closed_form_sigma(Phi, R, sigma_d2, N0)
    expected_scalar = sigma_d2 * np.trace(R).real + N0
    assert np.allclose(Sigma, expected_scalar * np.eye(N), atol=1e-10)


def test_mc_sigma_converges_to_closed_form():
    rng = np.random.default_rng(1)
    N, P = 8, 3
    Phi = _make_random_unitary_phis(P, N, rng)
    # Full-rank correlated R (simple BDLR with one cluster).
    rho, gain_var = 0.7, 0.25
    R = gain_var * ((1 - rho) * np.eye(P) + rho * np.ones((P, P)))
    sigma_d2, N0 = 1.0, 0.02
    Sigma_cf = closed_form_sigma(Phi, R, sigma_d2, N0)
    Sigma_mc = mc_sigma(Phi, R, sigma_d2, N0, n_trials=4000, rng=rng)
    rel_err = np.linalg.norm(Sigma_mc - Sigma_cf, "fro") / np.linalg.norm(Sigma_cf, "fro")
    assert rel_err < 0.08, f"MC/CF disagreement {rel_err:.3f} > 0.08"
