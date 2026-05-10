"""Position-dependent Rayleigh quotient P_res_eff(k) validation."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_closed_form import closed_form_sigma
from p_res_eff import p_res_eff_from_sigma, p_res_eff_mc


def _unit_cols(N, K, rng):
    M = rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))
    M /= np.linalg.norm(M, axis=0, keepdims=True)
    return M


def test_p_res_eff_scalar_sigma_is_constant():
    rng = np.random.default_rng(10)
    N, K = 12, 5
    Sigma = 0.7 * np.eye(N, dtype=complex)
    H = _unit_cols(N, K, rng)
    eff = p_res_eff_from_sigma(H, Sigma)
    assert np.allclose(eff, 0.7)


def test_p_res_eff_mc_matches_closed_form():
    rng = np.random.default_rng(11)
    N, K = 16, 4
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Sigma = A @ A.conj().T / N + 0.05 * np.eye(N)
    H = _unit_cols(N, K, rng)
    eff_cf = p_res_eff_from_sigma(H, Sigma)
    eff_mc = p_res_eff_mc(H, Sigma, n_trials=6000, rng=rng)
    rel_err = np.max(np.abs(eff_mc - eff_cf) / np.abs(eff_cf))
    assert rel_err < 0.08, f"max rel err {rel_err:.3f} > 0.08"
