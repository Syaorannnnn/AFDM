"""Sanity test that the new module imports and exposes the expected API."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_module_importable():
    import lrt_detector
    assert hasattr(lrt_detector, "lrt_statistic")
    assert hasattr(lrt_detector, "lrt_statistic_batch")
    assert hasattr(lrt_detector, "lrt_snr_gain")


from lrt_detector import lrt_statistic, lrt_statistic_batch, lrt_snr_gain


def _make_random_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_scalar_sigma_snr_gain_is_zero_db():
    rng = np.random.default_rng(100)
    N = 24
    Phi_k = _make_random_unitary(N, rng)
    alpha = 0.3
    Sigma = alpha * np.eye(N, dtype=complex)
    gain_db = lrt_snr_gain(Phi_k, Sigma)
    assert abs(gain_db) < 1e-10, f"expected 0 dB, got {gain_db:.3e} dB"


def test_lrt_statistic_batch_matches_single():
    rng = np.random.default_rng(101)
    N, M = 16, 5
    r = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    Phi_stack = np.stack([_make_random_unitary(N, rng) for _ in range(M)], axis=0)
    Sigma = 0.4 * np.eye(N, dtype=complex) + 0.1 * np.ones((N, N), dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    batch = lrt_statistic_batch(r, Phi_stack, Sigma_inv)
    singles = np.array([
        lrt_statistic(r, Phi_stack[m], Sigma_inv) for m in range(M)
    ])
    assert batch.shape == (M,)
    assert np.allclose(batch, singles, rtol=1e-12, atol=1e-14)


def test_lrt_snr_gain_positive_for_correlated_sigma():
    rng = np.random.default_rng(102)
    N = 32
    Phi_k = _make_random_unitary(N, rng)
    # Rank-1 coherent + identity floor
    v = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    Sigma = 0.05 * np.eye(N, dtype=complex) + np.outer(v, v.conj())
    gain_db = lrt_snr_gain(Phi_k, Sigma)
    assert gain_db > 0.5, f"expected >0.5 dB on correlated Sigma, got {gain_db:.3f}"
