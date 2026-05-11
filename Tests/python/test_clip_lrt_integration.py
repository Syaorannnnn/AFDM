"""clip_lrt end-to-end smoke tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from clip_lrt import clip_lrt_receive, ClipLrtConfig


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def _make_synthetic_observation(N, C, rng):
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.6 + 0.0j, 0.4 + 0.0j])
    N0 = 0.02
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    y = np.zeros(N, dtype=complex)
    for c in range(C):
        h = g_true[c] * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2))
        y = y + h * (Phi_per_cluster[c] @ x_unit)
    y = y + (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
    return y, Phi_per_cluster, N0


def test_clip_lrt_returns_expected_structure():
    rng = np.random.default_rng(500)
    y, Phi_per_cluster, N0 = _make_synthetic_observation(N=20, C=2, rng=rng)
    cfg = ClipLrtConfig(Phi_per_cluster=Phi_per_cluster, N0=N0)
    out = clip_lrt_receive(y, pilot_frame=None, config=cfg)
    for key in ["gains", "support", "sigma", "nll_trajectory", "n_outer", "stage_results"]:
        assert key in out, f"missing key {key}"
    # structural check: stage_results has 3 entries
    assert len(out["stage_results"]) == 3


def test_clip_lrt_is_seed_reproducible():
    rng_a = np.random.default_rng(501)
    rng_b = np.random.default_rng(501)
    y_a, Phi_a, N0_a = _make_synthetic_observation(N=16, C=2, rng=rng_a)
    y_b, Phi_b, N0_b = _make_synthetic_observation(N=16, C=2, rng=rng_b)
    cfg_a = ClipLrtConfig(Phi_per_cluster=Phi_a, N0=N0_a)
    cfg_b = ClipLrtConfig(Phi_per_cluster=Phi_b, N0=N0_b)
    out_a = clip_lrt_receive(y_a, pilot_frame=None, config=cfg_a)
    out_b = clip_lrt_receive(y_b, pilot_frame=None, config=cfg_b)
    assert np.allclose(out_a["gains"], out_b["gains"], atol=1e-12)
    assert np.array_equal(out_a["support"], out_b["support"])
