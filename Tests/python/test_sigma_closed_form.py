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


def test_bdlr_decompose_conditional_rank_equals_cluster_count():
    """Given a pilot vector x_p, the conditional-on-data coherent kernel
    sum_c P_c (G_c x_p) (G_c x_p)^H has rank exactly C at rho=1.

    Note: the ensemble-averaged operator `bdlr_decompose(...)[0]` is
    generically rank-N; the rank-C collapse is a *conditional* property
    that the Weyl bound and identifiability arguments rely on. This test
    builds the conditional kernel inline and checks rank.
    """
    rng = np.random.default_rng(2)
    N = 24
    cluster_sizes = [3, 2, 4]
    C = len(cluster_sizes)
    P = sum(cluster_sizes)
    Phi = _make_random_unitary_phis(P, N, rng)
    cluster_labels = np.repeat(np.arange(C), cluster_sizes)
    P_per = np.array([0.5, 0.3, 0.2])
    x_p = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
    coherent_cond = np.zeros((N, N), dtype=complex)
    for c in range(C):
        members = np.where(cluster_labels == c)[0]
        Gc = sum(Phi[i] for i in members)
        gc_x = Gc @ x_p
        coherent_cond += P_per[c] * np.outer(gc_x, gc_x.conj())
    eigvals = np.sort(np.linalg.eigvalsh(coherent_cond))[::-1]
    assert np.sum(eigvals > 1e-8) == C, (
        f"conditional rank={np.sum(eigvals > 1e-8)}, expected {C}"
    )


def test_bdlr_decompose_rho_zero_coherent_vanishes():
    """When rho_c = 0 for every cluster, coherent part is zero."""
    rng = np.random.default_rng(3)
    N, C, per = 16, 2, 3
    Phi = _make_random_unitary_phis(C * per, N, rng)
    from sigma_closed_form import bdlr_decompose
    coherent, scalar, N0_term = bdlr_decompose(
        Phi_list=Phi,
        cluster_labels=np.repeat(np.arange(C), per),
        P_per_cluster=np.array([0.5, 0.5]),
        rho_per_cluster=np.zeros(C),
        sigma_d2=1.0,
        N0=0.01,
    )
    assert np.allclose(coherent, 0.0, atol=1e-12)
    # scalar part is diagonal.
    assert np.allclose(scalar, np.trace(scalar).real / scalar.shape[0] * np.eye(scalar.shape[0]))


def test_bdlr_decompose_operator_rho_zero_gives_scalar_sigma():
    """At rho=0 the whole Sigma reconstructed from (coherent, scalar, N0_term)
    equals a scalar multiple of I."""
    rng = np.random.default_rng(4)
    N, C, per = 12, 3, 2
    Phi = _make_random_unitary_phis(C * per, N, rng)
    from sigma_closed_form import bdlr_decompose
    P_per = np.array([0.4, 0.35, 0.25])
    coherent, scalar, N0_term = bdlr_decompose(
        Phi_list=Phi,
        cluster_labels=np.repeat(np.arange(C), per),
        P_per_cluster=P_per,
        rho_per_cluster=np.zeros(C),
        sigma_d2=1.0,
        N0=0.05,
    )
    Sigma = coherent + scalar + N0_term
    expected = (np.sum(P_per) + 0.05) * np.eye(N)
    assert np.allclose(Sigma, expected, atol=1e-12)
