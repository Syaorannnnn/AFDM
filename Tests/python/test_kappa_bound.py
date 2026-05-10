"""Kappa upper bound from spec §3.3."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_closed_form import closed_form_sigma, bdlr_decompose
from kappa_bound import kappa_weyl_upper_bound


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def test_kappa_bound_dominates_empirical_kappa():
    rng = np.random.default_rng(40)
    N = 24
    cluster_sizes = [3, 2, 3]
    C = len(cluster_sizes)
    P = sum(cluster_sizes)
    Phi = _make_random_unitary_phis(P, N, rng)
    labels = np.repeat(np.arange(C), cluster_sizes)
    P_per = np.array([0.4, 0.3, 0.3])
    for rho_val in [0.0, 0.3, 0.7, 0.95]:
        rho_per = rho_val * np.ones(C)
        coherent, scalar, N0_term = bdlr_decompose(
            Phi, labels, P_per, rho_per, sigma_d2=1.0, N0=0.01,
        )
        Sigma = coherent + scalar + N0_term
        emp_kappa = np.linalg.cond(Sigma)
        bound = kappa_weyl_upper_bound(
            coherent_part=coherent, scalar_weight=scalar[0, 0].real, N0=0.01,
        )
        assert bound >= emp_kappa * (1.0 - 1e-6), (
            f"bound {bound:.2f} < empirical kappa {emp_kappa:.2f} at rho={rho_val}"
        )
