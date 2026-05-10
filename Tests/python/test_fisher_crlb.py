"""Fisher / CRLB sanity tests for M1-B claims."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from fisher_crlb import fisher_single_path, crlb_from_fisher
from fisher_crlb import fisher_two_paths_joint


def _make_delay_operator(N, delay_frac):
    """Simple shift-by-delay Phi: (N, N) circulant shift matrix sampled at frac delay."""
    k = np.arange(N)
    col = np.exp(-1j * 2 * np.pi * delay_frac * k / N)
    # Symbolic "derivative w.r.t. delay" is encoded in fisher_single_path test itself.
    return np.diag(col)


def test_fisher_matrix_hermitian_and_psd():
    rng = np.random.default_rng(20)
    N = 32
    Phi = _make_delay_operator(N, 0.3)
    dPhi = _make_delay_operator(N, 0.3001)
    dPhi = (dPhi - Phi) / 0.0001
    Sigma = 0.1 * np.eye(N) + 0.05 * np.ones((N, N))
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    h = 1.0 + 0j
    J = fisher_single_path(
        Phi=Phi, dPhi_dtheta=[dPhi], h=h, x_pilot=x_pilot, Sigma=Sigma
    )
    assert J.shape == (3, 3)  # 1 delay-like param + Re(h) + Im(h)
    assert np.allclose(J, J.conj().T, atol=1e-10)
    assert np.all(np.linalg.eigvalsh((J + J.conj().T) / 2) >= -1e-8)


def test_two_path_crlb_diverges_as_phis_become_coherent():
    """As Phi_j -> Phi_i, joint Fisher block near-singular -> CRLB -> inf."""
    rng = np.random.default_rng(30)
    N = 48
    base = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    U, _, _ = np.linalg.svd(base)
    Phi_i = U                              # unitary
    Sigma = 0.05 * np.eye(N, dtype=complex)
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    # Sweep similarity: Phi_j = cos(theta) Phi_i + sin(theta) Phi_perp
    V, _, _ = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )
    Phi_perp = V
    crlb_gain_i = []
    thetas = [1.0, 0.3, 0.1, 0.03, 0.01]
    for theta in thetas:
        Phi_j = np.cos(theta) * Phi_i + np.sin(theta) * Phi_perp
        # Renormalise so Phi_j is still approximately unitary for the test.
        J = fisher_two_paths_joint(
            Phi_i=Phi_i, Phi_j=Phi_j, h_i=1.0, h_j=1.0,
            x_pilot=x_pilot, Sigma=Sigma,
        )
        diag = crlb_from_fisher(J, eps_reg=1e-14)
        crlb_gain_i.append(diag[0])  # Re(h_i) variance
    # Expect monotone growth as theta shrinks.
    for k in range(len(thetas) - 1):
        assert crlb_gain_i[k + 1] > crlb_gain_i[k], (
            f"CRLB did not grow at theta step {k}: {crlb_gain_i}"
        )
    # Final (near-coherent) CRLB should be at least 10x the separated case.
    assert crlb_gain_i[-1] > 10.0 * crlb_gain_i[0]


from fisher_crlb import fisher_cluster_level


def test_cluster_level_crlb_bounded_across_similarity_sweep():
    rng = np.random.default_rng(31)
    N = 48
    U = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )[0]
    V = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )[0]
    Phi_i, Phi_perp = U, V
    Sigma = 0.05 * np.eye(N, dtype=complex)
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    cluster_crlbs = []
    for theta in [1.0, 0.3, 0.1, 0.03, 0.01]:
        Phi_j = np.cos(theta) * Phi_i + np.sin(theta) * Phi_perp
        Phi_cluster = Phi_i + Phi_j  # unnormalised coherent sum, matches spec §3.1
        J_cluster = fisher_cluster_level(
            Phi_cluster=Phi_cluster, g_cluster=1.0, x_pilot=x_pilot, Sigma=Sigma
        )
        diag = crlb_from_fisher(J_cluster, eps_reg=1e-14)
        cluster_crlbs.append(diag[0])
    # Cluster-level CRLB should stay within one order of magnitude across the sweep.
    ratio = max(cluster_crlbs) / min(cluster_crlbs)
    assert ratio < 10.0, f"Cluster CRLB varied by {ratio:.2f}x across similarity sweep"
