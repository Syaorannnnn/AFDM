"""Fisher / CRLB sanity tests for M1-B claims."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from fisher_crlb import fisher_single_path, crlb_from_fisher


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
