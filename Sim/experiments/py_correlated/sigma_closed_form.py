"""Closed-form ID2P covariance Sigma = sum R_ij Phi_i Phi_j^H * sigma_d2 + N0 I.

Implements M1-B spec §3 / 2026-05-07 §3. Monte Carlo validator averages the
same quantity over BDLR-sampled gains to confirm the closed form within
tolerance.
"""

import numpy as np


def closed_form_sigma(Phi_list, R, sigma_d2, N0):
    """
    Closed-form ID2P covariance Sigma under BDLR channel prior.

    Sigma = sigma_d2 * sum_{i,j} R[i,j] * Phi_i @ Phi_j.conj().T + N0 * I_N

    Parameters
    ----------
    Phi_list : list of (N, N) complex ndarray
        Per-path DAFT-domain operators; must satisfy Phi_i Phi_i^H = I_N.
    R : (P, P) complex ndarray
        Path-gain covariance (Hermitian PSD).
    sigma_d2 : float
        Data symbol power.
    N0 : float
        Noise power.

    Returns
    -------
    Sigma : (N, N) complex ndarray
    """
    P = len(Phi_list)
    assert R.shape == (P, P), f"R must be {P}x{P}, got {R.shape}"
    N = Phi_list[0].shape[0]
    Sigma = np.zeros((N, N), dtype=complex)
    for i in range(P):
        for j in range(P):
            if R[i, j] == 0:
                continue
            Sigma += R[i, j] * (Phi_list[i] @ Phi_list[j].conj().T)
    Sigma = sigma_d2 * Sigma + N0 * np.eye(N)
    return Sigma


def mc_sigma(Phi_list, R, sigma_d2, N0, n_trials, rng=None):
    """
    Monte Carlo estimate of the ID2P covariance by averaging h h^H and
    forming the outer channel response.
    """
    if rng is None:
        rng = np.random.default_rng()
    P = len(Phi_list)
    N = Phi_list[0].shape[0]
    L = np.linalg.cholesky(R + 1e-12 * np.eye(P))
    Sigma_accum = np.zeros((N, N), dtype=complex)
    for _ in range(n_trials):
        w = (rng.standard_normal(P) + 1j * rng.standard_normal(P)) / np.sqrt(2)
        h = L @ w
        x = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
        x = np.sqrt(sigma_d2) * x
        H_eff = sum(h[i] * Phi_list[i] for i in range(P))
        y = H_eff @ x
        Sigma_accum += np.outer(y, y.conj())
    Sigma_mc = Sigma_accum / n_trials + N0 * np.eye(N)
    return Sigma_mc
