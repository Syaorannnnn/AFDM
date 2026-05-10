"""Rayleigh-quotient effective residual power at each candidate position."""

import numpy as np


def p_res_eff_from_sigma(H, Sigma):
    """
    Parameters
    ----------
    H : (N, K) complex ndarray of candidate responses (columns = candidates)
    Sigma : (N, N) Hermitian PSD
    Returns
    -------
    (K,) float array: h_k^H Sigma h_k / ||h_k||^2
    """
    num = np.einsum("nk,nm,mk->k", H.conj(), Sigma, H).real
    den = np.einsum("nk,nk->k", H.conj(), H).real
    return num / den


def p_res_eff_mc(H, Sigma, n_trials, rng=None):
    """Monte Carlo: draw e ~ CN(0, Sigma) and average |h_k^H e|^2 / ||h_k||^2."""
    if rng is None:
        rng = np.random.default_rng()
    N, K = H.shape
    # Use eigendecomposition to sample CN(0, Sigma).
    w, V = np.linalg.eigh(Sigma)
    w_clip = np.clip(w, 0, None)
    L = V * np.sqrt(w_clip)
    accum = np.zeros(K)
    norms = np.einsum("nk,nk->k", H.conj(), H).real
    for _ in range(n_trials):
        z = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
        e = L @ z
        proj = H.conj().T @ e
        accum += (np.abs(proj) ** 2) / norms
    return accum / n_trials
