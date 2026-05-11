"""LRT-optimal detector: T(r;k) = |Phi_k^H Sigma^{-1} r|^2 / (Phi_k^H Sigma^{-1} Phi_k).

Implements spec §2.1 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

import numpy as np


def lrt_statistic(r, Phi_k, Sigma_inv):
    """Single-candidate LRT statistic.

    Parameters
    ----------
    r : (N,) complex ndarray
    Phi_k : (N, N) complex ndarray — candidate response operator
    Sigma_inv : (N, N) complex ndarray — precomputed inverse covariance
    """
    w = Phi_k.conj().T @ Sigma_inv @ r
    num = float(np.vdot(w, w).real)
    norm_term = Phi_k.conj().T @ Sigma_inv @ Phi_k
    den = float(np.trace(norm_term).real)
    if den <= 0:
        return 0.0
    return num / den


def lrt_statistic_batch(r, Phi_stack, Sigma_inv):
    """Batch over M candidate operators."""
    M = Phi_stack.shape[0]
    out = np.empty(M, dtype=float)
    for m in range(M):
        out[m] = lrt_statistic(r, Phi_stack[m], Sigma_inv)
    return out


def lrt_snr_gain(Phi_k, Sigma):
    """SNR gain in dB of LRT over naive matched filter at candidate k."""
    Sigma_inv = np.linalg.pinv(Sigma)
    num = np.trace(Phi_k.conj().T @ Sigma_inv @ Phi_k).real
    mid = np.trace(Phi_k.conj().T @ Sigma @ Phi_k).real
    den = np.linalg.norm(Phi_k, "fro") ** 4
    if den <= 0 or num <= 0 or mid <= 0:
        return 0.0
    ratio = (num * mid) / den
    return 10.0 * float(np.log10(max(ratio, 1e-300)))
