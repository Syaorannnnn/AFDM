"""BDLR-aware per-candidate CFAR threshold for LRT detector.

eta_new(k) = sigma_eff2 * Phi_k^H Sigma^{-1} Phi_k * ln(M_cand / P_fa)

Degenerates to P1 scalar cfar formula when Sigma = alpha * I.
See spec §2.2 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

import numpy as np


def cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0):
    """Per-candidate CFAR threshold for a single Phi_k."""
    if P_fa <= 0 or P_fa >= 1:
        raise ValueError(f"P_fa must be in (0,1), got {P_fa}")
    if M_cand < 1:
        raise ValueError(f"M_cand must be >=1, got {M_cand}")
    denom_term = Phi_k.conj().T @ Sigma_inv @ Phi_k
    rayleigh = float(np.trace(denom_term).real)
    if rayleigh <= 0:
        return 0.0
    return sigma_eff2 * rayleigh * float(np.log(M_cand / P_fa))


def cfar_threshold_lrt_batch(Phi_stack, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0):
    """Batch over M candidate operators."""
    M = Phi_stack.shape[0]
    out = np.empty(M, dtype=float)
    for m in range(M):
        out[m] = cfar_threshold_lrt(Phi_stack[m], Sigma_inv, M_cand, P_fa, sigma_eff2)
    return out
