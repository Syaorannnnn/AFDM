"""Closed-form condition-number upper bound for Sigma under BDLR."""

import numpy as np


def kappa_weyl_upper_bound(coherent_part, scalar_weight, N0):
    """
    Weyl-based upper bound on kappa(Sigma) when
       Sigma = coherent_part + scalar_weight * I + N0 * I.

    lambda_max(Sigma) <= lambda_max(coherent_part) + scalar_weight + N0
    lambda_min(Sigma) >= scalar_weight + N0

    so kappa(Sigma) <= (lambda_max(coherent_part) + scalar_weight + N0)
                        / (scalar_weight + N0).
    """
    eigvals = np.linalg.eigvalsh((coherent_part + coherent_part.conj().T) / 2)
    lam_max = float(max(eigvals.max(), 0.0))
    denom = scalar_weight + N0
    if denom <= 0:
        raise ValueError("scalar_weight + N0 must be positive")
    return (lam_max + scalar_weight + N0) / denom
