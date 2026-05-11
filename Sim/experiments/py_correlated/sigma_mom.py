"""Method-of-moments cluster-gain estimator and Sigma reconstructor.

Implements the inner block of the MM double loop in spec §2.3 of
2026-05-11-p2-lrt-clip-mm-framework-design.md.

Data model (coherent BDLR, rho=1 per cluster):
    r = sum_c g_c * (Phi_c @ x_pilot) + z,   z ~ CN(0, N0 I)

Given support mask, MoM estimator projects r onto each cluster response
and returns the closed-form magnitude-matched gain.
"""

import numpy as np


def estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support_mask):
    """Closed-form cluster-gain estimates via method of moments.

    For each active cluster c, solves the scalar LS
        g_hat_c = argmin_g || r - g * (sum x_pilot phased by Phi_c) ||^2
    using Phi_c @ 1 / sqrt(N) as the canonical response direction
    (pilot-symbol value drops into the 1/sqrt(N) unit vector).

    Parameters
    ----------
    r : (N,) complex ndarray
    Phi_per_cluster : length-C list of (N, N) complex ndarrays
    N0 : float
    support_mask : (C,) bool ndarray — True for active clusters

    Returns
    -------
    (C,) complex ndarray — zero for inactive clusters
    """
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    gains = np.zeros(C, dtype=complex)
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    for c in range(C):
        if not support_mask[c]:
            continue
        v_c = Phi_per_cluster[c] @ x_unit
        v_norm2 = float(np.vdot(v_c, v_c).real)
        if v_norm2 <= 0:
            continue
        gains[c] = np.vdot(v_c, r) / v_norm2
    return gains


def sigma_from_gains(gains, Phi_per_cluster, N0):
    """Reconstruct Sigma = N0 * I + sum_c |g_c|^2 * Phi_c Phi_c^H."""
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    Sigma = N0 * np.eye(N, dtype=complex)
    for c in range(C):
        w = abs(gains[c]) ** 2
        if w <= 0:
            continue
        Sigma = Sigma + w * (Phi_per_cluster[c] @ Phi_per_cluster[c].conj().T)
    return Sigma
