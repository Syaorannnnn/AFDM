"""Fisher matrix and CRLB for single-path parameters theta = [shape..., Re h, Im h]."""

import numpy as np


def fisher_single_path(Phi, dPhi_dtheta, h, x_pilot, Sigma):
    """
    Fisher information for theta = [theta_1, ..., theta_L, Re(h), Im(h)].

    The forward model is mu(theta) = h * Phi(theta) @ x_pilot, with observation
    y ~ CN(mu, Sigma). For Gaussian observation the standard form is
      J_ij = 2 Re[ (d mu / d theta_i)^H Sigma^{-1} (d mu / d theta_j) ].

    Parameters
    ----------
    Phi : (N, N) complex ndarray at the reference point.
    dPhi_dtheta : list of (N, N) complex ndarrays — partial derivatives wrt the
        shape parameters (delay, fractional Doppler, ...). Length L.
    h : complex scalar path gain.
    x_pilot : (N,) complex pilot symbol vector.
    Sigma : (N, N) Hermitian PSD residual covariance.

    Returns
    -------
    J : (L+2, L+2) real ndarray (Hermitian symmetric, PSD).
    """
    N = Phi.shape[0]
    L = len(dPhi_dtheta)
    Sigma_inv = np.linalg.pinv(Sigma)
    # Build (L+2) Jacobian columns.
    cols = []
    for d in dPhi_dtheta:
        cols.append(h * (d @ x_pilot))
    cols.append(Phi @ x_pilot)          # d mu / d Re(h)
    cols.append(1j * (Phi @ x_pilot))   # d mu / d Im(h)
    G = np.stack(cols, axis=1)           # (N, L+2)
    J_complex = G.conj().T @ Sigma_inv @ G
    J = 2.0 * J_complex.real
    return J


def crlb_from_fisher(J, eps_reg=1e-12):
    """Return diag(inv(J)) with tiny Tikhonov for near-singular J."""
    Jr = J + eps_reg * np.eye(J.shape[0])
    return np.diag(np.linalg.inv(Jr))


def fisher_two_paths_joint(Phi_i, Phi_j, h_i, h_j, x_pilot, Sigma):
    """
    Joint Fisher matrix for theta = [Re(h_i), Im(h_i), Re(h_j), Im(h_j)]
    with Phi_i, Phi_j fixed.
    """
    Sigma_inv = np.linalg.pinv(Sigma)
    a = Phi_i @ x_pilot
    b = Phi_j @ x_pilot
    cols = [a, 1j * a, b, 1j * b]
    G = np.stack(cols, axis=1)
    J = 2.0 * (G.conj().T @ Sigma_inv @ G).real
    return J
