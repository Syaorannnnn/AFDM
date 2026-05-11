# Sim/experiments/py_correlated/mm_lrt_loop.py
"""MM double-loop driver for LRT-CLIP.

Implements spec §2.3 + §2.4 of 2026-05-11-p2-lrt-clip-mm-framework-design.md:
outer loop updates cluster gains g via proximal gradient with group soft-threshold,
inner loop updates Sigma via MoM closed form. Stopping on relative NLL change.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np

from sigma_mom import estimate_cluster_gains_mom, sigma_from_gains
from kappa_bound import kappa_weyl_upper_bound  # reserved for diagnostics


@dataclass
class MMConfig:
    lambda_gl: float = 0.1
    T_in_max: int = 3
    T_out_max: int = 15
    nll_tol: float = 1e-4
    kappa_star: float = 100.0
    lipschitz_cap: float = 1e4


@dataclass
class MMResult:
    gains: np.ndarray
    support: np.ndarray
    sigma: np.ndarray
    sigma_inv: np.ndarray
    nll_trajectory: List[float] = field(default_factory=list)
    n_inner: int = 0
    n_outer: int = 0
    converged: bool = False
    fallback_count: int = 0


def tikhonov_inverse(Sigma, kappa_star, N0):
    """Eigendecomposition-based inverse with kappa-floor Tikhonov regularisation.

    Let lambda = eigenvalues of (Sigma + Sigma^H)/2 (clamped to >=0).
    eps = max(N0, lambda_max / kappa_star).
    Return V diag(1/(lambda_i + eps)) V^H.
    """
    H = (Sigma + Sigma.conj().T) / 2
    w, V = np.linalg.eigh(H)
    w_clip = np.clip(w, 0.0, None)
    lam_max = float(max(w_clip.max(), 0.0))
    eps = max(N0, lam_max / max(kappa_star, 1.0))
    inv_diag = 1.0 / (w_clip + eps)
    return V @ np.diag(inv_diag) @ V.conj().T


def _nll(r, Sigma, Sigma_inv):
    data = float(np.vdot(r, Sigma_inv @ r).real)
    # log det via eigvals of Sigma (PSD)
    w = np.linalg.eigvalsh((Sigma + Sigma.conj().T) / 2)
    w = np.clip(w, 1e-16, None)
    logdet = float(np.sum(np.log(w)))
    return data + logdet


def _group_soft_threshold(g, mu):
    """Apply soft-threshold per scalar g_c (cluster-level group lasso with group size 1)."""
    abs_g = abs(g)
    if abs_g <= mu:
        return 0.0 + 0.0j
    return (1.0 - mu / abs_g) * g


def _lipschitz_bound(Phi_per_cluster, Sigma_inv, cap):
    """Rough Lipschitz bound of grad_g NLL: L = 2 * max_c ||Phi_c^H Sigma^{-1} Phi_c||_F, capped."""
    L = 0.0
    for Phi_c in Phi_per_cluster:
        M = Phi_c.conj().T @ Sigma_inv @ Phi_c
        L = max(L, float(np.linalg.norm(M, ord="fro")))
    return min(2.0 * L, cap)


def mm_double_loop(r, Phi_per_cluster, N0, cfg):
    """Driver for spec §2.3 + §2.4."""
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    gains = np.zeros(C, dtype=complex)
    support = np.ones(C, dtype=bool)  # start with full support; prox will sparsify
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    Sigma_inv = tikhonov_inverse(Sigma, cfg.kappa_star, N0)
    traj = [_nll(r, Sigma, Sigma_inv)]
    fallback_count = 0
    n_inner_total = 0
    n_outer = 0
    converged = False

    for t in range(cfg.T_out_max):
        # --- Inner loop: MoM on current support ---
        g_prev = gains.copy()
        for _ in range(cfg.T_in_max):
            n_inner_total += 1
            gains = estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support)
            Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
            Sigma_inv = tikhonov_inverse(Sigma, cfg.kappa_star, N0)
            nll = _nll(r, Sigma, Sigma_inv)
            if traj and nll > traj[-1] + 1e-12:
                # fallback: revert this inner update
                gains = g_prev.copy()
                Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
                Sigma_inv = tikhonov_inverse(Sigma, cfg.kappa_star, N0)
                fallback_count += 1
                break
            g_prev = gains.copy()

        # --- Outer loop: proximal gradient step ---
        L = _lipschitz_bound(Phi_per_cluster, Sigma_inv, cfg.lipschitz_cap)
        if L <= 0:
            L = 1.0
        step = 1.0 / L
        # gradient of data term wrt g_c (treating g_c as complex):
        # d(data)/d g_c^* = - |g_c|^2 ... simplified proxy: use residual projection
        grads = np.zeros(C, dtype=complex)
        x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
        residual = r.copy()
        for c in range(C):
            v_c = Phi_per_cluster[c] @ x_unit
            grads[c] = -2.0 * np.vdot(v_c, residual - gains[c] * v_c) / max(float(np.vdot(v_c, v_c).real), 1e-16)
        gains_tentative = gains - step * grads
        new_gains = np.array([_group_soft_threshold(gains_tentative[c], cfg.lambda_gl * step) for c in range(C)], dtype=complex)
        new_support = np.abs(new_gains) > 1e-12
        Sigma_new = sigma_from_gains(new_gains, Phi_per_cluster, N0)
        Sigma_inv_new = tikhonov_inverse(Sigma_new, cfg.kappa_star, N0)
        nll_new = _nll(r, Sigma_new, Sigma_inv_new)

        if nll_new > traj[-1] + 1e-12:
            # outer step would increase NLL: fall back (skip update, still logged)
            fallback_count += 1
            traj.append(traj[-1])
            n_outer += 1
            # check stopping on flat trajectory
            if len(traj) >= 2:
                rel = abs(traj[-1] - traj[-2]) / max(abs(traj[-2]), 1e-16)
                if rel < cfg.nll_tol:
                    converged = True
                    break
            continue

        gains = new_gains
        support = new_support
        Sigma = Sigma_new
        Sigma_inv = Sigma_inv_new
        traj.append(nll_new)
        n_outer += 1

        # stopping rule
        rel = abs(traj[-1] - traj[-2]) / max(abs(traj[-2]), 1e-16)
        if rel < cfg.nll_tol:
            converged = True
            break

    return MMResult(
        gains=gains,
        support=support,
        sigma=Sigma,
        sigma_inv=Sigma_inv,
        nll_trajectory=traj,
        n_inner=n_inner_total,
        n_outer=n_outer,
        converged=converged,
        fallback_count=fallback_count,
    )
