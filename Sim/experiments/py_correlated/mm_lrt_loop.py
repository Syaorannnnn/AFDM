# Sim/experiments/py_correlated/mm_lrt_loop.py
"""MM double-loop driver for LRT-CLIP.

Implements spec §2.3 + §2.4 of 2026-05-11-p2-lrt-clip-mm-framework-design.md:
outer loop updates cluster gains g via proximal gradient with group soft-threshold,
inner loop updates Sigma via MoM closed form. Stopping on relative NLL change.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np


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
    raise NotImplementedError


def mm_double_loop(r, Phi_per_cluster, N0, cfg):
    raise NotImplementedError
