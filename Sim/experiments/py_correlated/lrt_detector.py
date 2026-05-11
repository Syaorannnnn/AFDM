# Sim/experiments/py_correlated/lrt_detector.py
"""LRT-optimal detector: T(r;k) = |Phi_k^H Sigma^{-1} r|^2 / (Phi_k^H Sigma^{-1} Phi_k).

Implements spec §2.1 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

import numpy as np


def lrt_statistic(r, Phi_k, Sigma_inv):
    raise NotImplementedError


def lrt_statistic_batch(r, Phi_stack, Sigma_inv):
    raise NotImplementedError


def lrt_snr_gain(Phi_k, Sigma):
    raise NotImplementedError
