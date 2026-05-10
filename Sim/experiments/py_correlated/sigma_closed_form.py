"""Closed-form ID2P covariance Sigma = sum R_ij Phi_i Phi_j^H * sigma_d2 + N0 I.

Implements M1-B spec §3 / 2026-05-07 §3. Monte Carlo validator averages the
same quantity over BDLR-sampled gains to confirm the closed form within
tolerance.
"""

import numpy as np


def closed_form_sigma(Phi_list, R, sigma_d2, N0):
    raise NotImplementedError


def mc_sigma(Phi_list, R, sigma_d2, N0, n_trials, rng=None):
    raise NotImplementedError
