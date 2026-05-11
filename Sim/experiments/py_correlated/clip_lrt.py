"""End-to-end LRT-CLIP three-stage receiver.

Stage 0 (Coarse-Lock diagnostic): short MM outer pass producing a diagnostic
    snapshot of the gains and support that a full MM sweep would lock onto.
    This stage is currently a diagnostic checkpoint only — its output is
    recorded in `stage_results[0]` but is NOT yet wired as a warm-start into
    Stage 1 because `mm_double_loop` re-initialises internally. True warm-
    start seeding is deferred to the P3 integration task.
Stage 1 (Iteration): full MM double loop, independent of Stage 0.
Stage 2 (Polishing): hold the support learned in Stage 1, re-estimate gains
    via MoM, and pair a freshly computed `sigma_inv` with the polished `sigma`
    for downstream consumers.

See spec §3 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np

from mm_lrt_loop import MMConfig, mm_double_loop, tikhonov_inverse
from sigma_mom import estimate_cluster_gains_mom, sigma_from_gains


@dataclass
class ClipLrtConfig:
    Phi_per_cluster: List[np.ndarray] = field(default_factory=list)
    N0: float = 0.01
    lambda_gl: float = 0.05
    stage0_T_out: int = 2
    stage1_T_out: int = 15
    stage2_T_out: int = 5
    T_in_max: int = 3
    nll_tol: float = 1e-4
    kappa_star: float = 100.0


def _run_mm_stage(y, Phi_per_cluster, N0, T_out, lambda_gl, nll_tol, kappa_star):
    cfg = MMConfig(
        lambda_gl=lambda_gl,
        T_in_max=3,
        T_out_max=T_out,
        nll_tol=nll_tol,
        kappa_star=kappa_star,
    )
    return mm_double_loop(y, Phi_per_cluster, N0, cfg)


def clip_lrt_receive(y, pilot_frame, config):
    """Three-stage LRT-CLIP receiver.

    Parameters
    ----------
    y : (N,) complex ndarray
        Received DAFT-domain signal.
    pilot_frame : placeholder — not used in this minimal integration; reserved
        for future extension into the full GiFree pipeline.
    config : ClipLrtConfig

    Returns
    -------
    dict with keys: gains, support, sigma, sigma_inv, nll_trajectory, n_outer,
    stage_results (list of length 3 with per-stage MMResult summaries).
    """
    stage_results = []

    # Stage 0: Coarse-Lock — short MM run to seed support.
    r0 = _run_mm_stage(
        y, config.Phi_per_cluster, config.N0,
        T_out=config.stage0_T_out,
        lambda_gl=config.lambda_gl,
        nll_tol=config.nll_tol,
        kappa_star=config.kappa_star,
    )
    stage_results.append({"gains": r0.gains.copy(), "support": r0.support.copy(), "nll_final": r0.nll_trajectory[-1]})

    # Stage 1: Iteration — full MM loop. `mm_double_loop` currently re-initialises
    # gains/support internally, so Stage 0's output is not a warm-start; that
    # rewire is a P3 integration task.
    r1 = _run_mm_stage(
        y, config.Phi_per_cluster, config.N0,
        T_out=config.stage1_T_out,
        lambda_gl=config.lambda_gl,
        nll_tol=config.nll_tol,
        kappa_star=config.kappa_star,
    )
    stage_results.append({"gains": r1.gains.copy(), "support": r1.support.copy(), "nll_final": r1.nll_trajectory[-1]})

    # Stage 2: Polishing — hold support from Stage 1, re-estimate via MoM once with tight tolerance.
    support_polished = r1.support.copy()
    gains_polished = estimate_cluster_gains_mom(y, config.Phi_per_cluster, config.N0, support_polished)
    Sigma_polished = sigma_from_gains(gains_polished, config.Phi_per_cluster, config.N0)
    # Keep sigma_inv paired with the polished Sigma so downstream consumers
    # (e.g. P3 unrolled-EM teacher signals) see a self-consistent (Σ, Σ⁻¹) pair.
    Sigma_inv_polished = tikhonov_inverse(Sigma_polished, config.kappa_star, config.N0)
    stage_results.append({"gains": gains_polished.copy(), "support": support_polished.copy(), "nll_final": r1.nll_trajectory[-1]})

    return {
        "gains": gains_polished,
        "support": support_polished,
        "sigma": Sigma_polished,
        "sigma_inv": Sigma_inv_polished,
        "nll_trajectory": r1.nll_trajectory,
        "n_outer": r1.n_outer,
        "stage_results": stage_results,
    }
