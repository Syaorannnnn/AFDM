# P1 M1-B Numerical Validation — Diagnostic Notes

**Date:** 2026-05-10
**Plan:** `docs/superpowers/plans/2026-05-10-p1-m1b-numerical-validation.md`
**Spec:** `docs/superpowers/specs/2026-05-10-afdm-correlated-clip-research-design.md` §3

## Gate Values

The latest `validate_m1b.py` sweep (see `Sim/experiments/Results/m1b_validation_<timestamp>.json`)
is considered a PASS iff every entry satisfies:

- `sigma_rel_err < 0.15` at `n_trials=2000`
- `kappa_bound >= kappa_empirical` (Weyl bound is a true upper bound)
- `fisher_single_path` returns Hermitian PSD matrix for all tested configurations
- Two-path CRLB grows monotonically as `cos(angle(Phi_i, Phi_j)) -> 1`
- Cluster-level CRLB stays within one order of magnitude over the same sweep

## Observed Values (2026-05-10 sweep)

- sigma_rel_err across 15 entries: min 0.1164, max 0.1368, mean 0.1292 (all below 0.15).
- kappa_bound dominates kappa_empirical at every entry. Bound becomes loose as rho approaches 1: at rho=0.9, SNR=20 dB, bound = 45.43 vs empirical = 4.79 (ratio ~9.5x). This is expected Weyl-chain looseness.
- Two-path CRLB at seed=30, N=48, theta sweep [1.0, 0.3, 0.1, 0.03, 0.01]: CRLB grows strictly monotonically from 4.62e-4 to 3.02, final/first ratio ~6534x (vs 10x threshold).
- Cluster-level CRLB at seed=31: varies only 1.21x across the same theta sweep, confirming identifiability recovery.
- Conditional rank kernel at seed=2, cluster_sizes=[3,2,4]: exactly 3 nonzero eigenvalues above 1e-8, remaining 21 at machine epsilon.

## Known Limitations

- `mc_sigma` convergence is slow at high `rho` (rank-1 coherent part); smoke test
  uses n_trials=2000 with threshold 0.15 (O(1/sqrt(T)) scaling from the main MC
  test which uses 4000 trials at threshold 0.08).
- `fisher_single_path` uses `np.linalg.pinv(Sigma)`; when `Sigma` is near-singular
  at rho -> 1 the CRLB is sensitive to `eps_reg`. This is captured by the
  kappa-bound test and handled explicitly in P2.
- Weyl kappa bound becomes loose at high rho / high SNR (ratio ~10x at rho=0.9,
  20 dB). Tighter bounds (e.g. interlacing) are an open question for the letter.
- `validate_m1b.py` module docstring mentions a PNG output that the current
  implementation does not produce; plotting deferred.

## Downstream Consumers

- P2 `clip_lrt.py` reuses `closed_form_sigma` and the kappa bound as inputs.
- P3 metric implementations consume the MC/closed-form parity as ground-truth
  reference for `Effective DoF ratio`.

## Regression Suite Snapshot

- New P1 tests: 13 (sigma 7, p_res_eff 2, fisher 3, kappa 1).
- Pre-existing tests in `Tests/python/`: 3 files (cfar, detection_mc, grid_stats).
- Combined `pytest Tests/python/ -v` count: **23 passed**.
