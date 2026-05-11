# P2 LRT-CLIP MM Framework — Diagnostic Notes

**Date:** 2026-05-11
**Plan:** `docs/superpowers/plans/2026-05-11-p2-lrt-clip-mm-framework.md`
**Spec:** `docs/superpowers/specs/2026-05-11-p2-lrt-clip-mm-framework-design.md`

## Gate Values

The benchmark sweep at `Sim/experiments/Results/p2/summary_<timestamp>.json` passes iff:

- `lrt` detector achieves lower BER than `naive` detector at (rho>=0.7, SNR>=15 dB) for at least 80% of cluster profiles.
- `lrt` NMSE (dB) is lower than `naive` NMSE at rho=0.9 for at least 2 of 3 cluster profiles.
- MM NLL trajectories are monotone non-increasing for all 20 test seeds in
  `test_mm_lrt_loop::test_nll_monotonically_non_increasing` (0 violations).
- `lrt_snr_gain` reported by `lrt_detector.lrt_snr_gain` is >= 0 dB for every tested Sigma.

## Observed Values (2026-05-11 sweep)

Artefact: `Sim/experiments/Results/p2/summary_20260511-115545.json` (180 entries,
5 rho x 6 SNR x 3 profile x 2 detector x 30 trials).

- LRT BER advantage at rho=0.9, SNR=15 dB, profile=[3,3,2]: delta_ber = +0.3889
  (naive BER 0.8667 vs LRT BER 0.4778; LRT improves by 38.9 pp, cluster Pd
  jumps from 0.333 to 1.000).
- BER-gate coverage (spec: "LRT beats naive at rho>=0.7, SNR>=15 dB for at least 80% of cluster profiles"):
  - Per-cell reading (15/18 = 83.3%): PASS
  - Per-profile reading (2/3 = 66.7%): FAIL — only profiles [2,2] and [3,3,2] show LRT winning at every cell in the rho>=0.7, SNR>=15 dB slice; profile [5] loses 3 of 6 cells
- Which reading is authoritative is ambiguous in the plan text; per-cell is reported for completeness, per-profile is the stricter interpretation
- LRT BER-gate coverage at rho>=0.7, SNR>=15 dB: LRT wins in 15 of 18 cells
  (83.3%), above the 80% target.
- LRT NMSE advantage at rho=0.9, averaged over profiles: delta_nmse_db = -0.2979
  (naive actually edges LRT on NMSE; per-profile deltas: [2,2]=-0.1267,
  [3,3,2]=-0.6850, [5]=-0.0821). LRT loses the 2-of-3 NMSE gate on this
  synthetic observation model even though it dominates BER and Pd. See "Known
  Limitations" below; the MoM Sigma estimator is unbiased but high-variance at
  N=24, C in {1,2,3}, which inflates LRT whitening residual relative to the
  naive projector that ignores correlation. Task 10 / P3 tuning will revisit.
- MM monotone NLL descent pairs observed: total_pairs = 43 / 0 violations
  (from `test_mm_lrt_loop::test_nll_monotonically_non_increasing`, 20 seeds x
  up to 3 consecutive-pair checks per trajectory).
- Fallback count (kappa-safeguard): not persisted (field exists in MMResult
  but not exported by exp_p2_benchmark).

Auxiliary empirical values already captured in unit tests:

- LRT SNR gain on rank-1-coherent + 0.05*I Sigma at N=32: 15.7 dB
  (`test_lrt_snr_gain_positive_for_correlated_sigma`, seed 102).
- MoM MC relative bias at N=24, C=2, n_trials=400: 0.15
  (`test_mom_estimator_recovers_gains_within_MC_bias`, seed 302).
- Oracle-Sigma |g| triple at high SNR: (0.718, 0.185, 0.521), ratio 0.355
  (`test_oracle_sigma_recovers_support_at_high_snr`, seed 422).

## Known Limitations

- Naive baseline here is a simplified top-half cluster selector, not the full
  P1 `clip_strict.py` pipeline. The spec 3.3 "byte-exact match" hook is deferred to
  the future integration task that wires `clip_lrt` into the real `runner.py`.
- Benchmark uses synthetic coherent-BDLR observations, not the full GiFree frame;
  downstream integration task will connect to `runner.py --detector lrt`.
- `lambda_gl` is fixed at 0.05 in all runs; the BIC-auto-select hook listed in
  spec 2.3 is reserved for P3 sensitivity ablation.
- NMSE regression vs. naive at rho=0.9 (delta_nmse_db = -0.30 dB average) is a
  known artefact of MoM Sigma variance on short coherent blocks; LRT still
  dominates on detection metrics (BER delta = +0.39 at the canonical cell,
  Pd 0.333 -> 1.000). Full integration into `clip_strict` where the LRT output
  feeds a downstream equaliser should recover NMSE parity; revisit in Task 10.
- Profile [5] (single cluster spanning all five paths, full coherence) is the consistent BER loser: every BER-gate failure above rho=0.7 is a [5] cell. Mechanism: with a single full-support cluster, the naive top-half selector reduces to the identity projector, so LRT's whitening advantage disappears and MoM variance dominates. Resolution is a P3 open question — either a different `_naive_receive` baseline for single-cluster profiles or a dedicated single-cluster branch in `clip_lrt`

## Downstream Consumers

- P3 (ML-enhanced LRT) will reuse `mm_lrt_loop` as the MM teacher and `clip_lrt`
  as the unrolled-EM structural template.
- P4 (3GPP CDL benchmark) will swap the synthetic observation generator for a
  CDL-backed channel sampler while keeping `clip_lrt_receive` unchanged.
