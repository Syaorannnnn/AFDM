# P2: LRT-CLIP Double-Loop MM Framework — Design Spec

**Date:** 2026-05-11
**Status:** Design spec (brainstorming approved, not yet planned)
**Scope:** Journal Section III (Algorithm) — LRT-optimal path detector + BDLR-aware CFAR threshold + MM double-loop convergence
**Branch:** `feat/dev` (milestone merge point for `feat/p2-lrt-mm`)
**Extends:** `docs/superpowers/specs/2026-05-10-afdm-correlated-clip-research-design.md` §3.1 (D2) and §3.3 (D4)
**Consumes:** P1 deliverables `sigma_closed_form.closed_form_sigma`, `sigma_closed_form.bdlr_decompose`, `p_res_eff.p_res_eff_from_sigma`, `fisher_crlb.*`, `kappa_bound.kappa_weyl_upper_bound`

---

## 1. Research Positioning

**Academic target:** P2 produces the Section III (Algorithm) chapter of a TWC/TSP journal paper on AFDM-GI-Free receivers over BDLR channels. P1 supplied the numerical validation layer; P2 supplies the **actual algorithm** that consumes those validated primitives and the convergence theorem that makes it reviewer-defensible.

**What P2 delivers:**

1. `LRT-optimal path detector` (spec §3.1 D2) — Capon-style generalised matched filter with $\Sigma^{-1}$ whitening.
2. `BDLR-aware CFAR threshold` (spec §3.1 D4) — per-candidate adaptive threshold derived from the LRT denominator.
3. `MM double-loop convergence framework` (new in P2) — global negative log-likelihood objective with monotone-decrease guarantee, carrying Proposition 1 for the journal.
4. `End-to-end LRT-CLIP pipeline` — three-stage receiver `clip_lrt.py` running the MM loop, parallel to P1 baseline `clip_strict.py`, selectable via a single CLI flag.
5. `Full benchmark matrix` vs Naive-MF and Oracle-Σ-LRT baselines — four-metric suite over a 5 × 6 × 3 grid, all journal-ready plots.

**What P2 explicitly does NOT do** (deferred to P3/P4/journal-write-up):

- Machine-learning enhancements (unrolled EM, γ threshold calibration) → P3.
- 3GPP CDL channel benchmark → P4 + M3 validation.
- Block-CPSBL / AMF / Kronecker-CFAR SOTA comparison → future work.
- MIMO-AFDM extension → future work.

---

## 2. Mathematical Framework

### 2.1 LRT-Optimal Detector (D2)

Observation model per candidate $k \in \{1, \ldots, M_\mathrm{cand}\}$:

$$\mathcal{H}_0: \mathbf{r} = \mathbf{i}_\mathrm{other} + \mathbf{z}, \qquad \mathcal{H}_1: \mathbf{r} = h_k \boldsymbol{\Phi}_k \mathbf{x}_\mathrm{p} + \mathbf{i}_\mathrm{other} + \mathbf{z}$$

where $\mathbf{i}_\mathrm{other}, \mathbf{z}$ jointly follow $\mathcal{CN}(\mathbf{0}, \boldsymbol{\Sigma})$ and $\boldsymbol{\Sigma}$ is the residual covariance under the current support estimate.

The generalised likelihood ratio test over unknown $h_k$ yields the **Capon-style whitened matched filter**:

$$\boxed{T(\mathbf{r}; k) = \frac{|\boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \mathbf{r}|^2}{\boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k}}$$

SNR gain vs the naive matched filter $T_\mathrm{naive}(\mathbf{r}; k) = |\boldsymbol{\Phi}_k^\mathsf{H} \mathbf{r}|^2$:

$$\mathrm{SNR}_\mathrm{gain}(k) = 10 \log_{10}\!\left( \frac{ \boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k \cdot \boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma} \boldsymbol{\Phi}_k }{ \|\boldsymbol{\Phi}_k\|^4 } \right)$$

By Cauchy-Schwarz this gain is non-negative, zero iff $\boldsymbol{\Sigma} \propto \mathbf{I}$ (verified in P1 `test_p_res_eff_scalar_sigma_is_constant`). When $\rho \to 1$ and $\boldsymbol{\Sigma}$ approaches rank-deficiency, gains of 5-10 dB at $\mathrm{SNR} \geq 15$ dB are the target.

### 2.2 BDLR-Aware CFAR Threshold (D4)

Under $\mathcal{H}_0$, $\boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \mathbf{r}$ is $\mathcal{CN}(0, \boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k)$, so $T(\mathbf{r}; k)$ conditional on $\boldsymbol{\Sigma}$ follows a scaled chi-square with 2 DoF. With Bonferroni correction over $M_\mathrm{cand}$ candidates:

$$\boxed{\eta_\mathrm{new}(k) = \sigma_\mathrm{eff}^2 \cdot \boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k \cdot \ln\!\left(\frac{M_\mathrm{cand}}{P_\mathrm{FA}}\right)}$$

where $\sigma_\mathrm{eff}^2 = 1$ because the Capon normalisation already absorbs the Rayleigh quotient scale.

**Degeneration** to P1 `cfar.py` formula: when $\boldsymbol{\Sigma} = \alpha \mathbf{I}$, $\boldsymbol{\Phi}_k^\mathsf{H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k = \|\boldsymbol{\Phi}_k\|^2 / \alpha$ and the threshold reduces to $\alpha \cdot \ln(M / P_\mathrm{FA})$, matching the scalar-residual form.

### 2.3 MM Double-Loop Optimisation (new)

Global objective (Gaussian NLL + cluster-level group lasso):

$$\boxed{\mathcal{L}(\boldsymbol{\Sigma}, \mathbf{g}) = \underbrace{\mathbf{r}^\mathsf{H} \boldsymbol{\Sigma}^{-1} \mathbf{r} + \log\det\boldsymbol{\Sigma}}_{\mathcal{L}_\mathrm{NLL}} + \underbrace{\lambda_\mathrm{gl} \|\mathbf{g}\|_{2,1}}_{\mathcal{L}_\mathrm{pen}}}$$

subject to the BDLR-operator reconstruction

$$\boldsymbol{\Sigma}(\mathbf{g}) = N_0 \mathbf{I} + \sum_{c=1}^C |g_c|^2 \boldsymbol{\Phi}_c \boldsymbol{\Phi}_c^\mathsf{H}, \qquad \boldsymbol{\Phi}_c := \sum_{i \in c} \boldsymbol{\Phi}_i$$

This is exactly the coherent part returned by P1 `bdlr_decompose` when $\rho_c = 1$ per cluster, with the simplifying identification $|g_c|^2 \leftrightarrow \pi_c \cdot \sigma_c^2$ (where $\pi_c$ is the cluster-power weight and $\sigma_c^2$ is the per-path variance; we use $\pi_c$ here to avoid the overloaded symbol $P_c$ that P1 uses for both quantities in different sections).

**Double loop:**

- **Inner loop** (outer $\mathbf{g}$ fixed): closed-form MoM for $\boldsymbol{\Sigma}$ given current support; equivalent to block-coordinate MLE under Gaussian observation.
- **Outer loop** (inner $\boldsymbol{\Sigma}$ fixed): proximal gradient on $\mathbf{g}$ with group soft-threshold operator

$$g_c^{(t+1)} = \mathrm{prox}_{\lambda_\mathrm{gl}/L}(g_c^{(t)} - L^{-1} \nabla_{g_c} \mathcal{L}_\mathrm{NLL}), \qquad \mathrm{prox}_\mu(x) = \max(0, 1 - \mu/|x|) \cdot x$$

with step size $L$ = Lipschitz bound of $\nabla \mathcal{L}_\mathrm{NLL}$ (closed form via eigendecomposition of $\boldsymbol{\Sigma}^{-1}$).

**Proposition 1 (monotone decrease):** Let $\mathcal{Q}_\mathrm{in}$ and $\mathcal{Q}_\mathrm{out}$ be standard Hong-Luo MM majorisers for the inner and outer problems respectively. Then
$$\mathcal{L}^{(t+1)} \leq \mathcal{Q}_\mathrm{out}(\mathbf{g}^{(t+1)} | \mathbf{g}^{(t)}, \boldsymbol{\Sigma}^{(t+1)}) \leq \mathcal{L}^{(t)}$$
with equality iff $\mathbf{g}^{(t+1)} = \mathbf{g}^{(t)}$. Proof sketched in Appendix A of the journal paper; numerical evidence provided by `test_mm_lrt_loop::test_nll_monotonically_non_increasing` (20 seeds × 15 outer steps = 300 descent pairs).

**Stopping rule:** $|\mathcal{L}^{(t+1)} - \mathcal{L}^{(t)}| / |\mathcal{L}^{(t)}| < 10^{-4}$ with safeguards $T_\mathrm{in}^\mathrm{max} = 3$, $T_\mathrm{out}^\mathrm{max} = 15$.

**Regularisation weight $\lambda_\mathrm{gl}$:** initialised from a Bayesian Information Criterion (BIC) heuristic $\lambda_\mathrm{gl}^{(0)} = \sigma_N \sqrt{2 \ln M_\mathrm{cand}}$ (matching P1 `cfar.py` tail probability), then kept fixed per trial. A grid sweep $\lambda_\mathrm{gl} \in \{0.5 \lambda^{(0)}, \lambda^{(0)}, 2 \lambda^{(0)}\}$ runs once in the benchmark as an ablation. The Oracle-Σ baseline exposes $\lambda^\star$ (closed form under known support) as a reference curve in the ablation plot.

### 2.4 Numerical Safeguard

$\boldsymbol{\Sigma}^{-1}$ is always computed via Tikhonov-regularised eigendecomposition

$$\hat{\boldsymbol{\Sigma}}^{-1} = \mathbf{V} (\boldsymbol{\Lambda} + \epsilon \mathbf{I})^{-1} \mathbf{V}^\mathsf{H}, \qquad \epsilon = \max(N_0, \lambda_\mathrm{min}(\boldsymbol{\Sigma}) / \kappa^\star)$$

where $\kappa^\star$ is the P1 `kappa_weyl_upper_bound` threshold (default 100, configurable). When $\hat\kappa(\boldsymbol{\Sigma}) > \kappa^\star$ the loop triggers a fallback to the scalar CFAR detector for that iteration, preventing NLL explosion.

---

## 3. Code Structure

Five new modules under `Sim/experiments/py_correlated/`. **No modification** to any existing file; P1 baseline (`clip_strict.py` and friends) stays byte-exact.

### 3.1 Module Map

| Module | Responsibility | P1 dependencies | Line budget |
| --- | --- | --- | --- |
| `lrt_detector.py` | Pure function $T(\mathbf{r}; k)$; accepts $(\mathbf{r}, \boldsymbol{\Phi}_k\mathrm{-or-batch}, \boldsymbol{\Sigma})$, returns per-candidate LRT statistics | `sigma_closed_form` for type contract | ~120 |
| `cfar_lrt.py` | Per-candidate threshold $\eta_\mathrm{new}(k)$; scalar-Σ degeneration test; Bonferroni hook | `p_res_eff` (Rayleigh quotient) | ~100 |
| `sigma_mom.py` | MoM closed-form $\hat{g}_c$ from support + residual; Σ reconstruction from $\mathbf{g}$ | `bdlr_decompose` (operator form) | ~180 |
| `mm_lrt_loop.py` | `mm_double_loop(r, Phi_list, clusters, N0, cfg) -> MMResult`; inner MoM + outer proximal gradient; adaptive stopping + safeguard | `kappa_weyl_upper_bound`, `sigma_closed_form` | ~300 |
| `clip_lrt.py` | Three-stage LRT-CLIP (Coarse-Lock / Iteration / Polishing); parallel to `clip_strict.py` | all of above + `fisher_crlb` for diagnostics | ~500 |

Total new code ≈ 1200 lines, comparable to P1's ≈ 950 lines. All files ≤ 500 lines (CLAUDE.md rule).

### 3.2 API Contract (frozen for parallel development)

```python
# lrt_detector.py
def lrt_statistic(r, Phi_k, Sigma_inv) -> float:
    """T(r;k) for a single candidate using precomputed Sigma_inv."""

def lrt_statistic_batch(r, Phi_stack, Sigma_inv) -> ndarray:
    """(M,) per-candidate statistics for a batch of M candidates."""

def lrt_snr_gain(Phi_k, Sigma) -> float:
    """SNR gain in dB vs naive matched filter at this candidate."""

# cfar_lrt.py
def cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0) -> float:
    """eta_new(k); degenerates to scalar CFAR when Sigma proportional to I."""

def cfar_threshold_lrt_batch(Phi_stack, Sigma_inv, M_cand, P_fa) -> ndarray:
    """(M,) per-candidate thresholds."""

# sigma_mom.py
def estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support_mask) -> ndarray:
    """(C,) complex gains via method of moments closed form."""

def sigma_from_gains(gains, Phi_per_cluster, N0) -> ndarray:
    """Build (N, N) Sigma from estimated cluster gains."""

# mm_lrt_loop.py
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
    gains: ndarray          # (C,) complex
    support: ndarray        # (M_cand,) bool
    sigma: ndarray          # (N, N) complex
    sigma_inv: ndarray      # (N, N) complex (Tikhonov-regularised)
    nll_trajectory: list[float]
    n_inner: int
    n_outer: int
    converged: bool
    fallback_count: int     # times the kappa-safeguard fired

def mm_double_loop(r, Phi_list, cluster_labels, N0, cfg: MMConfig) -> MMResult:
    """Main driver. MM theorem numerically enforced via nll_trajectory."""

# clip_lrt.py
def clip_lrt_receive(y, pilot_frame, config) -> dict:
    """End-to-end three-stage LRT receiver. Diagnostic dict matches clip_strict."""
```

### 3.3 Baseline Toggle

`runner.py` (existing) gains a `--detector {naive,lrt,oracle_lrt}` CLI flag. The baseline of record is **`naive` = byte-exact P1** (`clip_strict.receive`); `lrt` routes to `clip_lrt.clip_lrt_receive`; `oracle_lrt` injects ground-truth Σ from the channel sampler (used only in benchmarks, not a shipped algorithm). One regression test (`test_clip_lrt_integration::test_naive_mode_matches_p1_baseline`) locks this contract.

### 3.4 Test Layout (mirrors P1 convention)

| Test file | Count | Key assertions |
| --- | --- | --- |
| `test_lrt_detector.py` | 3 | scalar-Σ → naive MF equivalence (SNR gain 0 dB); rank-1 Σ produces expected whitening; Oracle-Σ matches hand-computed reference |
| `test_cfar_lrt.py` | 2 | scalar-Σ degenerates to P1 cfar formula; MC-measured $\hat P_\mathrm{FA}$ matches nominal at scalar Σ |
| `test_sigma_mom.py` | 3 | Σ = αI degeneration; MC bias/variance on $\hat\rho_c$ under oracle support; rank-collapse boundary does not explode |
| `test_mm_lrt_loop.py` | 4 | NLL monotonic non-increase (20 seeds × 15 outer steps); stopping rule fires correctly; adaptive safeguard engages on ill-conditioned Σ; Oracle-Σ initial condition recovers true support |
| `test_clip_lrt_integration.py` | 2 | single-trial end-to-end return structure; naive-mode byte-exact match with P1 baseline |

**Total: 14 new tests**, matching P1's cadence.

---

## 4. Experimental Matrix

### 4.1 Baselines

| Name | Definition | Role | Build cost |
| --- | --- | --- | --- |
| Naive-MF | $T = \lvert\boldsymbol{\Phi}_k^\mathsf{H} \mathbf{r}\rvert^2$, scalar CFAR — P1 `clip_strict` | Lower bound | 0 (reuse) |
| Oracle-Σ-LRT | LRT detector with ground-truth $\boldsymbol{\Sigma}^\star$ injected from channel sampler | Upper bound | ~40 lines fixture |
| **Proposed MM-LRT** | full `clip_lrt.py` with MM double loop | core deliverable | 1200 lines |

SOTA peers (Block-CPSBL, Kronecker-CFAR, AMF) are **not** built here; the journal paper will survey them in Section II and cite published numbers, deferring head-to-head benchmarks to future work.

### 4.2 Grid

- $\rho \in \{0.0, 0.3, 0.5, 0.7, 0.9\}$ (5 points)
- SNR $\in \{0, 5, 10, 15, 20, 25\}$ dB (6 points)
- Cluster profile $\in \{[2,2],\ [3,3,2],\ [5]\}$ (3 profiles: sparse / 3GPP-CDL-like / single-cluster-full-coherence)

$5 \times 6 \times 3 = 90$ grid cells × 3 detectors × 30 trials = **8100 trials**. Budget: ~6 h wall clock with 4-worker joblib on a commodity laptop.

### 4.3 Metrics (four-tuple)

1. **BER** (main headline metric).
2. **Cluster Pd** = $|\hat{\mathcal{S}}_c \cap \mathcal{S}_c^\star| / |\mathcal{S}_c^\star|$ averaged over clusters.
3. **NMSE** = $\|\hat{\mathbf{H}}_\mathrm{eff} - \mathbf{H}_\mathrm{eff}^\star\|_F^2 / \|\mathbf{H}_\mathrm{eff}^\star\|_F^2$.
4. **SNR gain** = decibel advantage of Proposed vs Naive-MF at matched BER target (or at matched detector threshold, whichever plot is clearer).

### 4.4 Output Plots

- Main plots (paper-ready, 4): BER vs SNR at $\rho=0.9$ / SNR gain vs $\rho$ at 15 dB / Cluster Pd vs SNR at $\rho=0.7$ / NMSE vs $\rho$ at 20 dB.
- Diagnostic plots (2): MM NLL trajectories (20 randomly chosen seeds from the 30-trial pool, overlaid) / $\kappa(\boldsymbol{\Sigma})$ distribution per grid cell.
- Per-cell JSON + NLL trajectory CSV saved to `Sim/experiments/Results/p2/<timestamp>/`.

---

## 5. Deliverables Checklist

| ID | Content | Spec tie-in |
| --- | --- | --- |
| E1 | `lrt_detector.py` + 3 tests | spec §3.1 D2 |
| E2 | `cfar_lrt.py` + 2 tests | spec §3.1 D4 |
| E3 | `sigma_mom.py` + 3 tests | P2 internal (MoM block) |
| E4 | `mm_lrt_loop.py` + 4 tests (incl. mandatory monotone-NLL test) | P2 internal (MM double loop) |
| E5 | `clip_lrt.py` + 2 integration tests (naive byte-exact lock included) | P2 end-to-end |
| E6 | `exp_p2_benchmark.py` sweep + JSON/CSV artefacts | P2 benchmark |
| E7 | `docs/diagnostics/2026-MM-DD-p2-journal-section-notes.md` (gate values, main/diag plots described, MM convergence report) | closing diagnostic |

---

## 6. Risks

| Risk | Likelihood | Mitigation |
| --- | :---: | --- |
| MM monotonicity violated by float roundoff at $\rho \to 1$ | Medium | `kappa_weyl_upper_bound` sentinel enlarges Tikhonov $\epsilon$; Prop 1 footnoted as "up to $O(\epsilon_\mathrm{mach})$" |
| $\lambda_\mathrm{gl}$ hard to tune → poor support recovery | Medium | BIC auto-select on initial call; Oracle-Σ baseline gives $\lambda^\star$ reference curve |
| Benchmark runtime overruns budget | Low | Two-stage sweep: halved $\rho$/profile in first pass, full grid once pipeline is frozen |
| Low-$\rho$ SNR gain < 0.5 dB (expected per P1 §3.1) | Low | Frame story around $\rho \geq 0.7$ regime; 0 dB at $\rho=0$ is a feature not a bug |
| Group-lasso prox op bugs | Medium | Oracle-Σ gives closed-form $\mathbf{g}^\star$ ground truth; unit test asserts recovery within 1e-6 |
| `clip_lrt.py` integration regresses P1 baseline | Medium | `--detector naive` regression test forces byte-exact match |
| MM double loop has bug that saves NLL but not support | Low | Test checks both NLL trajectory AND support stability post-stop |

---

## 7. Open Questions (for journal write-up, non-blocking)

1. Is there a closed-form optimal $\lambda_\mathrm{gl}$ (e.g. Stein-style BIC bound)?
2. Can the MM double loop be shown **globally** optimal in a low-$\rho$ / low-SNR regime (not just monotone)?
3. Does the LRT detector attain the P1 CRLB asymptotically? If yes, note as efficiency theorem; if not, quantify the gap.
4. Formal reduction relation to Hong-Luo 2017 SBL — is MM-LRT a special case or a proper generalisation?
5. Interaction with P1 Section §3.3 $\kappa^\star$ threshold: is there an optimal $\kappa^\star(\mathrm{SNR}, \rho)$ that minimises fallback count while preserving convergence?

---

## 8. Timeline

P2 decomposes into ~10 subagent-driven tasks, 2-5 min per step, following the P1 cadence. Estimated wall-clock 2-3 weeks including benchmark runs and two-stage review per task. Success criteria:

- All 14 tests green.
- MM NLL monotone across 300 descent pairs.
- `--detector naive` byte-exact match P1 baseline.
- Benchmark JSON + CSV artefacts produced for all 90 grid cells.
- Diagnostic note records gate values and four main plots.

On completion, P2 merges back to `feat/dev` via `--no-ff` as the second research milestone.
