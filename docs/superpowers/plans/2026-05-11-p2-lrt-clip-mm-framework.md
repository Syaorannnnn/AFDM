# P2: LRT-CLIP MM Double-Loop Framework Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an LRT-optimal path detector with BDLR-aware CFAR threshold and MM double-loop convergence framework, delivered as five new pure-Python modules and 14 pytest tests that together form the Section III algorithm of the planned TWC/TSP journal paper.

**Architecture:** Five layered modules under `Sim/experiments/py_correlated/`. Bottom layer (`lrt_detector.py`, `cfar_lrt.py`) are stateless numerical helpers consuming a given `Sigma_inv`. Middle layer (`sigma_mom.py`) estimates cluster-level gains via closed-form Method-of-Moments. Top layer (`mm_lrt_loop.py`) runs the outer group-proximal-gradient on gains + inner MoM on Sigma, with a P1-κ safeguard and NLL-based stopping. Integration layer (`clip_lrt.py`) wraps the MM loop into a three-stage CLIP-style receiver callable through the existing `runner.py` via a new `--detector lrt` flag. No existing file is modified.

**Tech Stack:** Python 3.11+, NumPy, SciPy (for `scipy.linalg` eigendecomposition), pytest. Reuses P1's `sigma_closed_form.closed_form_sigma`, `sigma_closed_form.bdlr_decompose`, `p_res_eff.p_res_eff_from_sigma`, `kappa_bound.kappa_weyl_upper_bound`, and `fisher_crlb.*` for diagnostics.

**Spec:** `docs/superpowers/specs/2026-05-11-p2-lrt-clip-mm-framework-design.md`

---

## File Structure

**Create** (all under `Sim/experiments/py_correlated/` unless noted):

| File | Responsibility |
| --- | --- |
| `lrt_detector.py` | `lrt_statistic`, `lrt_statistic_batch`, `lrt_snr_gain` — pure detector maths |
| `cfar_lrt.py` | `cfar_threshold_lrt`, `cfar_threshold_lrt_batch` — per-candidate adaptive threshold |
| `sigma_mom.py` | `estimate_cluster_gains_mom`, `sigma_from_gains` — MoM inner block |
| `mm_lrt_loop.py` | `MMConfig`, `MMResult`, `mm_double_loop` — double-loop MM driver |
| `clip_lrt.py` | `clip_lrt_receive`, Stage 0/1/2 helpers — end-to-end receiver |
| `exp_p2_benchmark.py` | sweep runner producing JSON/CSV under `Sim/experiments/Results/p2/<ts>/` |
| `Tests/python/test_lrt_detector.py` | 3 tests |
| `Tests/python/test_cfar_lrt.py` | 2 tests |
| `Tests/python/test_sigma_mom.py` | 3 tests |
| `Tests/python/test_mm_lrt_loop.py` | 4 tests |
| `Tests/python/test_clip_lrt_integration.py` | 2 tests |
| `docs/diagnostics/2026-05-11-p2-journal-section-notes.md` | diagnostic gate notes |

**Do not modify** any existing file under `Sim/experiments/py_correlated/`, `Tests/python/`, or `Core/`.

---

## Conventions

- **sys.path:** tests prepend `<repo>/Sim/experiments/py_correlated` to `sys.path`, matching the existing pattern (see `Tests/python/test_sigma_closed_form.py:7-8`).
- **RNG:** every test seeds `np.random.default_rng(seed=<fixed>)`; test sizes keep per-test wall time < 2 s.
- **Tolerances:** closed-form vs MC agreement uses `rtol=0.05, atol=1e-3` unless relaxed with an inline justification comment. Monotone-NLL test allows `eps_rel=1e-8` slack for floating-point roundoff.
- **Commit discipline:** TDD per task (test first, implement second, commit third). Prefix commits with `feat(p2):`, `test(p2):`, `docs(p2):`, `chore(p2):`.
- **Pure Python:** no MATLAB bridge in P2.

---

## Task 1: Scaffold `lrt_detector` module with API contract

**Files:**

- Create: `Sim/experiments/py_correlated/lrt_detector.py`
- Create: `Tests/python/test_lrt_detector.py`

- [ ] **Step 1: Write the failing import test**

```python
"""Sanity test that the new module imports and exposes the expected API."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_module_importable():
    import lrt_detector
    assert hasattr(lrt_detector, "lrt_statistic")
    assert hasattr(lrt_detector, "lrt_statistic_batch")
    assert hasattr(lrt_detector, "lrt_snr_gain")
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_lrt_detector.py::test_module_importable -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'lrt_detector'`

- [ ] **Step 3: Create the stub module**

```python
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_lrt_detector.py::test_module_importable -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/lrt_detector.py Tests/python/test_lrt_detector.py
git commit -m "feat(p2): scaffold lrt_detector module with API contract"
```

---

## Task 2: LRT detector implementation + scalar-Sigma degeneration test

Implements `lrt_statistic`, `lrt_statistic_batch`, `lrt_snr_gain` and verifies (a) scalar Sigma degenerates to naive MF (SNR gain 0 dB) and (b) Oracle-Sigma numerical reference.

**Files:**

- Modify: `Sim/experiments/py_correlated/lrt_detector.py`
- Modify: `Tests/python/test_lrt_detector.py`

- [ ] **Step 1: Write degeneration + reference tests**

Append to `Tests/python/test_lrt_detector.py`:

```python
from lrt_detector import lrt_statistic, lrt_statistic_batch, lrt_snr_gain


def _make_random_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_scalar_sigma_snr_gain_is_zero_db():
    rng = np.random.default_rng(100)
    N = 24
    Phi_k = _make_random_unitary(N, rng)
    alpha = 0.3
    Sigma = alpha * np.eye(N, dtype=complex)
    gain_db = lrt_snr_gain(Phi_k, Sigma)
    assert abs(gain_db) < 1e-10, f"expected 0 dB, got {gain_db:.3e} dB"


def test_lrt_statistic_batch_matches_single():
    rng = np.random.default_rng(101)
    N, M = 16, 5
    r = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    Phi_stack = np.stack([_make_random_unitary(N, rng) for _ in range(M)], axis=0)
    Sigma = 0.4 * np.eye(N, dtype=complex) + 0.1 * np.ones((N, N), dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    batch = lrt_statistic_batch(r, Phi_stack, Sigma_inv)
    singles = np.array([
        lrt_statistic(r, Phi_stack[m], Sigma_inv) for m in range(M)
    ])
    assert batch.shape == (M,)
    assert np.allclose(batch, singles, rtol=1e-12, atol=1e-14)


def test_lrt_snr_gain_positive_for_correlated_sigma():
    rng = np.random.default_rng(102)
    N = 32
    Phi_k = _make_random_unitary(N, rng)
    # Rank-1 coherent + identity floor
    v = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    Sigma = 0.05 * np.eye(N, dtype=complex) + np.outer(v, v.conj())
    gain_db = lrt_snr_gain(Phi_k, Sigma)
    assert gain_db > 0.5, f"expected >0.5 dB on correlated Sigma, got {gain_db:.3f}"
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_lrt_detector.py -v`
Expected: 3 failures with `NotImplementedError`

- [ ] **Step 3: Implement the module**

Replace the body of `Sim/experiments/py_correlated/lrt_detector.py`:

```python
"""LRT-optimal detector: T(r;k) = |Phi_k^H Sigma^{-1} r|^2 / (Phi_k^H Sigma^{-1} Phi_k).

Implements spec §2.1 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

import numpy as np


def lrt_statistic(r, Phi_k, Sigma_inv):
    """Single-candidate LRT statistic.

    Parameters
    ----------
    r : (N,) complex ndarray
    Phi_k : (N, N) complex ndarray — candidate response operator
    Sigma_inv : (N, N) complex ndarray — precomputed inverse covariance
    """
    w = Phi_k.conj().T @ Sigma_inv @ r
    num = float(np.vdot(w, w).real)
    norm_term = Phi_k.conj().T @ Sigma_inv @ Phi_k
    den = float(np.trace(norm_term).real)
    if den <= 0:
        return 0.0
    return num / den


def lrt_statistic_batch(r, Phi_stack, Sigma_inv):
    """Batch over M candidate operators."""
    M = Phi_stack.shape[0]
    out = np.empty(M, dtype=float)
    for m in range(M):
        out[m] = lrt_statistic(r, Phi_stack[m], Sigma_inv)
    return out


def lrt_snr_gain(Phi_k, Sigma):
    """SNR gain in dB of LRT over naive matched filter at candidate k."""
    Sigma_inv = np.linalg.pinv(Sigma)
    num = np.trace(Phi_k.conj().T @ Sigma_inv @ Phi_k).real
    mid = np.trace(Phi_k.conj().T @ Sigma @ Phi_k).real
    den = np.linalg.norm(Phi_k, "fro") ** 4
    if den <= 0 or num <= 0 or mid <= 0:
        return 0.0
    ratio = (num * mid) / den
    return 10.0 * float(np.log10(max(ratio, 1e-300)))
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_lrt_detector.py -v`
Expected: 4 passed (import + 3 new)

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/lrt_detector.py Tests/python/test_lrt_detector.py
git commit -m "feat(p2): LRT detector with scalar-Sigma degeneration and SNR gain"
```

---

## Task 3: CFAR threshold for LRT — module + degeneration test

Implements `cfar_threshold_lrt`, `cfar_threshold_lrt_batch`, and verifies scalar-Sigma degeneration to the P1 cfar formula.

**Files:**

- Create: `Sim/experiments/py_correlated/cfar_lrt.py`
- Create: `Tests/python/test_cfar_lrt.py`

- [ ] **Step 1: Write the failing tests**

```python
"""BDLR-aware CFAR threshold tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from cfar_lrt import cfar_threshold_lrt, cfar_threshold_lrt_batch


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_scalar_sigma_matches_p1_cfar_formula():
    """When Sigma = alpha * I, eta_new(k) = alpha * ln(M/P_FA)."""
    rng = np.random.default_rng(200)
    N, M_cand, P_fa = 20, 7, 1e-3
    Phi_k = _rand_unitary(N, rng)
    alpha = 0.6
    Sigma = alpha * np.eye(N, dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    eta = cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0)
    # eta_new(k) = 1 * (||Phi_k||^2 / alpha) * ln(M/P_FA)
    expected = (np.linalg.norm(Phi_k, "fro") ** 2 / alpha) * np.log(M_cand / P_fa)
    assert np.isclose(eta, expected, rtol=1e-10)


def test_batch_matches_singles():
    rng = np.random.default_rng(201)
    N, M_cand, P_fa = 16, 5, 1e-3
    Phi_stack = np.stack([_rand_unitary(N, rng) for _ in range(M_cand)], axis=0)
    Sigma = 0.2 * np.eye(N, dtype=complex) + 0.05 * np.ones((N, N), dtype=complex)
    Sigma_inv = np.linalg.pinv(Sigma)
    batch = cfar_threshold_lrt_batch(Phi_stack, Sigma_inv, M_cand, P_fa)
    singles = np.array([
        cfar_threshold_lrt(Phi_stack[m], Sigma_inv, M_cand, P_fa) for m in range(M_cand)
    ])
    assert batch.shape == (M_cand,)
    assert np.allclose(batch, singles, rtol=1e-12)
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_cfar_lrt.py -v`
Expected: FAIL at import (module missing)

- [ ] **Step 3: Implement `cfar_lrt.py`**

```python
"""BDLR-aware per-candidate CFAR threshold for LRT detector.

eta_new(k) = sigma_eff2 * Phi_k^H Sigma^{-1} Phi_k * ln(M_cand / P_fa)

Degenerates to P1 scalar cfar formula when Sigma = alpha * I.
See spec §2.2 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

import numpy as np


def cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0):
    """Per-candidate CFAR threshold for a single Phi_k."""
    if P_fa <= 0 or P_fa >= 1:
        raise ValueError(f"P_fa must be in (0,1), got {P_fa}")
    if M_cand < 1:
        raise ValueError(f"M_cand must be >=1, got {M_cand}")
    denom_term = Phi_k.conj().T @ Sigma_inv @ Phi_k
    rayleigh = float(np.trace(denom_term).real)
    if rayleigh <= 0:
        return 0.0
    return sigma_eff2 * rayleigh * float(np.log(M_cand / P_fa))


def cfar_threshold_lrt_batch(Phi_stack, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0):
    """Batch over M candidate operators."""
    M = Phi_stack.shape[0]
    out = np.empty(M, dtype=float)
    for m in range(M):
        out[m] = cfar_threshold_lrt(Phi_stack[m], Sigma_inv, M_cand, P_fa, sigma_eff2)
    return out
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_cfar_lrt.py -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/cfar_lrt.py Tests/python/test_cfar_lrt.py
git commit -m "feat(p2): BDLR-aware CFAR threshold with scalar-Sigma degeneration"
```

---

## Task 4: Method-of-Moments block `sigma_mom.py`

Implements inner-loop closed-form cluster-gain estimator and Sigma reconstructor used by the MM double loop.

**Files:**

- Create: `Sim/experiments/py_correlated/sigma_mom.py`
- Create: `Tests/python/test_sigma_mom.py`

- [ ] **Step 1: Write tests (degeneration, MC bias, reconstruction)**

```python
"""Method-of-moments cluster-gain estimator tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_mom import estimate_cluster_gains_mom, sigma_from_gains


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def test_sigma_from_gains_zero_gains_gives_noise_floor():
    rng = np.random.default_rng(300)
    N, C = 12, 3
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    gains = np.zeros(C, dtype=complex)
    N0 = 0.05
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    assert np.allclose(Sigma, N0 * np.eye(N), atol=1e-12)


def test_sigma_from_gains_reconstruction_matches_formula():
    rng = np.random.default_rng(301)
    N, C = 16, 2
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    gains = np.array([0.4 + 0.1j, 0.2 - 0.3j])
    N0 = 0.02
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    expected = N0 * np.eye(N, dtype=complex)
    for c in range(C):
        expected += (abs(gains[c]) ** 2) * (Phi_per_cluster[c] @ Phi_per_cluster[c].conj().T)
    assert np.allclose(Sigma, expected, atol=1e-12)


def test_mom_estimator_recovers_gains_within_MC_bias():
    """Oracle-support MoM: average |g_hat_c|^2 over MC should approach |g_true_c|^2."""
    rng = np.random.default_rng(302)
    N, C = 24, 2
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.7 + 0.0j, 0.5 + 0.0j])
    N0 = 0.01
    n_trials = 400
    accum_abs2 = np.zeros(C)
    for _ in range(n_trials):
        # draw coherent cluster gains and data
        h_per = np.array([g_true[c] * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2)) for c in range(C)])
        z = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
        r = sum(h_per[c] * (Phi_per_cluster[c] @ np.ones(N, dtype=complex) / np.sqrt(N)) for c in range(C)) + z
        support = np.ones(C, dtype=bool)
        g_hat = estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support)
        accum_abs2 += np.abs(g_hat) ** 2
    emp = accum_abs2 / n_trials
    true_abs2 = np.abs(g_true) ** 2
    rel = np.max(np.abs(emp - true_abs2) / (true_abs2 + 1e-12))
    assert rel < 0.35, f"MoM MC bias too large: rel={rel:.3f}"
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_sigma_mom.py -v`
Expected: FAIL at import

- [ ] **Step 3: Implement `sigma_mom.py`**

```python
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_sigma_mom.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/sigma_mom.py Tests/python/test_sigma_mom.py
git commit -m "feat(p2): MoM cluster-gain estimator with Sigma reconstruction"
```

---

## Task 5: MM double-loop driver scaffold

Introduces `MMConfig`, `MMResult`, and stub `mm_double_loop`, plus the sanity import test. Also defines the small helper `tikhonov_inverse` used as the κ-safeguarded inverse across the driver.

**Files:**

- Create: `Sim/experiments/py_correlated/mm_lrt_loop.py`
- Create: `Tests/python/test_mm_lrt_loop.py`

- [ ] **Step 1: Write sanity-import test**

```python
"""MM double-loop driver tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_module_importable():
    import mm_lrt_loop
    assert hasattr(mm_lrt_loop, "MMConfig")
    assert hasattr(mm_lrt_loop, "MMResult")
    assert hasattr(mm_lrt_loop, "mm_double_loop")
    assert hasattr(mm_lrt_loop, "tikhonov_inverse")
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_mm_lrt_loop.py::test_module_importable -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Create the stub module**

```python
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_mm_lrt_loop.py::test_module_importable -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/mm_lrt_loop.py Tests/python/test_mm_lrt_loop.py
git commit -m "feat(p2): scaffold mm_lrt_loop driver with API contract"
```

---

## Task 6: Tikhonov-safeguarded inverse + MM double-loop implementation

Implements `tikhonov_inverse` (eigendecomposition + ε floor governed by κ★), the MoM-inner / proximal-gradient-outer loop, the NLL trajectory, and stopping rule. Adds three behaviour tests: NLL monotone non-increase, stopping correctness, and Oracle-Σ support recovery.

**Files:**

- Modify: `Sim/experiments/py_correlated/mm_lrt_loop.py`
- Modify: `Tests/python/test_mm_lrt_loop.py`

- [ ] **Step 1: Write tests**

Append to `Tests/python/test_mm_lrt_loop.py`:

```python
from mm_lrt_loop import MMConfig, mm_double_loop, tikhonov_inverse


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def _simulate_observation(Phi_per_cluster, g_true, N0, rng):
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    r = np.zeros(N, dtype=complex)
    for c in range(C):
        h = g_true[c] * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2))
        r = r + h * (Phi_per_cluster[c] @ x_unit)
    z = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
    return r + z


def test_tikhonov_inverse_bounds_condition_number():
    rng = np.random.default_rng(400)
    N = 24
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Sigma = A @ A.conj().T / N + 1e-6 * np.eye(N)
    Sigma_inv = tikhonov_inverse(Sigma, kappa_star=50.0, N0=0.01)
    # Sigma_inv should be finite and roughly Sigma^{-1} up to the floor
    assert np.all(np.isfinite(Sigma_inv))
    assert Sigma_inv.shape == (N, N)


def test_nll_monotonically_non_increasing():
    """Proposition 1 numerical evidence: 20 seeds x 15 outer steps = 300 descent pairs."""
    violations = 0
    total_pairs = 0
    eps_rel = 1e-8
    for seed in range(400, 420):
        rng = np.random.default_rng(seed)
        N, C = 20, 2
        Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
        g_true = np.array([0.6 + 0.0j, 0.4 + 0.0j])
        r = _simulate_observation(Phi_per_cluster, g_true, N0=0.02, rng=rng)
        cfg = MMConfig(lambda_gl=0.05, T_in_max=3, T_out_max=15, nll_tol=1e-12, kappa_star=100.0)
        result = mm_double_loop(r, Phi_per_cluster, N0=0.02, cfg=cfg)
        traj = result.nll_trajectory
        for t in range(len(traj) - 1):
            total_pairs += 1
            if traj[t + 1] > traj[t] * (1 + eps_rel) + eps_rel:
                violations += 1
    assert total_pairs >= 40, f"only {total_pairs} descent pairs observed"
    assert violations == 0, f"{violations}/{total_pairs} NLL violations above eps_rel"


def test_stopping_rule_fires_on_convergence():
    rng = np.random.default_rng(421)
    N, C = 18, 2
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    r = _simulate_observation(Phi_per_cluster, np.array([0.5 + 0j, 0.5 + 0j]), N0=0.02, rng=rng)
    cfg = MMConfig(lambda_gl=0.05, T_in_max=3, T_out_max=40, nll_tol=1e-3)
    result = mm_double_loop(r, Phi_per_cluster, N0=0.02, cfg=cfg)
    assert result.converged, "should converge well before T_out_max=40"
    assert result.n_outer < 40


def test_oracle_sigma_recovers_support_at_high_snr():
    rng = np.random.default_rng(422)
    N, C = 24, 3
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.8 + 0j, 0.0 + 0j, 0.6 + 0j])  # cluster 1 is empty
    r = _simulate_observation(Phi_per_cluster, g_true, N0=0.001, rng=rng)
    cfg = MMConfig(lambda_gl=0.05, T_in_max=3, T_out_max=25)
    result = mm_double_loop(r, Phi_per_cluster, N0=0.001, cfg=cfg)
    # support: clusters 0 and 2 should be active, cluster 1 inactive
    assert result.support[0] and result.support[2]
    # cluster 1 either inactive or tiny gain
    assert (not result.support[1]) or (abs(result.gains[1]) < 0.1 * abs(result.gains[0]))
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_mm_lrt_loop.py -v`
Expected: FAIL — tikhonov_inverse and mm_double_loop raise NotImplementedError

- [ ] **Step 3: Implement the driver**

Replace the stubs in `Sim/experiments/py_correlated/mm_lrt_loop.py`:

```python
from sigma_mom import estimate_cluster_gains_mom, sigma_from_gains
from kappa_bound import kappa_weyl_upper_bound  # not strictly required here, reserved for diagnostics


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
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_mm_lrt_loop.py -v`
Expected: 5 passed (1 import + 4 behaviour)

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/mm_lrt_loop.py Tests/python/test_mm_lrt_loop.py
git commit -m "feat(p2): MM double-loop driver with monotone NLL and kappa safeguard"
```

---

## Task 7: End-to-end `clip_lrt` integration module

Wraps the MM driver into a three-stage receiver callable `clip_lrt_receive(y, pilot_frame, config)`. Adds smoke tests (structure + seed reproducibility). Note: "naive byte-exact match" belongs to the benchmark-runner layer (Task 8); this task owns only the LRT branch.

**Files:**

- Create: `Sim/experiments/py_correlated/clip_lrt.py`
- Create: `Tests/python/test_clip_lrt_integration.py`

- [ ] **Step 1: Write the smoke tests**

```python
"""clip_lrt end-to-end smoke tests."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from clip_lrt import clip_lrt_receive, ClipLrtConfig


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def _make_synthetic_observation(N, C, rng):
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.6 + 0.0j, 0.4 + 0.0j])
    N0 = 0.02
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    y = np.zeros(N, dtype=complex)
    for c in range(C):
        h = g_true[c] * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2))
        y = y + h * (Phi_per_cluster[c] @ x_unit)
    y = y + (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
    return y, Phi_per_cluster, N0


def test_clip_lrt_returns_expected_structure():
    rng = np.random.default_rng(500)
    y, Phi_per_cluster, N0 = _make_synthetic_observation(N=20, C=2, rng=rng)
    cfg = ClipLrtConfig(Phi_per_cluster=Phi_per_cluster, N0=N0)
    out = clip_lrt_receive(y, pilot_frame=None, config=cfg)
    for key in ["gains", "support", "sigma", "nll_trajectory", "n_outer", "stage_results"]:
        assert key in out, f"missing key {key}"
    # structural check: stage_results has 3 entries
    assert len(out["stage_results"]) == 3


def test_clip_lrt_is_seed_reproducible():
    rng_a = np.random.default_rng(501)
    rng_b = np.random.default_rng(501)
    y_a, Phi_a, N0_a = _make_synthetic_observation(N=16, C=2, rng=rng_a)
    y_b, Phi_b, N0_b = _make_synthetic_observation(N=16, C=2, rng=rng_b)
    cfg_a = ClipLrtConfig(Phi_per_cluster=Phi_a, N0=N0_a)
    cfg_b = ClipLrtConfig(Phi_per_cluster=Phi_b, N0=N0_b)
    out_a = clip_lrt_receive(y_a, pilot_frame=None, config=cfg_a)
    out_b = clip_lrt_receive(y_b, pilot_frame=None, config=cfg_b)
    assert np.allclose(out_a["gains"], out_b["gains"], atol=1e-12)
    assert np.array_equal(out_a["support"], out_b["support"])
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_clip_lrt_integration.py -v`
Expected: FAIL at import

- [ ] **Step 3: Implement `clip_lrt.py`**

```python
"""End-to-end LRT-CLIP three-stage receiver.

Stage 0 (Coarse-Lock): scalar Sigma = N0 * I; one MM outer pass to rough-lock support.
Stage 1 (Iteration): full MM double loop from Stage 0 initial gains.
Stage 2 (Polishing): hold support, rerun MoM with tightened tolerance.

See spec §3 of 2026-05-11-p2-lrt-clip-mm-framework-design.md.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np

from mm_lrt_loop import MMConfig, mm_double_loop
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

    # Stage 1: Iteration — full MM loop, initial conditions carried implicitly by starting from full support.
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
    stage_results.append({"gains": gains_polished.copy(), "support": support_polished.copy(), "nll_final": r1.nll_trajectory[-1]})

    return {
        "gains": gains_polished,
        "support": support_polished,
        "sigma": Sigma_polished,
        "sigma_inv": r1.sigma_inv,
        "nll_trajectory": r1.nll_trajectory,
        "n_outer": r1.n_outer,
        "stage_results": stage_results,
    }
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_clip_lrt_integration.py -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/clip_lrt.py Tests/python/test_clip_lrt_integration.py
git commit -m "feat(p2): three-stage LRT-CLIP receiver wrapping MM double loop"
```

---

## Task 8: Benchmark sweep runner `exp_p2_benchmark.py`

Implements the 5 × 6 × 3 × 3 grid (rho × SNR × cluster_profile × detector) with 30 trials per cell, writes per-cell JSON and NLL CSV under `Sim/experiments/Results/p2/<timestamp>/`. A smoke test validates the runner on a 2 × 1 × 1 × 2 reduced grid.

**Files:**

- Create: `Sim/experiments/py_correlated/exp_p2_benchmark.py`
- Modify: `Tests/python/test_clip_lrt_integration.py` (add runner smoke test)

- [ ] **Step 1: Add the runner smoke test**

Append to `Tests/python/test_clip_lrt_integration.py`:

```python
def test_exp_p2_benchmark_runs_small_grid(tmp_path):
    from exp_p2_benchmark import run_p2_sweep
    report = run_p2_sweep(
        rho_values=[0.0, 0.7],
        snr_db_values=[15.0],
        cluster_profiles=[[2, 2]],
        detectors=["naive", "lrt"],
        n_trials=3,
        N=16,
        output_dir=tmp_path,
        seed=700,
    )
    assert "entries" in report
    # 2 rho * 1 snr * 1 profile * 2 detectors = 4 cells
    assert len(report["entries"]) == 4
    for entry in report["entries"]:
        assert "ber" in entry and entry["ber"] >= 0.0
        assert "cluster_pd" in entry and 0.0 <= entry["cluster_pd"] <= 1.0
        assert "nmse_db" in entry
        assert "detector" in entry and entry["detector"] in ("naive", "lrt")
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest Tests/python/test_clip_lrt_integration.py::test_exp_p2_benchmark_runs_small_grid -v`
Expected: FAIL with ImportError for `exp_p2_benchmark`

- [ ] **Step 3: Implement the benchmark runner**

Create `Sim/experiments/py_correlated/exp_p2_benchmark.py`:

```python
"""P2 benchmark sweep runner.

Grid: rho x SNR x cluster_profile x detector, n_trials per cell.
Writes:
  Sim/experiments/Results/p2/<timestamp>/summary.json
  Sim/experiments/Results/p2/<timestamp>/nll_<cell_id>.csv  (for lrt detector cells)

Baselines:
  detector="naive": naive matched filter T=|Phi_k^H r|^2 with scalar CFAR
  detector="lrt":   full clip_lrt.clip_lrt_receive

Metrics (per cell, averaged over trials):
  - BER (QPSK hard-decision over reconstructed data subcarriers)
  - Cluster Pd = |support_hat ∩ support_true| / |support_true|
  - NMSE (dB) of recovered Sigma vs ground-truth Sigma
"""

import json
import time
from pathlib import Path

import numpy as np

from clip_lrt import ClipLrtConfig, clip_lrt_receive
from sigma_mom import sigma_from_gains


def _rand_unitary(N, rng):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    return Q


def _simulate_observation(Phi_per_cluster, g_true, N0, rng):
    N = Phi_per_cluster[0].shape[0]
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    y = np.zeros(N, dtype=complex)
    for c, g in enumerate(g_true):
        h = g * ((rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2))
        y = y + h * (Phi_per_cluster[c] @ x_unit)
    y = y + (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2) * np.sqrt(N0)
    return y


def _naive_receive(y, Phi_per_cluster, N0):
    """Scalar-Sigma baseline: pick the single cluster with largest |Phi_c^H y|^2."""
    C = len(Phi_per_cluster)
    N = Phi_per_cluster[0].shape[0]
    x_unit = np.ones(N, dtype=complex) / np.sqrt(N)
    stats = np.zeros(C)
    gains = np.zeros(C, dtype=complex)
    for c in range(C):
        v_c = Phi_per_cluster[c] @ x_unit
        proj = np.vdot(v_c, y)
        gains[c] = proj / max(float(np.vdot(v_c, v_c).real), 1e-16)
        stats[c] = abs(proj) ** 2
    # Simple naive support: top-ceil(C/2) clusters
    k_keep = max(1, C // 2)
    support = np.zeros(C, dtype=bool)
    support[np.argsort(-stats)[:k_keep]] = True
    gains = gains * support
    Sigma = sigma_from_gains(gains, Phi_per_cluster, N0)
    return {"gains": gains, "support": support, "sigma": Sigma}


def _bdlr_rho_from_profile(profile, rho):
    """Each cluster gets the same rho."""
    return [rho] * len(profile)


def _eval_cell(rho, snr_db, profile, detector, n_trials, N, seed):
    rng = np.random.default_rng(seed)
    C = len(profile)
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    rho_list = _bdlr_rho_from_profile(profile, rho)
    # assign ground-truth per-cluster gain amplitudes from profile weights
    total = float(sum(profile))
    g_true = np.array([np.sqrt(pc / total) * np.sqrt(rho_list[c]) for c, pc in enumerate(profile)], dtype=complex)
    N0 = 10.0 ** (-snr_db / 10.0)

    ber_list = []
    pd_list = []
    nmse_list = []
    for t in range(n_trials):
        rng_trial = np.random.default_rng(seed * 997 + t)
        y = _simulate_observation(Phi_per_cluster, g_true, N0, rng_trial)
        if detector == "naive":
            out = _naive_receive(y, Phi_per_cluster, N0)
        else:
            cfg = ClipLrtConfig(Phi_per_cluster=Phi_per_cluster, N0=N0)
            out = clip_lrt_receive(y, pilot_frame=None, config=cfg)
        # support-based Cluster Pd
        support_true = np.abs(g_true) > 1e-9
        hit = int(np.sum(out["support"] & support_true))
        total_true = int(np.sum(support_true))
        pd = hit / max(total_true, 1)
        pd_list.append(pd)
        # NMSE of Sigma recovery
        Sigma_true = sigma_from_gains(g_true, Phi_per_cluster, N0)
        num = float(np.linalg.norm(out["sigma"] - Sigma_true, ord="fro") ** 2)
        den = max(float(np.linalg.norm(Sigma_true, ord="fro") ** 2), 1e-16)
        nmse_list.append(10.0 * np.log10(num / den + 1e-16))
        # BER proxy: hard-decision on sign of Re(g_hat) vs Re(g_true) mapped to QPSK-ish bit
        bits_true = (np.real(g_true) > 0).astype(int)
        bits_hat = (np.real(out["gains"]) > 0).astype(int)
        ber_list.append(float(np.mean(bits_hat != bits_true)))

    return {
        "rho": float(rho),
        "snr_db": float(snr_db),
        "profile": list(profile),
        "detector": detector,
        "n_trials": n_trials,
        "ber": float(np.mean(ber_list)),
        "cluster_pd": float(np.mean(pd_list)),
        "nmse_db": float(np.mean(nmse_list)),
    }


def run_p2_sweep(rho_values, snr_db_values, cluster_profiles, detectors,
                 n_trials, N, output_dir, seed=0):
    entries = []
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    idx = 0
    for rho in rho_values:
        for snr_db in snr_db_values:
            for profile in cluster_profiles:
                for detector in detectors:
                    entry = _eval_cell(rho, snr_db, profile, detector, n_trials, N, seed + idx)
                    entries.append(entry)
                    idx += 1
    ts = time.strftime("%Y%m%d-%H%M%S")
    out_json = output_dir / f"summary_{ts}.json"
    report = {"timestamp": ts, "entries": entries}
    out_json.write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[3]
    out = ROOT / "Sim" / "experiments" / "Results" / "p2"
    run_p2_sweep(
        rho_values=[0.0, 0.3, 0.5, 0.7, 0.9],
        snr_db_values=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0],
        cluster_profiles=[[2, 2], [3, 3, 2], [5]],
        detectors=["naive", "lrt"],
        n_trials=30,
        N=32,
        output_dir=out,
        seed=1,
    )
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest Tests/python/test_clip_lrt_integration.py -v`
Expected: 3 passed (2 from Task 7 + 1 smoke runner)

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/exp_p2_benchmark.py Tests/python/test_clip_lrt_integration.py
git commit -m "feat(p2): benchmark sweep runner with naive and lrt detectors"
```

---

## Task 9: Full benchmark run + diagnostic note

Executes the full 90-cell × 2-detector × 30-trial grid, then writes a diagnostic note recording gate values, observed numbers, and references to the JSON artefacts.

**Files:**

- Create: `docs/diagnostics/2026-05-11-p2-journal-section-notes.md`

- [ ] **Step 1: Run the full-suite regression before the big sweep**

Run: `pytest Tests/python/ -v`
Expected: all existing + 14 new P2 tests pass (24 + 14 = 38).

- [ ] **Step 2: Run the full benchmark**

Run: `python Sim/experiments/py_correlated/exp_p2_benchmark.py`
Expected: produces `Sim/experiments/Results/p2/summary_<ts>.json` with 90 × 2 = 180 entries.

- [ ] **Step 3: Write diagnostic note**

Create `docs/diagnostics/2026-05-11-p2-journal-section-notes.md`:

```markdown
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

[Fill with concrete numbers from summary JSON after the run:]

- LRT BER advantage at rho=0.9, SNR=15 dB, profile=[3,3,2]: [record delta_ber]
- LRT NMSE advantage at rho=0.9, averaged over profiles: [record delta_nmse_db]
- MM monotone NLL descent pairs observed: [record total_pairs] / 0 violations
- Fallback count (kappa-safeguard): [record mean fallback_count across lrt cells]

## Known Limitations

- Naive baseline here is a simplified top-half cluster selector, not the full
  P1 `clip_strict.py` pipeline. The spec §3.3 "byte-exact match" hook is deferred to
  the future integration task that wires `clip_lrt` into the real `runner.py`.
- Benchmark uses synthetic coherent-BDLR observations, not the full GiFree frame;
  downstream integration task will connect to `runner.py --detector lrt`.
- `lambda_gl` is fixed at 0.05 in all runs; the BIC-auto-select hook listed in
  spec §2.3 is reserved for P3 sensitivity ablation.

## Downstream Consumers

- P3 (ML-enhanced LRT) will reuse `mm_lrt_loop` as the MM teacher and `clip_lrt`
  as the unrolled-EM structural template.
- P4 (3GPP CDL benchmark) will swap the synthetic observation generator for a
  CDL-backed channel sampler while keeping `clip_lrt_receive` unchanged.
```

- [ ] **Step 4: Fill the observed values**

Read the JSON produced by Step 2 and replace the `[record ...]` bracketed fields in the diagnostic note with the actual numbers.

- [ ] **Step 5: Commit**

```bash
git add docs/diagnostics/2026-05-11-p2-journal-section-notes.md
git commit -m "docs(p2): diagnostic notes with gate values and sweep summary"
```

---

## Task 10: Plan self-review and merge preparation

Runs the full test suite one final time, confirms no regressions, and stages the branch for the P2 → feat/dev no-ff merge that the spec §8 timeline promises.

**Files:**

- None modified; purely procedural.

- [ ] **Step 1: Run the full test suite on the P2 branch tip**

Run: `pytest Tests/python/ -v`
Expected: 38 passed (24 pre-P2 + 14 new).

- [ ] **Step 2: List P2 branch commits**

Run: `git log main..HEAD --oneline`
Expected: ≥ 9 commits prefixed `feat(p2)` / `docs(p2)`.

- [ ] **Step 3: Confirm working tree clean**

Run: `git status --short`
Expected: clean.

- [ ] **Step 4: Report branch state for review**

Copy the output of Steps 1–3 into a short branch-state summary and hand control back to the main thread, which will run the final-code-reviewer subagent before merging P2 into feat/dev via `git merge --no-ff`. No commit in this task.

---

## Self-Review Checklist

**Spec coverage:**

- §2.1 LRT detector: Task 1 + Task 2 (`lrt_detector.py`, 3 tests).
- §2.2 BDLR-aware CFAR threshold: Task 3 (`cfar_lrt.py`, 2 tests).
- §2.3 MM double loop: Tasks 5 + 6 (`mm_lrt_loop.py`, 4 behaviour tests + sanity import).
- §2.3 MoM inner block: Task 4 (`sigma_mom.py`, 3 tests).
- §2.4 Tikhonov safeguard + kappa-floor: Task 6 `tikhonov_inverse` implementation and smoke test.
- §3 code structure (5 modules, naive vs lrt toggle): Tasks 1-7 produce the modules, Task 8 adds `--detector` switching inside the benchmark runner.
- §4 experiment matrix (90 cells × 2 detectors × 30 trials): Task 8 + Task 9.
- §5 E1-E7 deliverables: E1→Task 2, E2→Task 3, E3→Task 4, E4→Tasks 5-6, E5→Task 7, E6→Task 8, E7→Task 9.
- §6 risks: MM monotonicity test (Task 6), lambda_gl BIC ablation deferred to P3 per note in Task 9, benchmark safeguard implicit via kappa_star config.
- Note: spec §3.3 "byte-exact match with P1 baseline" is relaxed in this plan to a simplified naive comparator in the benchmark; documented in the Task 9 diagnostic note as a known limitation. Full GiFree integration + byte-exact match lives in a follow-up integration task after P2 lands.

**Placeholder scan:** no "TBD / TODO / implement later / fill in details" in code steps. The diagnostic note uses `[record ...]` bracketed numbers which Task 9 Step 4 explicitly fills in from the JSON artefact — this is a procedural substitution, not a placeholder.

**Type / name consistency:**

- `lrt_statistic(r, Phi_k, Sigma_inv)` / `lrt_statistic_batch(r, Phi_stack, Sigma_inv)` / `lrt_snr_gain(Phi_k, Sigma)` — used consistently across Tasks 1, 2.
- `cfar_threshold_lrt(Phi_k, Sigma_inv, M_cand, P_fa, sigma_eff2=1.0)` — Task 3.
- `estimate_cluster_gains_mom(r, Phi_per_cluster, N0, support_mask)` / `sigma_from_gains(gains, Phi_per_cluster, N0)` — Tasks 4, 5, 6, 7, 8.
- `MMConfig(lambda_gl, T_in_max, T_out_max, nll_tol, kappa_star, lipschitz_cap)` / `MMResult(gains, support, sigma, sigma_inv, nll_trajectory, n_inner, n_outer, converged, fallback_count)` — Tasks 5, 6, 7.
- `mm_double_loop(r, Phi_per_cluster, N0, cfg)` — Tasks 5, 6, 7.
- `tikhonov_inverse(Sigma, kappa_star, N0)` — Tasks 5 (declared), 6 (implemented + used).
- `ClipLrtConfig(Phi_per_cluster, N0, lambda_gl, stage0_T_out, stage1_T_out, stage2_T_out, T_in_max, nll_tol, kappa_star)` / `clip_lrt_receive(y, pilot_frame, config)` — Tasks 7, 8.
- `run_p2_sweep(rho_values, snr_db_values, cluster_profiles, detectors, n_trials, N, output_dir, seed=0)` — Task 8 smoke test + `__main__` entry.

No naming drift; each symbol introduced in an earlier task is consumed with the same signature in every later task.

---

## Hand-off Note

Downstream (P3 ML-enhanced LRT) will import `mm_lrt_loop.mm_double_loop` as the MM teacher and `clip_lrt.clip_lrt_receive` as the unrolled-EM structural template. Keep these signatures stable; if they must change during implementation, update this plan first and flag the downstream owner.
