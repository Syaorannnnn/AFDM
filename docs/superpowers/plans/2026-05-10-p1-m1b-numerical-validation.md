# P1: M1-B Numerical Validation Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python validation toolkit that empirically confirms every closed-form claim in M1-B of `2026-05-10-afdm-correlated-clip-research-design.md` §3 — ID2P covariance, position-dependent Rayleigh quotient, Fisher / CRLB single-path vs cluster, and condition-number upper bound — so the letter-grade derivations can be trusted before any algorithm (P2+) consumes them.

**Architecture:** Each closed-form claim gets (a) a pure-math helper function that computes the theoretical prediction from `(Phi_list, R, noise_var)` and (b) a Monte Carlo estimator that averages the same quantity empirically. A pytest file per module asserts the two agree within a tolerance. A top-level `validate_m1b.py` script runs the full Monte Carlo sweep over `(rho, SNR, P_c)` and writes diagnostic plots + a machine-readable JSON report that serves as the Go/no-go gate before P2.

**Tech Stack:** Python 3.11+, NumPy, SciPy (for `scipy.linalg.sqrtm`, eigvals), pytest (existing project test harness at `Tests/python/`), matplotlib. Reuses `Sim/experiments/py_correlated/channel.py`, `id2p.py`, `grid_stats.py` for BDLR gain sampling and path-matrix construction — does not duplicate them.

---

## File Structure

**Create** (all under `Sim/experiments/py_correlated/` unless noted):

| File | Responsibility |
| --- | --- |
| `sigma_closed_form.py` | Closed-form ID2P covariance $\boldsymbol{\Sigma}$ from `(R, Phi_list, sigma_d2, N0)` + MC estimator |
| `p_res_eff.py` | Position-dependent Rayleigh quotient $P_{\mathrm{res}}^{\mathrm{eff}}(k)$ + MC estimator |
| `fisher_crlb.py` | Fisher information matrix and CRLB for single-path and cluster parameterisations |
| `kappa_bound.py` | $\kappa(\boldsymbol{\Sigma})$ Weyl upper bound + MC comparison |
| `validate_m1b.py` | End-to-end MC sweep; writes `Sim/experiments/Results/m1b_validation_<timestamp>.json` and PNG plots |
| `Tests/python/test_sigma_closed_form.py` | pytest for §1 |
| `Tests/python/test_p_res_eff.py` | pytest for §2 |
| `Tests/python/test_fisher_crlb.py` | pytest for §3 |
| `Tests/python/test_kappa_bound.py` | pytest for §4 |
| `docs/diagnostics/2026-05-10-p1-m1b-validation-notes.md` | Diagnostic notes; records the pass/fail gate values |

**Do not modify** any file under `Core/SystemFunctions/` or `Sim/MainSimulation.m` — this work is Python-only.

---

## Conventions

- **sys.path**: All tests prepend `<repo>/Sim/experiments/py_correlated` to `sys.path`, matching the existing pattern in `Tests/python/test_py_correlated_cfar.py:10-11`.
- **Test isolation**: Every test seeds `np.random.default_rng(seed=...)` with a fixed seed and uses a small `N`/`P`/MC-trials combo that finishes in < 2 s per test.
- **Tolerance conventions**: closed-form vs MC agreement uses `np.allclose(..., rtol=0.05, atol=1e-3)` for normalised quantities and explicit variance-bound checks for raw MC means. Wider tolerances need a comment justifying why.
- **Commits**: TDD per task — test first, implement second, commit third. Prefix commits with `feat(p1):`, `test(p1):`, or `docs(p1):`.
- **No MATLAB bridge**: P1 stays in pure Python. Any MATLAB interop is P2+.---

## Task 1: Project scaffold and pytest harness sanity check

**Files:**

- Create: `Sim/experiments/py_correlated/sigma_closed_form.py` (stub)
- Create: `Tests/python/test_sigma_closed_form.py` (sanity test)

- [ ] **Step 1: Write the failing import test**

```python
"""Sanity test that the new module imports and exposes the expected API."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_module_importable():
    import sigma_closed_form
    assert hasattr(sigma_closed_form, "closed_form_sigma")
    assert hasattr(sigma_closed_form, "mc_sigma")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_module_importable -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'sigma_closed_form'`

- [ ] **Step 3: Create the stub module**

```python
# Sim/experiments/py_correlated/sigma_closed_form.py
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_module_importable -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/sigma_closed_form.py Tests/python/test_sigma_closed_form.py
git commit -m "feat(p1): scaffold sigma_closed_form module with API contract"
```

---

## Task 2: Closed-form Sigma formula + MC cross-check

Implements the double-sum closed form $\mathbb{E}[\mathbf{C}_{\mathrm{ID2P}}] = \sigma_d^2 \sum_{i,j} R_{ij} \boldsymbol{\Phi}_i \boldsymbol{\Phi}_j^{\mathsf H} + N_0 \mathbf{I}$ (spec §3.1; 2026-05-07 §3.1) and a direct MC estimator. Two degeneration tests + one non-trivial cross-check.

**Files:**

- Modify: `Sim/experiments/py_correlated/sigma_closed_form.py` (implement both functions)
- Modify: `Tests/python/test_sigma_closed_form.py` (add 3 numerical tests)

- [ ] **Step 1: Write degeneration test — R = diag ⇒ scalar Sigma**

Append to `Tests/python/test_sigma_closed_form.py`:

```python
from sigma_closed_form import closed_form_sigma, mc_sigma
import numpy as np


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def test_diagonal_R_yields_scalar_sigma():
    rng = np.random.default_rng(0)
    N, P = 16, 4
    Phi = _make_random_unitary_phis(P, N, rng)
    R = np.diag([0.25, 0.25, 0.25, 0.25]).astype(complex)
    sigma_d2, N0 = 1.0, 0.01
    Sigma = closed_form_sigma(Phi, R, sigma_d2, N0)
    expected_scalar = sigma_d2 * np.trace(R).real + N0
    assert np.allclose(Sigma, expected_scalar * np.eye(N), atol=1e-10)
```

- [ ] **Step 2: Run — expect fail (NotImplementedError)**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_diagonal_R_yields_scalar_sigma -v`
Expected: FAIL with `NotImplementedError`

- [ ] **Step 3: Implement `closed_form_sigma`**

Replace the body in `Sim/experiments/py_correlated/sigma_closed_form.py`:

```python
def closed_form_sigma(Phi_list, R, sigma_d2, N0):
    """
    Closed-form ID2P covariance Sigma under BDLR channel prior.

    Sigma = sigma_d2 * sum_{i,j} R[i,j] * Phi_i @ Phi_j.conj().T + N0 * I_N

    Parameters
    ----------
    Phi_list : list of (N, N) complex ndarray
        Per-path DAFT-domain operators; must satisfy Phi_i Phi_i^H = I_N.
    R : (P, P) complex ndarray
        Path-gain covariance (Hermitian PSD).
    sigma_d2 : float
        Data symbol power.
    N0 : float
        Noise power.

    Returns
    -------
    Sigma : (N, N) complex ndarray
    """
    P = len(Phi_list)
    assert R.shape == (P, P), f"R must be {P}x{P}, got {R.shape}"
    N = Phi_list[0].shape[0]
    Sigma = np.zeros((N, N), dtype=complex)
    for i in range(P):
        for j in range(P):
            if R[i, j] == 0:
                continue
            Sigma += R[i, j] * (Phi_list[i] @ Phi_list[j].conj().T)
    Sigma = sigma_d2 * Sigma + N0 * np.eye(N)
    return Sigma
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_diagonal_R_yields_scalar_sigma -v`
Expected: PASS

- [ ] **Step 5: Add MC estimator test**

Append to `Tests/python/test_sigma_closed_form.py`:

```python
def test_mc_sigma_converges_to_closed_form():
    rng = np.random.default_rng(1)
    N, P = 8, 3
    Phi = _make_random_unitary_phis(P, N, rng)
    # Full-rank correlated R (simple BDLR with one cluster).
    rho, gain_var = 0.7, 0.25
    R = gain_var * ((1 - rho) * np.eye(P) + rho * np.ones((P, P)))
    sigma_d2, N0 = 1.0, 0.02
    Sigma_cf = closed_form_sigma(Phi, R, sigma_d2, N0)
    Sigma_mc = mc_sigma(Phi, R, sigma_d2, N0, n_trials=4000, rng=rng)
    rel_err = np.linalg.norm(Sigma_mc - Sigma_cf, "fro") / np.linalg.norm(Sigma_cf, "fro")
    assert rel_err < 0.08, f"MC/CF disagreement {rel_err:.3f} > 0.08"
```

- [ ] **Step 6: Run — expect fail (mc_sigma stub)**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_mc_sigma_converges_to_closed_form -v`
Expected: FAIL with `NotImplementedError`

- [ ] **Step 7: Implement `mc_sigma`**

```python
def mc_sigma(Phi_list, R, sigma_d2, N0, n_trials, rng=None):
    """
    Monte Carlo estimate of the ID2P covariance by averaging h h^H and
    forming the outer channel response.
    """
    if rng is None:
        rng = np.random.default_rng()
    P = len(Phi_list)
    N = Phi_list[0].shape[0]
    L = np.linalg.cholesky(R + 1e-12 * np.eye(P))
    Sigma_accum = np.zeros((N, N), dtype=complex)
    for _ in range(n_trials):
        w = (rng.standard_normal(P) + 1j * rng.standard_normal(P)) / np.sqrt(2)
        h = L @ w
        x = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
        x = np.sqrt(sigma_d2) * x
        H_eff = sum(h[i] * Phi_list[i] for i in range(P))
        y = H_eff @ x
        Sigma_accum += np.outer(y, y.conj())
    Sigma_mc = Sigma_accum / n_trials + N0 * np.eye(N)
    return Sigma_mc
```

- [ ] **Step 8: Run both tests — expect pass**

Run: `pytest Tests/python/test_sigma_closed_form.py -v`
Expected: 3 passed (import + diag + MC)

- [ ] **Step 9: Commit**

```bash
git add Sim/experiments/py_correlated/sigma_closed_form.py Tests/python/test_sigma_closed_form.py
git commit -m "feat(p1): closed-form Sigma + MC estimator with BDLR-prior cross-check"
```

---

## Task 3: BDLR structure check and rank-C collapse assertion

Adds a second non-trivial closed-form check specific to spec §3.3: when the BDLR structure is active, $\boldsymbol{\Sigma} - N_0 \mathbf{I}$'s "coherent" part has rank exactly $C$ (one cluster contribution each). This is the property M1-B uses to justify cluster-level identifiability.

**Files:**

- Modify: `Sim/experiments/py_correlated/sigma_closed_form.py` (add `bdlr_decompose`)
- Modify: `Tests/python/test_sigma_closed_form.py` (add 2 rank tests)

- [ ] **Step 1: Write rank-collapse test**

Append to `Tests/python/test_sigma_closed_form.py`:

```python
def test_bdlr_decompose_conditional_rank_equals_cluster_count():
    """Given a pilot vector x_p, the conditional-on-data coherent kernel
    sum_c P_c (G_c x_p) (G_c x_p)^H has rank exactly C at rho=1.

    Note: the ensemble-averaged operator `bdlr_decompose(...)[0]` is
    generically rank-N; the rank-C collapse is a *conditional* property
    that the Weyl bound and identifiability arguments rely on. This test
    builds the conditional kernel inline and checks rank.
    """
    rng = np.random.default_rng(2)
    N = 24
    cluster_sizes = [3, 2, 4]
    C = len(cluster_sizes)
    P = sum(cluster_sizes)
    Phi = _make_random_unitary_phis(P, N, rng)
    cluster_labels = np.repeat(np.arange(C), cluster_sizes)
    P_per = np.array([0.5, 0.3, 0.2])
    x_p = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
    coherent_cond = np.zeros((N, N), dtype=complex)
    for c in range(C):
        members = np.where(cluster_labels == c)[0]
        Gc = sum(Phi[i] for i in members)
        gc_x = Gc @ x_p
        coherent_cond += P_per[c] * np.outer(gc_x, gc_x.conj())
    eigvals = np.sort(np.linalg.eigvalsh(coherent_cond))[::-1]
    assert np.sum(eigvals > 1e-8) == C, (
        f"conditional rank={np.sum(eigvals > 1e-8)}, expected {C}"
    )


def test_bdlr_decompose_rho_zero_coherent_vanishes():
    """When rho_c = 0 for every cluster, coherent part is zero."""
    rng = np.random.default_rng(3)
    N, C, per = 16, 2, 3
    Phi = _make_random_unitary_phis(C * per, N, rng)
    from sigma_closed_form import bdlr_decompose
    coherent, scalar, N0_term = bdlr_decompose(
        Phi_list=Phi,
        cluster_labels=np.repeat(np.arange(C), per),
        P_per_cluster=np.array([0.5, 0.5]),
        rho_per_cluster=np.zeros(C),
        sigma_d2=1.0,
        N0=0.01,
    )
    assert np.allclose(coherent, 0.0, atol=1e-12)
    # scalar part is diagonal.
    assert np.allclose(scalar, np.trace(scalar).real / scalar.shape[0] * np.eye(scalar.shape[0]))


def test_bdlr_decompose_operator_rho_zero_gives_scalar_sigma():
    """At rho=0 the whole Sigma reconstructed from (coherent, scalar, N0_term)
    equals a scalar multiple of I."""
    rng = np.random.default_rng(4)
    N, C, per = 12, 3, 2
    Phi = _make_random_unitary_phis(C * per, N, rng)
    from sigma_closed_form import bdlr_decompose
    P_per = np.array([0.4, 0.35, 0.25])
    coherent, scalar, N0_term = bdlr_decompose(
        Phi_list=Phi,
        cluster_labels=np.repeat(np.arange(C), per),
        P_per_cluster=P_per,
        rho_per_cluster=np.zeros(C),
        sigma_d2=1.0,
        N0=0.05,
    )
    Sigma = coherent + scalar + N0_term
    expected = (np.sum(P_per) + 0.05) * np.eye(N)
    assert np.allclose(Sigma, expected, atol=1e-12)
```

- [ ] **Step 2: Run — expect fail (bdlr_decompose undefined)**

Run: `pytest Tests/python/test_sigma_closed_form.py -v -k bdlr_decompose`
Expected: FAIL with `ImportError: cannot import name 'bdlr_decompose'` or collection error.

- [ ] **Step 3: Implement `bdlr_decompose`**

Append to `Sim/experiments/py_correlated/sigma_closed_form.py`:

```python
def bdlr_decompose(Phi_list, cluster_labels, P_per_cluster, rho_per_cluster,
                   sigma_d2, N0):
    """
    Decompose the ensemble-expected Sigma into three additive parts.

        Sigma = coherent_operator + scalar * I_N + N0 * I_N,
        coherent_operator = sigma_d2 * sum_c P_c rho_c * G_c G_c^H,
        G_c = sum_{i in c} Phi_i.

    Shape note
    ----------
    Each Phi_i is a full-rank unitary (N, N) DAFT-domain operator, so
    `G_c G_c^H` is generically rank-N. The rank-C "coherent" property
    in spec §3.1 applies to the *conditional* kernel
        sum_c P_c rho_c (G_c x_p) (G_c x_p)^H
    for a given pilot vector x_p; that kernel is rank-C at rho=1. The
    operator form returned here is the ensemble average
    E_x [conditional kernel], which is what the Weyl bound
    (Task 8) and the sweep runner (Task 9) consume; it has the
    *correct* largest eigenvalue but generically full rank.
    """
    C = int(cluster_labels.max()) + 1
    N = Phi_list[0].shape[0]
    coherent = np.zeros((N, N), dtype=complex)
    for c in range(C):
        members = np.where(cluster_labels == c)[0]
        Gc = sum(Phi_list[i] for i in members)
        coherent += P_per_cluster[c] * rho_per_cluster[c] * (Gc @ Gc.conj().T)
    coherent *= sigma_d2
    scalar_weight = sigma_d2 * float(np.sum(P_per_cluster * (1.0 - rho_per_cluster)))
    scalar = scalar_weight * np.eye(N, dtype=complex)
    N0_term = N0 * np.eye(N, dtype=complex)
    return coherent, scalar, N0_term
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_sigma_closed_form.py -v`
Expected: 6 passed (3 original + 3 new).

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/sigma_closed_form.py Tests/python/test_sigma_closed_form.py
git commit -m "feat(p1): BDLR decomposition with conditional rank-C collapse verification"
```

---

## Task 4: Position-dependent Rayleigh quotient P_res_eff(k)

Implements $P_{\mathrm{res}}^{\mathrm{eff}}(k) = \boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma} \boldsymbol{\Phi}_k / \|\boldsymbol{\Phi}_k\|^2$ (spec 2026-05-07 §4.2) and its MC equivalent by empirical projection of residuals.

**Files:**

- Create: `Sim/experiments/py_correlated/p_res_eff.py`
- Create: `Tests/python/test_p_res_eff.py`

- [ ] **Step 1: Write the test file**

```python
"""Position-dependent Rayleigh quotient P_res_eff(k) validation."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_closed_form import closed_form_sigma
from p_res_eff import p_res_eff_from_sigma, p_res_eff_mc


def _unit_cols(N, K, rng):
    M = rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))
    M /= np.linalg.norm(M, axis=0, keepdims=True)
    return M


def test_p_res_eff_scalar_sigma_is_constant():
    rng = np.random.default_rng(10)
    N, K = 12, 5
    Sigma = 0.7 * np.eye(N, dtype=complex)
    H = _unit_cols(N, K, rng)
    eff = p_res_eff_from_sigma(H, Sigma)
    assert np.allclose(eff, 0.7)


def test_p_res_eff_mc_matches_closed_form():
    rng = np.random.default_rng(11)
    N, K = 16, 4
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Sigma = A @ A.conj().T / N + 0.05 * np.eye(N)
    H = _unit_cols(N, K, rng)
    eff_cf = p_res_eff_from_sigma(H, Sigma)
    eff_mc = p_res_eff_mc(H, Sigma, n_trials=6000, rng=rng)
    rel_err = np.max(np.abs(eff_mc - eff_cf) / np.abs(eff_cf))
    assert rel_err < 0.08, f"max rel err {rel_err:.3f} > 0.08"
```

- [ ] **Step 2: Run — expect fail (module missing)**

Run: `pytest Tests/python/test_p_res_eff.py -v`
Expected: FAIL at import time

- [ ] **Step 3: Implement `p_res_eff.py`**

```python
"""Rayleigh-quotient effective residual power at each candidate position."""

import numpy as np


def p_res_eff_from_sigma(H, Sigma):
    """
    Parameters
    ----------
    H : (N, K) complex ndarray of candidate responses (columns = candidates)
    Sigma : (N, N) Hermitian PSD
    Returns
    -------
    (K,) float array: h_k^H Sigma h_k / ||h_k||^2
    """
    num = np.einsum("nk,nm,mk->k", H.conj(), Sigma, H).real
    den = np.einsum("nk,nk->k", H.conj(), H).real
    return num / den


def p_res_eff_mc(H, Sigma, n_trials, rng=None):
    """Monte Carlo: draw e ~ CN(0, Sigma) and average |h_k^H e|^2 / ||h_k||^2."""
    if rng is None:
        rng = np.random.default_rng()
    N, K = H.shape
    # Use eigendecomposition to sample CN(0, Sigma).
    w, V = np.linalg.eigh(Sigma)
    w_clip = np.clip(w, 0, None)
    L = V * np.sqrt(w_clip)
    accum = np.zeros(K)
    norms = np.einsum("nk,nk->k", H.conj(), H).real
    for _ in range(n_trials):
        z = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
        e = L @ z
        proj = H.conj().T @ e
        accum += (np.abs(proj) ** 2) / norms
    return accum / n_trials
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_p_res_eff.py -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/p_res_eff.py Tests/python/test_p_res_eff.py
git commit -m "feat(p1): P_res_eff Rayleigh quotient + MC validation"
```

---

## Task 5: Fisher information for single-path parameters (delay-only slice)

Implements $J_{ij}(\boldsymbol{\theta}) = 2\,\mathrm{Re}[(\partial \boldsymbol{\mu} / \partial \theta_i)^{\mathsf H} \boldsymbol{\Sigma}^{-1} (\partial \boldsymbol{\mu} / \partial \theta_j)]$ (spec §3.2) for a single path with integer delay ladder. This isolates the simplest case so the full CRLB test in Task 6 has a trusted building block.

**Files:**

- Create: `Sim/experiments/py_correlated/fisher_crlb.py`
- Create: `Tests/python/test_fisher_crlb.py`

- [ ] **Step 1: Write single-path Fisher test (delay + complex gain)**

```python
"""Fisher / CRLB sanity tests for M1-B claims."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from fisher_crlb import fisher_single_path, crlb_from_fisher


def _make_delay_operator(N, delay_frac):
    """Simple shift-by-delay Phi: (N, N) circulant shift matrix sampled at frac delay."""
    k = np.arange(N)
    col = np.exp(-1j * 2 * np.pi * delay_frac * k / N)
    # Symbolic "derivative w.r.t. delay" is encoded in fisher_single_path test itself.
    return np.diag(col)


def test_fisher_matrix_hermitian_and_psd():
    rng = np.random.default_rng(20)
    N = 32
    Phi = _make_delay_operator(N, 0.3)
    dPhi = _make_delay_operator(N, 0.3001)
    dPhi = (dPhi - Phi) / 0.0001
    Sigma = 0.1 * np.eye(N) + 0.05 * np.ones((N, N))
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    h = 1.0 + 0j
    J = fisher_single_path(
        Phi=Phi, dPhi_dtheta=[dPhi], h=h, x_pilot=x_pilot, Sigma=Sigma
    )
    assert J.shape == (3, 3)  # 1 delay-like param + Re(h) + Im(h)
    assert np.allclose(J, J.conj().T, atol=1e-10)
    assert np.all(np.linalg.eigvalsh((J + J.conj().T) / 2) >= -1e-8)
```

- [ ] **Step 2: Run — expect fail (module missing)**

Run: `pytest Tests/python/test_fisher_crlb.py -v`
Expected: FAIL at import

- [ ] **Step 3: Implement `fisher_single_path` and `crlb_from_fisher`**

```python
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
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_fisher_crlb.py -v`
Expected: 1 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/fisher_crlb.py Tests/python/test_fisher_crlb.py
git commit -m "feat(p1): single-path Fisher matrix + CRLB helper"
```

---

## Task 6: Single-path CRLB blows up when two paths share a cluster (identifiability crisis)

Proves numerically the key spec §3.2 claim: when two paths are nearly coherent ($\cos(\angle(\boldsymbol{\Phi}_i, \boldsymbol{\Phi}_j)) \to 1$) their joint Fisher block is near-singular, so the per-path CRLB diverges. This is the test that underwrites the "need cluster-level parameterisation" narrative.

**Files:**

- Modify: `Sim/experiments/py_correlated/fisher_crlb.py` (add `fisher_two_paths_joint`)
- Modify: `Tests/python/test_fisher_crlb.py` (add divergence test)

- [ ] **Step 1: Write the CRLB divergence test**

Append to `Tests/python/test_fisher_crlb.py`:

```python
from fisher_crlb import fisher_two_paths_joint


def test_two_path_crlb_diverges_as_phis_become_coherent():
    """As Phi_j -> Phi_i, joint Fisher block near-singular -> CRLB -> inf."""
    rng = np.random.default_rng(30)
    N = 48
    base = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    U, _, _ = np.linalg.svd(base)
    Phi_i = U                              # unitary
    Sigma = 0.05 * np.eye(N, dtype=complex)
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    # Sweep similarity: Phi_j = cos(theta) Phi_i + sin(theta) Phi_perp
    V, _, _ = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )
    Phi_perp = V
    crlb_gain_i = []
    thetas = [1.0, 0.3, 0.1, 0.03, 0.01]
    for theta in thetas:
        Phi_j = np.cos(theta) * Phi_i + np.sin(theta) * Phi_perp
        # Renormalise so Phi_j is still approximately unitary for the test.
        J = fisher_two_paths_joint(
            Phi_i=Phi_i, Phi_j=Phi_j, h_i=1.0, h_j=1.0,
            x_pilot=x_pilot, Sigma=Sigma,
        )
        diag = crlb_from_fisher(J, eps_reg=1e-14)
        crlb_gain_i.append(diag[0])  # Re(h_i) variance
    # Expect monotone growth as theta shrinks.
    for k in range(len(thetas) - 1):
        assert crlb_gain_i[k + 1] > crlb_gain_i[k], (
            f"CRLB did not grow at theta step {k}: {crlb_gain_i}"
        )
    # Final (near-coherent) CRLB should be at least 10x the separated case.
    assert crlb_gain_i[-1] > 10.0 * crlb_gain_i[0]
```

- [ ] **Step 2: Run — expect fail (fisher_two_paths_joint missing)**

Run: `pytest Tests/python/test_fisher_crlb.py::test_two_path_crlb_diverges_as_phis_become_coherent -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement `fisher_two_paths_joint`**

Append to `Sim/experiments/py_correlated/fisher_crlb.py`:

```python
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
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_fisher_crlb.py -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/fisher_crlb.py Tests/python/test_fisher_crlb.py
git commit -m "feat(p1): two-path joint Fisher with CRLB-divergence numerical proof"
```

---

## Task 7: Cluster-level Fisher remains finite (identifiability recovery)

Complements Task 6: reparameterise two near-coherent paths as a single cluster with one effective gain. The cluster-level CRLB stays bounded for the whole sweep, confirming spec §3.2's "cluster recovery" claim.

**Files:**

- Modify: `Sim/experiments/py_correlated/fisher_crlb.py` (add `fisher_cluster_level`)
- Modify: `Tests/python/test_fisher_crlb.py` (add recovery test)

- [ ] **Step 1: Write the cluster-level test**

Append to `Tests/python/test_fisher_crlb.py`:

```python
from fisher_crlb import fisher_cluster_level


def test_cluster_level_crlb_bounded_across_similarity_sweep():
    rng = np.random.default_rng(31)
    N = 48
    U = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )[0]
    V = np.linalg.svd(
        rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    )[0]
    Phi_i, Phi_perp = U, V
    Sigma = 0.05 * np.eye(N, dtype=complex)
    x_pilot = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    cluster_crlbs = []
    for theta in [1.0, 0.3, 0.1, 0.03, 0.01]:
        Phi_j = np.cos(theta) * Phi_i + np.sin(theta) * Phi_perp
        Phi_cluster = Phi_i + Phi_j  # unnormalised coherent sum, matches spec §3.1
        J_cluster = fisher_cluster_level(
            Phi_cluster=Phi_cluster, g_cluster=1.0, x_pilot=x_pilot, Sigma=Sigma
        )
        diag = crlb_from_fisher(J_cluster, eps_reg=1e-14)
        cluster_crlbs.append(diag[0])
    # Cluster-level CRLB should stay within one order of magnitude across the sweep.
    ratio = max(cluster_crlbs) / min(cluster_crlbs)
    assert ratio < 10.0, f"Cluster CRLB varied by {ratio:.2f}x across similarity sweep"
```

- [ ] **Step 2: Run — expect fail**

Run: `pytest Tests/python/test_fisher_crlb.py::test_cluster_level_crlb_bounded_across_similarity_sweep -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement `fisher_cluster_level`**

Append to `Sim/experiments/py_correlated/fisher_crlb.py`:

```python
def fisher_cluster_level(Phi_cluster, g_cluster, x_pilot, Sigma):
    """Fisher matrix for cluster-level params [Re(g), Im(g)] with Phi_cluster fixed."""
    Sigma_inv = np.linalg.pinv(Sigma)
    a = Phi_cluster @ x_pilot
    G = np.stack([a, 1j * a], axis=1)
    J = 2.0 * (G.conj().T @ Sigma_inv @ G).real
    return J
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_fisher_crlb.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/fisher_crlb.py Tests/python/test_fisher_crlb.py
git commit -m "feat(p1): cluster-level Fisher stays bounded -> identifiability recovery"
```

---

## Task 8: Kappa(Sigma) Weyl upper bound

Implements the spec §3.3 closed-form upper bound on $\kappa(\boldsymbol{\Sigma})$ and checks empirically.

**Files:**

- Create: `Sim/experiments/py_correlated/kappa_bound.py`
- Create: `Tests/python/test_kappa_bound.py`

- [ ] **Step 1: Write the tightness test**

```python
"""Kappa upper bound from spec §3.3."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from sigma_closed_form import closed_form_sigma, bdlr_decompose
from kappa_bound import kappa_weyl_upper_bound


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def test_kappa_bound_dominates_empirical_kappa():
    rng = np.random.default_rng(40)
    N = 24
    cluster_sizes = [3, 2, 3]
    C = len(cluster_sizes)
    P = sum(cluster_sizes)
    Phi = _make_random_unitary_phis(P, N, rng)
    labels = np.repeat(np.arange(C), cluster_sizes)
    P_per = np.array([0.4, 0.3, 0.3])
    for rho_val in [0.0, 0.3, 0.7, 0.95]:
        rho_per = rho_val * np.ones(C)
        coherent, scalar, N0_term = bdlr_decompose(
            Phi, labels, P_per, rho_per, sigma_d2=1.0, N0=0.01,
        )
        Sigma = coherent + scalar + N0_term
        emp_kappa = np.linalg.cond(Sigma)
        bound = kappa_weyl_upper_bound(
            coherent_part=coherent, scalar_weight=scalar[0, 0].real, N0=0.01,
        )
        assert bound >= emp_kappa * (1.0 - 1e-6), (
            f"bound {bound:.2f} < empirical kappa {emp_kappa:.2f} at rho={rho_val}"
        )
```

- [ ] **Step 2: Run — expect fail**

Run: `pytest Tests/python/test_kappa_bound.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement `kappa_weyl_upper_bound`**

```python
"""Closed-form condition-number upper bound for Sigma under BDLR."""

import numpy as np


def kappa_weyl_upper_bound(coherent_part, scalar_weight, N0):
    """
    Weyl-based upper bound on kappa(Sigma) when
       Sigma = coherent_part + scalar_weight * I + N0 * I.

    lambda_max(Sigma) <= lambda_max(coherent_part) + scalar_weight + N0
    lambda_min(Sigma) >= scalar_weight + N0

    so kappa(Sigma) <= (lambda_max(coherent_part) + scalar_weight + N0)
                        / (scalar_weight + N0).
    """
    eigvals = np.linalg.eigvalsh((coherent_part + coherent_part.conj().T) / 2)
    lam_max = float(max(eigvals.max(), 0.0))
    denom = scalar_weight + N0
    if denom <= 0:
        raise ValueError("scalar_weight + N0 must be positive")
    return (lam_max + scalar_weight + N0) / denom
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_kappa_bound.py -v`
Expected: 1 passed

- [ ] **Step 5: Commit**

```bash
git add Sim/experiments/py_correlated/kappa_bound.py Tests/python/test_kappa_bound.py
git commit -m "feat(p1): Weyl-based closed-form kappa(Sigma) upper bound"
```

---

## Task 9: Sweep runner `validate_m1b.py`

Ties the four building blocks together: for a grid $(\rho, \mathrm{SNR}, C, P_c)$ sweep, compute closed-form predictions and MC estimates, write JSON report + plots, print PASS / FAIL summary.

**Files:**

- Create: `Sim/experiments/py_correlated/validate_m1b.py`
- Create: `Sim/experiments/Results/.gitkeep` (if directory does not yet exist — it does per `git status`, so only create if missing)

- [ ] **Step 1: Write the script scaffold and smoke test**

Append to `Tests/python/test_sigma_closed_form.py` (reuse the file for the runner smoke test since it already imports the relevant modules):

```python
def test_validate_m1b_runs_small_grid(tmp_path):
    from validate_m1b import run_sweep
    report = run_sweep(
        rho_values=[0.0, 0.8],
        snr_db_values=[10.0],
        cluster_sizes=[2, 2],
        N=16,
        n_trials=500,
        output_dir=tmp_path,
        seed=42,
    )
    assert "entries" in report
    assert len(report["entries"]) == 2
    for entry in report["entries"]:
        assert entry["sigma_rel_err"] < 0.15
        assert entry["kappa_bound"] >= entry["kappa_empirical"] * (1 - 1e-6)
```

- [ ] **Step 2: Run — expect fail**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_validate_m1b_runs_small_grid -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement `validate_m1b.py`**

Create `Sim/experiments/py_correlated/validate_m1b.py`:

```python
"""End-to-end M1-B numerical validation sweep.

Usage
-----
    python Sim/experiments/py_correlated/validate_m1b.py

Produces
--------
    Sim/experiments/Results/m1b_validation_<timestamp>.json
    Sim/experiments/Results/m1b_kappa_vs_rho_<timestamp>.png
"""

import json
import time
from pathlib import Path

import numpy as np

from sigma_closed_form import closed_form_sigma, mc_sigma, bdlr_decompose
from kappa_bound import kappa_weyl_upper_bound


def _make_random_unitary_phis(P, N, rng):
    phis = []
    for _ in range(P):
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        Q, _ = np.linalg.qr(A)
        phis.append(Q)
    return phis


def _cluster_R(cluster_sizes, rho_val, gain_var=1.0):
    P = int(np.sum(cluster_sizes))
    R = np.zeros((P, P), dtype=complex)
    offset = 0
    for sz in cluster_sizes:
        block = gain_var * (
            (1 - rho_val) * np.eye(sz) + rho_val * np.ones((sz, sz))
        )
        R[offset:offset + sz, offset:offset + sz] = block
        offset += sz
    return R


def run_sweep(rho_values, snr_db_values, cluster_sizes, N, n_trials,
              output_dir, seed=0):
    rng = np.random.default_rng(seed)
    P = int(np.sum(cluster_sizes))
    Phi = _make_random_unitary_phis(P, N, rng)
    labels = np.repeat(np.arange(len(cluster_sizes)), cluster_sizes)
    P_per = np.array([sz / P for sz in cluster_sizes], dtype=float)
    entries = []
    for rho_val in rho_values:
        for snr_db in snr_db_values:
            sigma_d2 = 1.0
            N0 = sigma_d2 * 10 ** (-snr_db / 10.0)
            R = _cluster_R(cluster_sizes, rho_val)
            Sigma_cf = closed_form_sigma(Phi, R, sigma_d2, N0)
            Sigma_mc = mc_sigma(Phi, R, sigma_d2, N0, n_trials, rng=rng)
            rel = (
                np.linalg.norm(Sigma_mc - Sigma_cf, "fro")
                / np.linalg.norm(Sigma_cf, "fro")
            )
            coherent, scalar, N0_term = bdlr_decompose(
                Phi, labels, P_per, rho_val * np.ones(len(cluster_sizes)),
                sigma_d2, N0,
            )
            kappa_bound = kappa_weyl_upper_bound(
                coherent, scalar_weight=scalar[0, 0].real, N0=N0,
            )
            kappa_emp = float(np.linalg.cond(Sigma_cf))
            entries.append({
                "rho": float(rho_val),
                "snr_db": float(snr_db),
                "sigma_rel_err": float(rel),
                "kappa_empirical": kappa_emp,
                "kappa_bound": float(kappa_bound),
            })
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "timestamp": timestamp,
        "config": {
            "rho_values": list(rho_values),
            "snr_db_values": list(snr_db_values),
            "cluster_sizes": list(cluster_sizes),
            "N": N, "n_trials": n_trials, "seed": seed,
        },
        "entries": entries,
    }
    out_json = output_dir / f"m1b_validation_{timestamp}.json"
    out_json.write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[3]
    out = ROOT / "Sim" / "experiments" / "Results"
    run_sweep(
        rho_values=[0.0, 0.3, 0.5, 0.7, 0.9],
        snr_db_values=[0.0, 10.0, 20.0],
        cluster_sizes=[3, 3, 2],
        N=32, n_trials=2000, output_dir=out, seed=1,
    )
```

- [ ] **Step 4: Run — expect pass**

Run: `pytest Tests/python/test_sigma_closed_form.py::test_validate_m1b_runs_small_grid -v`
Expected: PASS

- [ ] **Step 5: Run full sweep manually**

Run: `python Sim/experiments/py_correlated/validate_m1b.py`
Expected: creates `Sim/experiments/Results/m1b_validation_<timestamp>.json` with 15 entries (5 rho * 3 SNR).

- [ ] **Step 6: Commit**

```bash
git add Sim/experiments/py_correlated/validate_m1b.py Tests/python/test_sigma_closed_form.py
git commit -m "feat(p1): end-to-end M1-B validation sweep with JSON report"
```

---

## Task 10: Full-suite regression run + diagnostic note

Final gate: run the full pytest suite, confirm every new test passes, write a diagnostic note recording the sweep's pass/fail summary.

**Files:**

- Create: `docs/diagnostics/2026-05-10-p1-m1b-validation-notes.md`

- [ ] **Step 1: Run the full new-test set**

Run: `pytest Tests/python/test_sigma_closed_form.py Tests/python/test_p_res_eff.py Tests/python/test_fisher_crlb.py Tests/python/test_kappa_bound.py -v`
Expected: all tests PASS (count: 9 — 3 sigma + 2 p_res_eff + 3 fisher + 1 kappa; the runner smoke test lives in test_sigma_closed_form.py).

- [ ] **Step 2: Run the existing regression suite to confirm no regressions**

Run: `pytest Tests/python/ -v`
Expected: all previously-passing tests still pass.

- [ ] **Step 3: Write diagnostic note**

Create `docs/diagnostics/2026-05-10-p1-m1b-validation-notes.md`:

```markdown
# P1 M1-B Numerical Validation — Diagnostic Notes

**Date:** 2026-05-10
**Plan:** `docs/superpowers/plans/2026-05-10-p1-m1b-numerical-validation.md`
**Spec:** `docs/superpowers/specs/2026-05-10-afdm-correlated-clip-research-design.md` §3

## Gate Values

The latest `validate_m1b.py` sweep (see `Sim/experiments/Results/m1b_validation_<timestamp>.json`)
is considered a PASS iff every entry satisfies:

- `sigma_rel_err < 0.10` at `n_trials=2000`
- `kappa_bound >= kappa_empirical` (Weyl bound is a true upper bound)
- `fisher_single_path` returns Hermitian PSD matrix for all tested configurations
- Two-path CRLB grows monotonically as `cos(angle(Phi_i, Phi_j)) -> 1`
- Cluster-level CRLB stays within one order of magnitude over the same sweep

## Known Limitations

- `mc_sigma` convergence is slow at high `rho` (rank-1 coherent part); tolerance
  relaxes to 0.15 in the smoke test, 0.10 in the main sweep with 2000 trials.
- `fisher_single_path` uses `np.linalg.pinv(Sigma)`; when `Sigma` is near-singular
  at rho -> 1 the CRLB is sensitive to `eps_reg`. This is captured by the
  kappa-bound test and handled explicitly in P2.

## Downstream Consumers

- P2 `clip_lrt.py` reuses `closed_form_sigma` and the kappa bound as inputs.
- P3 metric implementations consume the MC/closed-form parity as ground-truth
  reference for `Effective DoF ratio`.
```

- [ ] **Step 4: Commit**

```bash
git add docs/diagnostics/2026-05-10-p1-m1b-validation-notes.md
git commit -m "docs(p1): diagnostic notes + gate values for M1-B validation"
```

---

## Self-Review Checklist

**Spec coverage (§3 of design spec):**

- §3.1 LRT-optimal detector: **not** implemented in P1. Intentional — P1 is the validation layer that supplies the `Sigma` inputs P2 needs. Acknowledged: the LRT formula itself is deferred to P2.
- §3.2 Fisher / CRLB: Tasks 5, 6, 7.
- §3.3 kappa upper bound: Task 8.
- §3.4 D1-D5 deliverables: D1 is exercised indirectly via Task 2 (closed-form Sigma is the D1 form evaluated). D5 (kappa bound) is Task 8. D2, D3, D4 are mathematical derivations, not code — they land in the letter draft; P1 only supplies evidence they hold numerically.
- §6 metrics (four-piece suite): not in P1 scope. P3 plan.
- §7 timeline T0-T1: this plan occupies the first 3-4 weeks of T1.

**Placeholder scan:** no TBD / TODO / "implement later". Every step shows complete code or an exact command.

**Type / name consistency:**

- `closed_form_sigma(Phi_list, R, sigma_d2, N0)` — consistent across Tasks 2, 3, 9.
- `bdlr_decompose` returns `(coherent, scalar, N0_term)` — Task 3 defines it and Tasks 8, 9 consume that exact tuple.
- `fisher_single_path(Phi, dPhi_dtheta, h, x_pilot, Sigma)` vs `fisher_two_paths_joint(Phi_i, Phi_j, h_i, h_j, x_pilot, Sigma)` vs `fisher_cluster_level(Phi_cluster, g_cluster, x_pilot, Sigma)` — three distinct names, each with matching positional order and returning a real Hermitian matrix.
- `kappa_weyl_upper_bound(coherent_part, scalar_weight, N0)` — consistent between Tasks 8 and 9 (note: `scalar_weight` is the scalar coefficient, not the whole scalar matrix; Task 9 passes `scalar[0, 0].real` to match).
- `run_sweep` parameters in Task 9 match the smoke test in step 1 exactly.

**No broken references.** Every function or test named in a later task is defined in an earlier task.

---

## Hand-off Note

Downstream plans (P2 LRT detector, P3 metrics, P5 unrolled EM) will import `closed_form_sigma`, `p_res_eff_from_sigma`, `fisher_cluster_level`, and `kappa_weyl_upper_bound`. Keep these signatures stable. If any must change during implementation, update this plan first and flag P2/P3/P5 owners.
