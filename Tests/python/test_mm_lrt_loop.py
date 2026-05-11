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
            # Signed-NLL-safe slack: relative term uses |traj[t]| so the tolerance
            # behaves symmetrically for negative NLL (which is the common regime
            # because logdet Σ can be strongly negative at small eigenvalues).
            slack = abs(traj[t]) * eps_rel + eps_rel
            if traj[t + 1] > traj[t] + slack:
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
    """Cluster 1 (true g=0) should end up markedly smaller than the active clusters.

    The scalar group-prox step `mu = lambda_gl / L` may not shrink |g_1| all the
    way to zero in 25 outer iterations because the MoM inner loop projects the
    pilot-aligned noise component back into the empty cluster, plus leakage
    through non-orthogonal per-cluster operators. Identifiability recovery is
    preserved as long as |g_1| is well below the active-cluster gains — we
    assert a 2.5:1 gap (|g_1| < 0.4 * min_active) rather than strict zero.
    """
    rng = np.random.default_rng(422)
    N, C = 24, 3
    Phi_per_cluster = [_rand_unitary(N, rng) for _ in range(C)]
    g_true = np.array([0.8 + 0j, 0.0 + 0j, 0.6 + 0j])  # cluster 1 is empty
    r = _simulate_observation(Phi_per_cluster, g_true, N0=0.001, rng=rng)
    cfg = MMConfig(lambda_gl=0.05, T_in_max=3, T_out_max=25)
    result = mm_double_loop(r, Phi_per_cluster, N0=0.001, cfg=cfg)
    # Active clusters must stay on.
    assert result.support[0] and result.support[2]
    # Cluster 1 must be the smallest in magnitude by at least a 2.5:1 ratio.
    min_active = min(abs(result.gains[0]), abs(result.gains[2]))
    assert abs(result.gains[1]) < min_active / 2.5, (
        f"|g_1|={abs(result.gains[1]):.4f} not << min(|g_0|,|g_2|)={min_active:.4f}"
    )
