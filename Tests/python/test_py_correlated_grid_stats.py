"""候选网格统计量模块的单元测试。"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_grid_stats_module_imports():
    import grid_stats
    assert hasattr(grid_stats, "build_candidate_grid")
    assert hasattr(grid_stats, "estimate_grid_stat_covariance")


def test_build_candidate_grid_size_and_columns():
    from grid_stats import build_candidate_grid
    grid = build_candidate_grid(max_delay=3, max_doppler=2)
    assert grid.shape == (20, 2)
    assert np.array_equal(grid[0], np.array([0, -2]))
    assert np.array_equal(grid[-1], np.array([3, 2]))


def test_grid_covariance_is_reproducible_for_fixed_seed():
    from grid_stats import estimate_grid_stat_covariance
    cfg = {
        "n_sc": 64,
        "n_paths": 6,
        "n_clusters": 2,
        "max_delay": 3,
        "max_doppler": 2,
        "doppler_guard": 3,
        "dirichlet_r": 3,
        "rho": 0.5,
        "mc_samples": 40,
        "seed": 123,
    }
    first = estimate_grid_stat_covariance(cfg)
    second = estimate_grid_stat_covariance(cfg)
    assert first["cov"].shape == (20, 20)
    assert np.allclose(first["cov"], second["cov"])
    assert np.all(first["power_samples"] >= 0.0)
