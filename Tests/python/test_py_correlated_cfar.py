"""CFAR 有效自由度与门限修正比的单元测试。"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))

from cfar import compute_effective_candidates, compute_threshold_ratio, bootstrap_m_eff


# ---- 基础退化测试 ----

def test_effective_candidates_independent_covariance_equals_raw_size():
    cov = np.eye(7)
    assert compute_effective_candidates(cov) == 7.0


def test_effective_candidates_rank_one_covariance_equals_one():
    cov = np.ones((7, 7))
    assert compute_effective_candidates(cov) == 1.0


def test_threshold_ratio_uses_effective_candidate_count():
    cov = np.ones((7, 7))
    ratio, m_eff, m_raw = compute_threshold_ratio(cov, p_fa=1e-3)
    expected = np.log(1.0 / 1e-3) / np.log(7.0 / 1e-3)
    assert m_eff == 1.0
    assert m_raw == 7
    assert np.isclose(ratio, expected)


# ---- 输入验证测试 ----

def test_effective_candidates_rejects_non_square_matrix():
    with pytest.raises(ValueError, match="square"):
        compute_effective_candidates(np.ones((2, 3)))


# ---- bootstrap 区间测试 ----

def test_bootstrap_m_eff_returns_ordered_interval():
    rng = np.random.default_rng(123)
    samples = rng.exponential(scale=1.0, size=(60, 5))
    summary = bootstrap_m_eff(samples, n_boot=20, seed=456)
    assert summary["p05"] <= summary["p50"] <= summary["p95"]
    assert 1.0 <= summary["mean"] <= 5.0
