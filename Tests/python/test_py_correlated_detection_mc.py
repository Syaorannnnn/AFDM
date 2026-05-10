"""经验 Pfa/Pd 验证模块的单元测试。"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Sim" / "experiments" / "py_correlated"))


def test_detection_mc_module_imports():
    import detection_mc
    assert hasattr(detection_mc, "run_detection_threshold_mc")


def test_detection_threshold_mc_outputs_expected_fields():
    from detection_mc import run_detection_threshold_mc
    cfg = {
        "n_sc": 64,
        "n_paths": 6,
        "n_clusters": 2,
        "max_delay": 3,
        "max_doppler": 2,
        "doppler_guard": 3,
        "dirichlet_r": 3,
        "rho": 0.5,
        "mc_samples": 60,
        "seed": 321,
        "p_fa": 1e-3,
    }
    result = run_detection_threshold_mc(cfg)
    assert 1.0 <= result["m_eff_grid"] <= result["m_raw_grid"]
    assert 0.0 <= result["empirical_pfa_independent"] <= 1.0
    assert 0.0 <= result["empirical_pfa_corrected"] <= 1.0
    assert 0.0 <= result["empirical_pd_independent"] <= 1.0
    assert 0.0 <= result["empirical_pd_corrected"] <= 1.0
