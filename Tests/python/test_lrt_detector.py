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
