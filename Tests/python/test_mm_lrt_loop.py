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
