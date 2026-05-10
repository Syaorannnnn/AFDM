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
