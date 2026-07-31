"""Central path resolution for the module-aging toolkit.

All locations are overridable by environment variable so the same code runs
from a git checkout, a CI runner, or a relocated data drive.

Env overrides
-------------
MAT_ROOT       repo root (default: two levels above this file)
MAT_DATA       raw/derived data dir              (default: $MAT_ROOT/data)
MAT_RESULTS    analysis outputs / source CSVs    (default: $MAT_ROOT/results)
MAT_FIGURES    figure outputs                    (default: $MAT_ROOT/figures)
TAGE_DIR       Gladyshev-Lab/tAge checkout        (default: $MAT_ROOT/third_party/tAge)
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(key: str, default: Path) -> Path:
    val = os.environ.get(key)
    return Path(val) if val else default


ROOT: Path = _env_path("MAT_ROOT", Path(__file__).resolve().parents[2])

DATA: Path = _env_path("MAT_DATA", ROOT / "data")
RESULTS: Path = _env_path("MAT_RESULTS", ROOT / "results")
FIGURES: Path = _env_path("MAT_FIGURES", ROOT / "figures")

MODULE_CLOCKS: Path = DATA / "module_clocks"
CLOCK_MODELS: Path = DATA / "clock_models"
MIL: Path = RESULTS / "mil"

TAGE_DIR: Path = _env_path("TAGE_DIR", ROOT / "third_party" / "tAge")
TAGE_PY: Path = TAGE_DIR / "inst" / "python"
TAGE_EXTDATA: Path = TAGE_DIR / "inst" / "extdata"

THIRD_PARTY: Path = ROOT / "third_party"


def add_import_paths() -> None:
    """Make the tAge python API and our third_party port importable."""
    import sys

    for p in (str(TAGE_PY), str(THIRD_PARTY)):
        if p not in sys.path:
            sys.path.insert(0, p)
