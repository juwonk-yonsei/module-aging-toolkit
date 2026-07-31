"""module_aging -- module- and cell-type-resolved aging-clock toolkit.

Reusable, clock-agnostic primitives for the manuscript's attribution ladder
(Table 1), plus path configuration. The transcriptomic clocks themselves come
from the tAge dependency (see THIRD_PARTY_NOTICES.md); this package holds the
original analysis method contributed by the paper.
"""
from __future__ import annotations

from . import paths
from .attribution_ladder import (
    Concordance,
    LadderVerdict,
    classify,
    concordance_permutation_p,
    effect_concordance,
    split_half_ceiling,
)

__all__ = [
    "paths",
    "Concordance",
    "LadderVerdict",
    "classify",
    "concordance_permutation_p",
    "effect_concordance",
    "split_half_ceiling",
]

__version__ = "0.1.0"
