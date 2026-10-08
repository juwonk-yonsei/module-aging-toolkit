"""module_aging -- module- and cell-type-resolved aging-clock toolkit.

Reusable, clock-agnostic primitives: the diagnostic checklist for cross-cohort
non-replication (revised manuscript, Table 2) and, for reproducing v0.1.0, the
attribution ladder of the original submission. The transcriptomic clocks themselves
come from the tAge dependency (see THIRD_PARTY_NOTICES.md).
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
from .checklist import (
    ChecklistResult,
    checklist,
    cross_cohort_test,
    effect_vector,
    split_half_reliability,
)

__all__ = [
    "paths",
    "ChecklistResult",
    "checklist",
    "cross_cohort_test",
    "effect_vector",
    "split_half_reliability",
    "Concordance",
    "LadderVerdict",
    "classify",
    "concordance_permutation_p",
    "effect_concordance",
    "split_half_ceiling",
]

__version__ = "0.2.0"
