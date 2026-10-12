"""The murine organ cohort.

A separate experiment from the human cohort, kept in separate tables throughout: five organs,
three sections each, imaged whole with no pathologist's annotation.

**Why the section, not the animal, is the held-out unit.** Each mouse contributed one
section per organ, so within any one organ no two sections come from the same animal: a
held-out section shares no animal with any training section of its own class, and
individual-level information cannot leak within a class. That is the opposite of the human
cohort, where two or three fields of view come from a single section and the patient has to
be the unit. The difference between the two schemes is correct, not an inconsistency.
"""

from __future__ import annotations

import pandas as pd

from . import modelling, schema
from .cohorts import MOUSE, MOUSE_ORGANS, SECTIONS_PER_ORGAN
from .config import GLYCANS, PATCHES_PER_VIEW

ANALYSIS = "mouse"


def read_wide() -> pd.DataFrame:
    return schema.read_wide(MOUSE)


def read_cases() -> pd.DataFrame:
    return schema.read_cases(MOUSE)


def load(wide: pd.DataFrame | None = None) -> pd.DataFrame:
    """The patches the classifier sees: every one of them.

    There is no annotation to mask by and no dim-patch rule, so a section contributes all
    256 of its patches.
    """
    return (read_wide() if wide is None else wide).dropna(subset=GLYCANS)


def cross_validate(wide: pd.DataFrame | None = None) -> modelling.CVResult:
    """Leave-one-section-out, pooling the out-of-fold predictions."""
    return modelling.cross_validate_frame(load(wide), ANALYSIS)


def fit_final(wide: pd.DataFrame | None = None, *, save: bool = True) -> dict:
    """Fit on all fifteen sections and serialise, for the projections to read."""
    bundle = modelling.fit_bundle(load(wide), ANALYSIS)
    if save:
        modelling.save_model(bundle)
    return bundle


def check(wide: pd.DataFrame, cases: pd.DataFrame) -> list[str]:
    """Compare what was extracted against what the tree should contain."""
    problems = []
    if len(cases) != MOUSE.expected_specimens["section"]:
        problems.append(f"{len(cases)} sections, expected {MOUSE.expected_specimens['section']}")
    per_organ = cases.groupby("label").size().to_dict()
    for organ in MOUSE_ORGANS:
        got = per_organ.get(organ, 0)
        if got != SECTIONS_PER_ORGAN:
            problems.append(f"{organ}: {got} sections, expected {SECTIONS_PER_ORGAN}")
    expected_rows = MOUSE.expected_specimens["section"] * PATCHES_PER_VIEW
    if len(wide) != expected_rows:
        problems.append(f"{len(wide):,} patches, expected {expected_rows:,}")
    if wide[GLYCANS].isna().any(axis=None):
        problems.append(f"{int(wide[GLYCANS[0]].isna().sum())} patches with no signal at all")
    return problems
