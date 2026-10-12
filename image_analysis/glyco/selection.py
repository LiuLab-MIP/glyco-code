"""Which rows belong to which analysis.

**The one rule.** Membership of an analysis is decided by `sample_type`, and by nothing
else. Never select rows by a `label` value or by a string match on `specimen_id`. The
metastases carry `label = "CRC"` and `label = "BC"`, the same values as two of the nine
training classes, because that is their correct tissue of origin and is what makes the
origin prediction meaningful. Filtering on `label` would silently fold them back into
training. The two sub-analyses below narrow *within* the tumours by label, which is a
different thing from deciding membership.
"""

from __future__ import annotations

import pandas as pd

from .config import GLYCANS
from .schema import read_long, read_wide

ANALYSES = {
    "primary":   dict(where='sample_type == "tumour"', target="label"),
    "subtype":   dict(where='sample_type == "tumour" and label == "BC"', target="subtype"),
    "phenotype": dict(where='sample_type == "tumour" and label == "PDAC"', target="subtype"),
    "origin":    dict(where='sample_type == "metastasis"', target="label"),
    "adjacent":  dict(where='sample_type in ["tumour", "normal"]', target="label"),
}


def load(analysis: str, *, tumour_only: bool = True, long: bool = False, patches=None) -> pd.DataFrame:
    """Patches of the human cohort for one analysis, with its target renamed to `label`.

    `analysis` is required and has no default, so no call site can receive everything by
    omission.

    `tumour_only` uses the annotated tumour region where an annotation exists and the whole
    section where there is none. The adjacent-normal blocks are unannotated, so their
    `is_tumour` is NA and `fillna(True)` retains them, matching the Methods, which omit the
    masking step for non-malignant tissue. Keeping those blocks out of the models is the job
    of `sample_type`, not of this flag.
    """
    spec = ANALYSES[analysis]
    if patches is None:
        patches = read_long() if long else read_wide()
    df = patches.query(spec["where"])
    df = df.dropna(subset=["intensity"] if long else GLYCANS)   # blank patches
    if tumour_only:
        df = df[df["is_tumour"].fillna(True)]
    return df.assign(label=df[spec["target"]])


def features(df: pd.DataFrame):
    """The (X, y, groups) triple every model is fitted on."""
    return df[GLYCANS].to_numpy(dtype="float64"), df["label"].to_numpy(), df["specimen_id"].to_numpy()
