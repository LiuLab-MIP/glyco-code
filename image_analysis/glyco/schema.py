"""Table columns, dtypes, and the only readers and writers.

CSV carries no dtypes, so every read passes `DTYPES`. Restoring them rather than letting
pandas infer is what keeps missing values missing: an inferred integer column cannot hold
NaN, and an inferred numpy bool turns NA into False. Either would silently turn "no
measurement" into "a measurement of zero" or "not tumour".
"""

from __future__ import annotations

import pandas as pd

from .cohorts import HUMAN, Cohort
from .config import GLYCANS

METADATA_COLUMNS = [
    "specimen_id",
    "label",
    "sample_type",
    "subtype",
    "view_id",
    "grid_row",
    "grid_col",
    "n_valid_px",
    "is_tumour",
]

WIDE_COLUMNS = METADATA_COLUMNS + GLYCANS
LONG_COLUMNS = METADATA_COLUMNS + ["glycan", "intensity"]

CASE_COLUMNS = [
    "specimen_id",
    "cohort",
    "sample_type",
    "label",
    "subtype",
    "n_views",
    "n_patches",
    "n_tumour_patches",
    "n_blank_patches",
    "has_annotation",
]

DTYPES: dict[str, str] = {
    "specimen_id": "string",
    "label": "string",
    "sample_type": "string",
    "subtype": "string",
    "view_id": "string",
    "grid_row": "int16",
    "grid_col": "int16",
    "n_valid_px": "int32",
    "is_tumour": "boolean",
    **{g: "float64" for g in GLYCANS},
}

LONG_DTYPES: dict[str, str] = {
    **{k: v for k, v in DTYPES.items() if k not in GLYCANS},
    "glycan": "string",
    "intensity": "float64",
}

CASE_DTYPES: dict[str, str] = {
    "specimen_id": "string",
    "cohort": "string",
    "sample_type": "string",
    "label": "string",
    "subtype": "string",
    "n_views": "int64",
    "n_patches": "int64",
    "n_tumour_patches": "Int64",
    "n_blank_patches": "int64",
    "has_annotation": "boolean",
}


def read_wide(cohort: Cohort = HUMAN) -> pd.DataFrame:
    return pd.read_csv(cohort.table("wide"), dtype=DTYPES)


def read_long(cohort: Cohort = HUMAN) -> pd.DataFrame:
    return pd.read_csv(cohort.table("long"), dtype=LONG_DTYPES)


def read_cases(cohort: Cohort = HUMAN) -> pd.DataFrame:
    return pd.read_csv(cohort.table("cases"), dtype=CASE_DTYPES)


def to_long(wide: pd.DataFrame) -> pd.DataFrame:
    """Reshape the patch table to one row per patch per channel.

    The long table is always derived from the wide one, never measured separately: that is
    what guarantees the violin summaries and the classifier read the same numbers.
    """
    long = wide.melt(
        id_vars=METADATA_COLUMNS,
        value_vars=GLYCANS,
        var_name="glycan",
        value_name="intensity",
    )
    return long.astype(LONG_DTYPES)[LONG_COLUMNS]


def write_tables(cohort: Cohort, wide: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    """Write the wide, long and case tables of a cohort and return the long one.

    No `na_rep`: a missing value must round-trip as an empty field.
    """
    cohort.table("wide").parent.mkdir(parents=True, exist_ok=True)
    long = to_long(wide)
    wide.to_csv(cohort.table("wide"), index=False)
    long.to_csv(cohort.table("long"), index=False)
    cases.to_csv(cohort.table("cases"), index=False)
    return long
