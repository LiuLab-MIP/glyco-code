"""Missing values survive the round trip."""

import numpy as np
import pandas as pd

from glyco import config, schema
from glyco.selection import load


def test_dtypes_and_nan_counts_survive_a_write_read_round_trip(wide, long, tmp_path):
    for g in config.GLYCANS:
        assert wide[g].dtype == np.dtype("float64"), g
    assert isinstance(wide.is_tumour.dtype, pd.BooleanDtype)

    p = tmp_path / "wide.csv"
    wide.to_csv(p, index=False)          # no na_rep: NaN must round-trip as an empty field
    back = pd.read_csv(p, dtype=schema.DTYPES)
    assert back.dtypes.equals(wide.dtypes)
    assert back.isna().sum().equals(wide.isna().sum())
    pd.testing.assert_frame_equal(back, wide.reset_index(drop=True))

    q = tmp_path / "long.csv"
    long.to_csv(q, index=False)
    back_long = pd.read_csv(q, dtype=schema.LONG_DTYPES)
    assert back_long.dtypes.equals(long.dtypes)
    assert back_long.isna().sum().equals(long.isna().sum())


def test_blank_patches_are_nan_in_all_five_channels_never_zero(wide):
    minimum = wide.sample_type.map({k: v.min_valid_px for k, v in config.SETTINGS.items()})
    blank = wide.n_valid_px < minimum
    assert wide.loc[blank, config.GLYCANS].isna().all(axis=None), "a blank patch kept a value"
    assert wide.loc[~blank, config.GLYCANS].notna().all(axis=None), "a non-blank patch lost its value"
    # the signature of an integer cast: a patch whose five channels are all exactly zero
    assert not (wide[config.GLYCANS] == 0).all(axis=1).any()
    # NaN is all-or-nothing across the five channels
    assert set(wide[config.GLYCANS].isna().sum(axis=1).unique()) <= {0, 5}


def test_is_tumour_is_na_exactly_where_there_is_no_mask(wide):
    normal = wide.sample_type == "normal"
    assert wide.loc[normal, "is_tumour"].isna().all()
    assert wide.loc[~normal, "is_tumour"].notna().all()
    # tumour_only must RETAIN the normal patches: fillna(True), not fillna(False)
    kept = load("adjacent", tumour_only=True, patches=wide)
    n_normal_nonblank = int((normal & wide[config.GLYCANS[0]].notna()).sum())
    assert int((kept.sample_type == "normal").sum()) == n_normal_nonblank > 0
    # and it must still mask the annotated tissue
    assert kept[kept.sample_type == "tumour"].is_tumour.all()


def test_long_melts_exactly_from_wide_and_back(wide, long):
    pd.testing.assert_frame_equal(
        schema.to_long(wide).reset_index(drop=True), long.reset_index(drop=True)
    )
    assert len(long) == 5 * len(wide) == 209_920
    # NaN survives the melt: a blank patch is five NaN rows
    assert int(long.intensity.isna().sum()) == 5 * int(wide[config.GLYCANS[0]].isna().sum())

    # unstack, not pivot_table: a pivot drops rows whose index contains NaN, and `subtype`
    # is NaN for seven of the nine cancer types
    key = ["view_id", "grid_row", "grid_col"]
    back = long.set_index(key + ["glycan"]).intensity.unstack("glycan")[config.GLYCANS]
    front = wide.set_index(key)[config.GLYCANS]
    assert len(back) == len(front) == len(wide)
    pd.testing.assert_frame_equal(
        back.sort_index(), front.sort_index(), check_names=False, check_column_type=False
    )
