"""Leakage and isolation.

Two kinds of check. The behavioural ones assert which rows actually reach a fit; the source
checks assert that the code cannot reintroduce the known defects by hand: a random split, a
seed loop, or selecting rows by a label or an identifier string.
"""

import re

import numpy as np
import pytest
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import LeaveOneGroupOut

from glyco import config
from glyco.selection import load

#: Modules that fit models or compute results; they must select rows through `load()`.
SELECTING_MODULES = ("modelling.py", "results.py")


def test_no_specimen_appears_in_both_train_and_test_in_any_fold(wide):
    df = load("primary", patches=wide)
    groups = df.specimen_id.to_numpy()
    n = 0
    for train_idx, test_idx in LeaveOneGroupOut().split(df, df.label, groups):
        assert not set(groups[train_idx]) & set(groups[test_idx])
        assert len(set(groups[test_idx])) == 1
        n += 1
    assert n == df.specimen_id.nunique() == 72
    # no field of view is split across the boundary either, which patient grouping implies
    assert df.groupby("view_id").specimen_id.nunique().eq(1).all()


def test_no_metastasis_or_normal_row_reaches_a_fit(wide):
    for analysis in ("primary", "subtype", "phenotype"):
        df = load(analysis, patches=wide)
        assert set(df.sample_type) == {"tumour"}, analysis
        assert df[config.GLYCANS].notna().all(axis=None), f"{analysis}: NaN would reach the fit"
        LinearDiscriminantAnalysis(solver="svd").fit(df[config.GLYCANS], df.label)

    origin = load("origin", patches=wide)
    assert set(origin.sample_type) == {"metastasis"}
    assert origin.specimen_id.nunique() == 4
    # external by construction: the training selection asks for tumours
    assert not set(origin.specimen_id) & set(load("primary", patches=wide).specimen_id)
    # the metastases keep their true origin labels, which collide with training classes:
    # that is correct, and is exactly why filtering on `label` would fold them back in
    assert set(origin.label) == {"CRC", "BC"}


def test_sklearn_refuses_a_blank_patch_rather_than_training_on_a_zero(wide):
    df = load("primary", patches=wide).copy()
    X = df[config.GLYCANS].to_numpy(dtype=float)
    X[0, 0] = np.nan
    with pytest.raises(ValueError):
        LinearDiscriminantAnalysis(solver="svd").fit(X, df.label)


def test_no_row_is_selected_by_a_label_value_or_a_specimen_id_string_match(source_files):
    """Membership of an analysis is decided by `sample_type`, and by nothing else.

    `ANALYSES` in selection.py narrows WITHIN tumours by label (the breast and pancreas
    sub-analyses). What is banned is deciding membership by label, and any string match on
    `specimen_id`.
    """
    banned_id = re.compile(r"specimen_id[^\n]*\.str\.(contains|startswith|endswith|match)")
    banned_prefix = re.compile(r"\.str\.(contains|startswith|endswith|match)\([^)]*specimen_id")
    for path in source_files:
        text = path.read_text()
        assert not banned_id.search(text), f"{path.name}: string match on specimen_id"
        assert not banned_prefix.search(text), f"{path.name}: string match on specimen_id"
        for literal in ("CRC-LM", "BC-LM"):
            assert literal not in text, f"{path.name}: metastasis selected by the literal {literal}"
        if path.name in SELECTING_MODULES:
            assert "sample_type ==" not in text, f"{path.name}: hand-rolled selection, use load()"


def test_no_random_split_anywhere(source_files):
    banned = [
        "train_test_split", "ShuffleSplit", "StratifiedShuffleSplit",
        "np.random", "numpy.random", "random.seed", "random_state=",
        "shuffle=True", "set.seed",
    ]
    for path in source_files:
        text = path.read_text()
        for token in banned:
            assert token not in text, f"{path.name}: {token}"
        # no seed loop of the kind that searches for a favourable split
        assert not re.search(r"seed\s*(\+=|=\s*seed\s*\+)", text), f"{path.name}: seed loop"


def test_results_module_fits_no_classifier(source_files):
    """Result tables read the serialised models; none may fit one of its own."""
    text = next(p for p in source_files if p.name == "results.py").read_text()
    assert "LinearDiscriminantAnalysis" not in text
    assert "new_estimator" not in text and "fit_bundle" not in text and ".fit(" not in text
