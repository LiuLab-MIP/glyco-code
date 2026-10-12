"""The model sees what it should and nothing else, and reproduces the reference numbers."""

import numpy as np
import pytest

from glyco import config, modelling
from glyco.cohorts import NEGATIVE_CONTROL
from glyco.selection import load


@pytest.fixture(scope="module")
def primary(wide):
    return modelling.cross_validate("primary", patches=wide)


def test_every_fold_holds_out_exactly_one_patient(primary, wide):
    """One fold per patient, and no patient on both sides of any fold."""
    assert primary.n_folds == primary.n_specimens == 72
    assert set(primary.groups) == set(load("primary", patches=wide).specimen_id)


def test_cross_validation_is_deterministic(wide):
    a = modelling.cross_validate("phenotype", patches=wide)
    b = modelling.cross_validate("phenotype", patches=wide)
    assert np.array_equal(a.y_pred, b.y_pred)


def test_confusion_matrix_rows_are_the_reference(primary):
    cm = primary.confusion()
    assert list(cm.index) == list(cm.columns) == primary.labels
    assert (cm.sum(axis=1).to_numpy() == primary.per_class()["support"].to_numpy()).all()
    assert np.allclose(primary.confusion(normalise=True).sum(axis=1), 100)


def test_the_final_model_is_fitted_on_every_tumour_specimen_and_only_those(wide):
    bundle = modelling.fit_final("primary", patches=wide, save=False)
    tumour = wide[wide.sample_type == "tumour"]
    assert set(bundle["training_specimens"]) == set(tumour.specimen_id)
    assert len(bundle["training_specimens"]) == 72
    assert bundle["features"] == config.GLYCANS
    assert len(bundle["classes"]) == 9


def test_no_metastatic_or_normal_specimen_is_ever_fitted(wide):
    for analysis in modelling.HUMAN_ANALYSES:
        bundle = modelling.fit_final(analysis, patches=wide, save=False)
        seen = wide[wide.specimen_id.isin(bundle["training_specimens"])]
        assert set(seen.sample_type) == {"tumour"}, analysis


def test_metastases_are_external_to_the_model_that_scores_them(wide):
    bundle = modelling.fit_final("primary", patches=wide, save=False)
    per_patch, summary = modelling.predict_origin(bundle, patches=wide)
    assert not set(per_patch.specimen_id) & set(bundle["training_specimens"])
    assert per_patch.specimen_id.nunique() == 4
    # every class appears in the summary, so a zero is visible rather than absent
    for label in per_patch.label.unique():
        assert set(summary[summary.label == label].predicted) == set(bundle["classes"])


def test_the_negative_control_covers_every_lesion(wide):
    """Each lesion is reported against the organ it was resected from, and none is assigned there."""
    bundle = modelling.fit_final("primary", patches=wide, save=False)
    per_patch, summary = modelling.predict_origin(bundle, patches=wide)
    control = modelling.negative_control(summary)
    assert set(control.origin) == set(per_patch.label)
    assert (control.predicted == NEGATIVE_CONTROL).all()
    assert control.patches.sum() == 0


def test_projection_has_one_axis_per_discriminant(wide):
    bundle = modelling.fit_final("primary", patches=wide, save=False)
    df = load("primary", patches=wide)
    scores = modelling.project(bundle, df)
    assert len(scores) == len(df)
    assert list(scores.columns) == [f"LD{i + 1}" for i in range(scores.shape[1])]
    assert scores.shape[1] == min(len(bundle["classes"]) - 1, len(config.GLYCANS))


def test_thin_classes_are_flagged_not_silently_reported(wide):
    """13 breast and 7 pancreatic specimens is few; the report must say so itself."""
    for analysis in ("subtype", "phenotype"):
        cv = modelling.cross_validate(analysis, patches=wide)
        per_class = load(analysis, patches=wide).groupby("label").specimen_id.nunique()
        if (per_class < 3).any():
            assert cv.warnings, f"{analysis}: thin class not flagged"


def test_headline_numbers_of_the_reference_run(wide):
    """Leave-one-patient-out results recorded in results/RESULTS.md."""
    expected = {"primary": (0.8924, 0.8768), "subtype": (0.6893, 0.5825), "phenotype": (0.9976, 0.9952)}
    for analysis, (accuracy, kappa) in expected.items():
        cv = modelling.cross_validate(analysis, patches=wide)
        assert round(cv.accuracy, 4) == accuracy, analysis
        assert round(cv.kappa, 4) == kappa, analysis


def test_per_case_table_has_one_row_per_specimen(primary):
    table = primary.per_case()
    assert list(table.columns) == ["specimen_id", "label", "n_patches", "n_correct", "accuracy"]
    assert len(table) == 72 and table.specimen_id.is_unique
    assert np.allclose(table.accuracy, table.n_correct / table.n_patches)


def test_specificity_is_one_versus_rest(primary):
    """Specificity of a class = of the patches not in it, the share not assigned to it."""
    table = primary.per_class()
    assert list(table.columns) == ["precision", "recall", "f1-score", "support", "specificity"]
    for label in primary.labels:
        truth, predicted = primary.y_true == label, primary.y_pred == label
        expected = (~truth & ~predicted).sum() / (~truth).sum()
        assert np.isclose(table.loc[label, "specificity"], expected, rtol=0, atol=1e-12), label
    assert table.specificity.between(0, 1).all()
