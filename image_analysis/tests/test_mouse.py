"""The murine organ cohort: its own tables, its own held-out unit."""

import numpy as np

from glyco import config, mouse
from glyco.cohorts import HUMAN, MOUSE, MOUSE_ORGANS, SECTIONS_PER_ORGAN


def test_five_organs_three_sections_each(mouse_wide, mouse_cases):
    assert len(mouse_cases) == len(MOUSE_ORGANS) * SECTIONS_PER_ORGAN == 15
    per_organ = mouse_cases.groupby("label").size()
    assert set(per_organ.index) == set(MOUSE_ORGANS)
    assert (per_organ == SECTIONS_PER_ORGAN).all()
    assert len(mouse_wide) == 15 * config.PATCHES_PER_VIEW == 3840
    assert (mouse_wide.groupby("view_id").size() == config.PATCHES_PER_VIEW).all()


def test_sections_are_imaged_whole_so_nothing_is_masked_or_missing(mouse_wide):
    """No annotation and no dim-patch rule, so every patch carries a measurement."""
    assert mouse_wide.is_tumour.isna().all()
    assert mouse_wide.subtype.isna().all()
    assert mouse_wide[config.GLYCANS].notna().all(axis=None)
    assert set(mouse_wide.sample_type) == set(MOUSE.sample_types)


def test_the_two_cohorts_never_share_a_table_or_an_identifier(mouse_wide, wide):
    """A murine organ must not be reachable from a human analysis, or the reverse."""
    assert MOUSE.table("wide") != HUMAN.table("wide")
    assert not set(mouse_wide.specimen_id) & set(wide.specimen_id)
    assert not set(mouse_wide.label) & set(wide.label)
    assert not set(mouse_wide.sample_type) & set(wide.sample_type)


def test_every_fold_holds_out_exactly_one_section(mouse_wide):
    cv = mouse.cross_validate(mouse_wide)
    assert cv.n_folds == cv.n_specimens == 15
    assert set(cv.groups) == set(mouse_wide.specimen_id)
    # each section is one whole image, so a fold is a section and an animal at once
    assert (mouse_wide.groupby("specimen_id").view_id.nunique() == 1).all()


def test_cross_validation_is_deterministic(mouse_wide):
    a, b = mouse.cross_validate(mouse_wide), mouse.cross_validate(mouse_wide)
    assert np.array_equal(a.y_pred, b.y_pred)


def test_reference_accuracy(mouse_wide):
    cv = mouse.cross_validate(mouse_wide)
    assert (round(cv.accuracy, 4), round(cv.kappa, 4)) == (0.9974, 0.9967)


def test_the_final_model_sees_every_section_and_only_sections(mouse_wide):
    bundle = mouse.fit_final(mouse_wide, save=False)
    assert set(bundle["training_specimens"]) == set(mouse_wide.specimen_id)
    assert len(bundle["training_specimens"]) == 15
    assert set(bundle["classes"]) == set(MOUSE_ORGANS)
    assert bundle["features"] == config.GLYCANS


def test_saturated_pixels_are_dropped_in_every_tissue():
    """Every rule, murine and human, drops a pixel saturated in all five channels."""
    assert config.MOUSE_SECTION.saturated_above == 250
    assert all(rule.saturated_above == 250 for rule in config.SETTINGS.values())
