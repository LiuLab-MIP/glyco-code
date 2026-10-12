"""Extraction: the pixel rule on synthetic images, and the raw-image tree where it is present."""

import re

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from glyco import config, extraction, schema
from glyco.cohorts import COHORTS, HUMAN, MOUSE, PDAC_TMA
from glyco.config import GLYCANS, ExtractionSettings

RULE = ExtractionSettings("test", valid_threshold=10, saturated_above=250, min_valid_px=512)


def _write_view(directory, channels):
    directory.mkdir(parents=True, exist_ok=True)
    for glycan, value in zip(GLYCANS, channels):
        Image.fromarray(np.asarray(value, dtype=np.uint8)).convert("RGB").save(directory / f"{glycan}.jpg", quality=100)


def _flat(value):
    return np.full((config.IMAGE_SIZE, config.IMAGE_SIZE), value, dtype=np.uint8)


def test_a_pixel_is_kept_when_any_channel_reaches_the_threshold(tmp_path):
    """One bright channel keeps the pixel for all five, so background cannot dilute the mean."""
    dim = _flat(2)
    bright = _flat(2)
    bright[:32, :] = 100                     # half of the first row of patches carries signal
    _write_view(tmp_path, [bright, dim, dim, dim, dim])
    n_valid, means = extraction.patch_statistics(tmp_path, RULE)
    assert n_valid[0, 0] == 32 * config.BLOCK and n_valid[5, 5] == 0
    # JPEG is lossy at the edges, so compare loosely: the GlcNAc mean is about 100, Gal about 2
    assert abs(means["GlcNAc"][0, 0] - 100) < 5 and abs(means["Gal"][0, 0] - 2) < 2


def test_a_patch_with_too_little_signal_is_missing_in_all_five_channels_never_zero(tmp_path):
    _write_view(tmp_path, [_flat(0)] * 5)
    n_valid, means = extraction.patch_statistics(tmp_path, RULE)
    assert (n_valid == 0).all()
    for glycan in GLYCANS:
        assert np.isnan(means[glycan]).all(), glycan


def test_a_pixel_saturated_in_every_channel_is_dropped(tmp_path):
    _write_view(tmp_path, [_flat(255)] * 5)
    n_valid, means = extraction.patch_statistics(tmp_path, RULE)
    assert (n_valid == 0).all() and np.isnan(means["Fuc"]).all()


def test_unrecognised_annotation_colour_stops_the_run(tmp_path):
    arr = np.zeros((264, 264, 3), dtype=np.uint8)
    Image.fromarray(arr).save(tmp_path / "annotation.png")
    with pytest.raises(ValueError, match="Unrecognised colour in annotation"):
        extraction.read_annotation(tmp_path / "annotation.png")


# --------------------------------------------------------------------------------------
# The raw-image tree (skipped when the images are not present)
# --------------------------------------------------------------------------------------

needs_images = pytest.mark.skipif(
    not HUMAN.image_root.is_dir(), reason="raw images not present (see raw_images/README.md)"
)


@needs_images
def test_the_tree_holds_exactly_the_expected_cohorts():
    for cohort in COHORTS.values():
        views = list(extraction.iter_views(cohort))
        for sample_type in cohort.sample_types:
            of_type = [v for v in views if v.sample_type == sample_type]
            assert len({v.specimen_id for v in of_type}) == cohort.expected_specimens[sample_type]
            assert len(of_type) == cohort.expected_views[sample_type]
        for view in views:
            for glycan in GLYCANS:
                assert (view.directory / f"{glycan}.jpg").is_file(), (view.specimen_id, glycan)
            annotated = view.sample_type in cohort.annotated_sample_types
            assert (view.directory / config.ANNOTATION_FILE).is_file() == annotated


@needs_images
def test_the_manifest_agrees_with_the_tree():
    manifest = extraction.read_manifest()
    in_tree = {(c.name, v.specimen_id, v.sample_type, v.label)
               for c in COHORTS.values() for v in extraction.iter_views(c)}
    listed = {(r.cohort, r.specimen_id, r.sample_type, r.label) for r in manifest.itertuples()}
    assert in_tree == listed


@needs_images
def test_no_path_in_the_image_tree_contains_non_latin_characters_or_a_pathology_number():
    for path in config.RAW_IMAGES.rglob("*"):
        relative = str(path.relative_to(config.RAW_IMAGES))
        assert relative.isascii(), relative
    for view in extraction.iter_views(HUMAN):
        assert re.fullmatch(r"[A-Za-z]+(-[A-Z]+)?-P\d{2}", view.specimen_id), view.specimen_id


@needs_images
def test_re_extracting_a_cohort_reproduces_the_stored_tables_exactly():
    """The mouse cohort is small enough to redo from the images on every test run."""
    wide, cases = extraction.build(MOUSE)
    stored = schema.read_wide(MOUSE)
    pd.testing.assert_frame_equal(wide.reset_index(drop=True), stored.reset_index(drop=True),
                                  check_exact=False, rtol=1e-12)
    assert extraction.check_cohort(MOUSE, wide, cases) == []


@needs_images
def test_re_extracting_the_second_pancreatic_cohort_reproduces_the_stored_measurements():
    wide, _ = extraction.build(PDAC_TMA)
    stored = schema.read_wide(PDAC_TMA)
    columns = [c for c in schema.WIDE_COLUMNS if c != "subtype"]   # the phenotype is assigned afterwards
    pd.testing.assert_frame_equal(
        wide[columns].reset_index(drop=True), stored[columns].reset_index(drop=True),
        check_exact=False, rtol=1e-12)
