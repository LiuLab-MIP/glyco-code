"""Identity and counts of the human tables."""

from glyco import config
from glyco.cohorts import EXPECTED_PRIMARY, HUMAN


def test_specimen_ids_are_unique_and_complete(wide, cases):
    assert cases.specimen_id.is_unique
    assert set(cases.specimen_id) == set(wide.specimen_id)
    assert (wide.groupby("specimen_id").size() > 0).all()
    # the case table's own patch counts must agree with the patch table
    counted = {str(k): int(v) for k, v in wide.groupby("specimen_id").size().items()}
    declared = {str(k): int(v) for k, v in cases.set_index("specimen_id").n_patches.items()}
    assert counted == declared


def test_sample_type_has_exactly_three_values_with_the_expected_sizes(wide, long, cases):
    assert set(wide.sample_type) == set(HUMAN.sample_types)
    got = cases.groupby("sample_type").specimen_id.nunique().to_dict()
    assert got == HUMAN.expected_specimens
    assert sum(got.values()) == 101
    assert set(wide.specimen_id) == set(long.specimen_id)


def test_row_and_view_counts(wide):
    per_view = wide.groupby("view_id").size()
    assert (per_view == config.PATCHES_PER_VIEW).all(), per_view[per_view != config.PATCHES_PER_VIEW].to_dict()
    assert len(wide) == sum(HUMAN.expected_views.values()) * config.PATCHES_PER_VIEW == 41_984
    assert wide.groupby("sample_type").view_id.nunique().to_dict() == HUMAN.expected_views

    tumour = wide[wide.sample_type == "tumour"]
    for label, (patients, views) in EXPECTED_PRIMARY.items():
        sub = tumour[tumour.label == label]
        assert (sub.specimen_id.nunique(), sub.view_id.nunique()) == (patients, views), label

    # every patch sits at a distinct grid position within its view
    assert not wide.duplicated(["view_id", "grid_row", "grid_col"]).any()
    assert wide.grid_row.between(1, config.GRID).all() and wide.grid_col.between(1, config.GRID).all()


def test_subtype_is_non_null_only_for_bc_and_pdac(wide):
    labelled = wide[wide.subtype.notna()]
    assert set(labelled.label) <= {"BC", "PDAC"}
    assert set(labelled[labelled.label == "BC"].subtype) == {"HER2", "LuminalA", "LuminalB", "TNBC"}
    assert set(labelled[labelled.label == "PDAC"].subtype) == {"TypeI", "TypeII"}
    for label, n in (("BC", 13), ("PDAC", 7)):
        sub = wide[(wide.sample_type == "tumour") & (wide.label == label)]
        assert sub.subtype.notna().all(), f"{label} specimens without a subtype"
        assert sub.groupby("specimen_id").subtype.nunique().eq(1).all()
        assert sub.specimen_id.nunique() == n
    # a metastasis or normal block never carries one
    assert wide[wide.sample_type != "tumour"].subtype.isna().all()


def test_every_patch_appears_exactly_once(wide):
    """Breast and pancreatic patches must not be duplicated by a second filing of the images."""
    for label, views in (("BC", 29), ("PDAC", 10)):
        sub = wide[(wide.sample_type == "tumour") & (wide.label == label)]
        assert len(sub) == views * config.PATCHES_PER_VIEW
    assert not wide.duplicated(["view_id", "grid_row", "grid_col"]).any()


def test_identifiers_carry_no_pathology_numbers(wide):
    """Specimen ids are pseudonyms of the form <prefix>-P<nn>."""
    assert wide.specimen_id.str.fullmatch(r"[A-Za-z]+(-[A-Z]+)?-P\d{2}").all()
