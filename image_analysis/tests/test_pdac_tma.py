"""The second pancreatic cohort: one annotated field of view per specimen."""


from glyco import config, pdac_tma
from glyco.cohorts import HUMAN, PDAC_TMA, PHENOTYPE_FUC_HIGH, PHENOTYPE_MAN_HIGH


def test_every_specimen_is_one_annotated_field_of_view(tma):
    assert tma.specimen_id.nunique() == PDAC_TMA.expected_specimens["tumour"] == 15
    assert len(tma) == 15 * config.PATCHES_PER_VIEW
    assert (tma.groupby("specimen_id").view_id.nunique() == 1).all()
    assert tma.is_tumour.notna().all()
    assert tma.is_tumour.any() and not tma.is_tumour.all()


def test_phenotype_follows_the_stated_rule(tma):
    """Fucose above mannose, on the specimen's own annotated tumour patches."""
    tumour = tma[tma.is_tumour.fillna(False)]
    means = tumour.groupby("specimen_id")[["Fuc", "Man"]].mean()
    declared = tumour.groupby("specimen_id").subtype.first()
    expected = means.Fuc.gt(means.Man).map({True: PHENOTYPE_FUC_HIGH, False: PHENOTYPE_MAN_HIGH})
    assert (declared == expected).all()
    assert set(declared) == {PHENOTYPE_FUC_HIGH, PHENOTYPE_MAN_HIGH}


def test_the_two_pancreatic_cohorts_stay_apart(tma, wide):
    """Different cohorts, different labelling rules, so never one table."""
    assert PDAC_TMA.table("wide") != HUMAN.table("wide")
    assert not set(tma.specimen_id) & set(wide.specimen_id)
    # the phenotype vocabularies must not be confusable either
    assert not set(tma.subtype.dropna()) & set(wide.subtype.dropna())


def test_the_two_channel_rule_recovers_the_first_cohort_labels(wide):
    """If it did not, the two cohorts would not be describing the same distinction."""
    agreement = pdac_tma.compare_with_first_cohort()
    assert agreement is not None and (agreement.labelled == agreement.rule_would_say).all()


def test_the_separation_is_carried_by_fucose_and_mannose(tma):
    stats = pdac_tma.compare_phenotypes(tma)
    assert stats.loc["Fuc", "rank_biserial_r"] > 0.5
    assert stats.loc["Man", "rank_biserial_r"] < -0.5
    for glycan in ("GlcNAc", "Gal", "SA"):
        assert abs(stats.loc[glycan, "rank_biserial_r"]) < 0.2, glycan


def test_the_tma_cohort_is_measured_like_the_rest_of_the_malignant_tissue():
    """Including the saturation cut, which applies to every tissue."""
    assert config.SETTINGS["tumour"] is config.CANCER
    assert config.CANCER.saturated_above == config.SATURATED_ABOVE == 250


def test_the_specimen_level_comparison_uses_one_value_per_specimen(tma):
    by_case = pdac_tma.compare_phenotypes_by_case(tma)
    assert (by_case.n_fuc_high == 9).all() and (by_case.n_man_high == 6).all()
    # a two-sided exact test on 9 against 6 cannot go below 2 / C(15, 6)
    assert by_case.p_value.min() >= 2 / 5005 - 1e-12
    # fucose and mannose define the phenotype, so they separate the specimens; the rest do not
    assert by_case.loc["Fuc", "p_value"] < 0.01 and by_case.loc["Man", "p_value"] < 0.01
    for glycan in ("GlcNAc", "Gal", "SA"):
        assert by_case.loc[glycan, "p_value"] > 0.5, glycan
