"""Result tables: they must carry the right labels, and describe the same numbers as the models."""

import numpy as np
import pandas as pd
import pytest

from glyco import config, modelling, mouse, pdac_tma, results, schema
from glyco.cohorts import HUMAN
from glyco.selection import load

pytestmark = pytest.mark.skipif(
    not (config.MODEL_DIR / modelling.MODEL_FILES["primary"]).exists(),
    reason="models not built yet; run `python -m glyco analyse`",
)


@pytest.fixture(scope="module")
def tables(wide, long, mouse_wide, tma):
    human = modelling.run_human(wide, save=False)
    return results.build_all(wide, long, human, mouse_wide, mouse.cross_validate(mouse_wide), tma)


def test_profile_rows_keep_their_own_label(tables, wide, mouse_wide):
    """Each row of a profile table carries the class whose specimen it holds."""
    for name, source in (("profile_specimen_primary", wide), ("profile_section_mouse", mouse_wide)):
        table = tables[name].set_index("specimen_id")
        truth = source.groupby("specimen_id").label.first()
        assert (table.label == truth.reindex(table.index)).all(), name
    order = tables["profile_section_mouse"].label.drop_duplicates().tolist()
    assert order == list(results.ORGAN_PROFILE_ORDER)


def test_cluster_positions_are_a_permutation(tables):
    for name in ("profile_cancer_type", "profile_specimen_pdac_tma"):
        positions = tables[name].cluster_position
        assert sorted(positions) == list(range(len(positions))), name


def test_one_scores_row_per_patch_and_axis_per_discriminant(tables, wide):
    primary = tables["lda_scores_primary"]
    assert len(primary) == len(load("primary", patches=wide))
    axes = [c for c in primary.columns if c.startswith("LD")]
    assert axes == [f"LD{i + 1}" for i in range(len(axes))] and len(axes) == 5


def test_mouse_discriminant_axes_keep_a_fixed_orientation(tables):
    """A discriminant axis is defined only up to direction, so pin which way each one runs.

    Without this the projection could come out mirrored between library versions. The
    extremes below are those of the reference run (`MOUSE_AXIS_ORIENTATION` applied).
    """
    scores = tables["lda_scores_mouse"]
    centroid = scores.groupby("label")[["LD1", "LD2", "LD3"]].mean()
    assert centroid.LD1.idxmin() == "heart" and centroid.LD1.idxmax() == "kidney"
    assert centroid.LD2.idxmin() == "heart" and centroid.LD2.idxmax() == "spleen"
    assert centroid.LD3.idxmin() == "liver" and centroid.LD3.idxmax() == "lung"


def test_per_case_tables_cover_every_held_out_specimen(tables):
    for name, n in (("primary", 72), ("subtype", 13), ("phenotype", 7), ("mouse", 15)):
        table = tables[f"per_case_accuracy_{name}"]
        assert len(table) == n and table.specimen_id.is_unique, name
        assert table.accuracy.between(0, 1).all()


def test_trimming_only_removes_points_and_never_adds_them(tables):
    summary = tables["violin_group_summary"]
    assert set(summary.trim) == {"none", "iqr_1.5"}
    both = summary.pivot_table(index=["panel", "group", "glycan"], columns="trim", values="n")
    both = both.dropna()
    assert (both["iqr_1.5"] <= both["none"]).all()
    assert set(summary.panel) == {"tumour_vs_adjacent", "metastasis_profile", "pdac_by_patient",
                                  "breast_subtype", "pdac_tma_phenotype"}


def test_phenotype_statistics_are_untrimmed(tables, tma):
    """The reported effect must not depend on what a display summary chose to trim."""
    stats = tables["pdac_tma_phenotype_stats"].set_index("glycan")
    direct = pdac_tma.compare_phenotypes(tma)
    pd.testing.assert_frame_equal(stats, direct, check_names=False)


def test_pca_summary_covers_the_seven_first_cohort_patients(tables):
    summary = tables["pdac_pca_specimen_summary"]
    assert len(summary) == 7 and summary.share_above_boundary.between(0, 1).all()
    assert np.isclose(tables["pdac_pca_variance"].explained_variance_ratio.sum(), 1.0)


def test_origin_summary_lists_every_class_for_every_lesion(tables):
    summary = tables["origin_assignment"]
    assert summary.groupby("label").predicted.nunique().eq(9).all()
    assert tables["origin_negative_control"].patches.sum() == 0


def test_tables_round_trip_through_csv(tables, tmp_path):
    results.write_all(tables, tmp_path)
    for name in ("loso_summary", "per_case_accuracy_primary", "origin_assignment"):
        back = pd.read_csv(tmp_path / f"{name}.csv")
        assert list(back.columns) == list(tables[name].columns), name
        assert len(back) == len(tables[name]), name
    assert schema.read_wide(HUMAN) is not None


def test_clustering_is_ward_on_euclidean_distance_everywhere(tables):
    """One clustering method for every heatmap: Ward linkage on Euclidean distance."""
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import pdist

    for name, linkage_name in (("profile_cancer_type", "profile_cancer_type_linkage"),
                               ("profile_specimen_pdac_tma", "profile_specimen_pdac_tma_linkage")):
        profile, merges = tables[name], tables[linkage_name]
        z = profile[[f"z_{g}" for g in config.GLYCANS]].to_numpy()
        # `profile` is already in clustered order; clustering it again must give a tree whose
        # merge heights are those of the stored linkage
        stored = np.sort(merges.height.to_numpy())
        again = np.sort(hierarchy.linkage(pdist(z, metric="euclidean"), method="ward")[:, 2])
        assert np.allclose(stored, again, rtol=0, atol=1e-9), name
        assert len(merges) == len(profile) - 1 and merges["size"].iloc[-1] == len(profile)
        assert (np.diff(merges.height.to_numpy()) >= -1e-12).all()   # Ward heights never decrease
