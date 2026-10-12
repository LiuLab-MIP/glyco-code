"""The second pancreatic cohort (tissue microarray).

A tissue microarray (TMA) of annotated pancreatic cases, one field of view each, measured with the
same rule as the rest of the malignant tissue and kept in its own tables.

**How the glycophenotype is assigned here.** Per specimen, over the annotated tumour patches:
mean fucose above mean mannose is Fuc-high, otherwise Man-high. That is a rule on two
channels, and it is not how the first pancreatic cohort was labelled: those labels came from a
principal-component analysis of all five channels. The two cohorts are therefore kept apart
rather than pooled, and `compare_with_first_cohort` reports whether the simple rule would
have recovered the first cohort's labels.
"""

from __future__ import annotations

import pandas as pd
from scipy.stats import mannwhitneyu

from . import schema
from .cohorts import PDAC_TMA, PHENOTYPE_FUC_HIGH, PHENOTYPE_MAN_HIGH
from .config import GLYCANS, PATCHES_PER_VIEW
from .selection import load

PHENOTYPES = (PHENOTYPE_FUC_HIGH, PHENOTYPE_MAN_HIGH)


# --------------------------------------------------------------------------------------
# Phenotype
# --------------------------------------------------------------------------------------


def assign_phenotype(tumour_patches: pd.DataFrame) -> pd.Series:
    """specimen_id -> phenotype, from the specimen's own mean fucose and mannose.

    Taken over the annotated tumour patches only, so the surrounding tissue cannot tip a
    specimen from one phenotype into the other.
    """
    means = tumour_patches.groupby("specimen_id", observed=True)[["Fuc", "Man"]].mean()
    phenotype = means.Fuc.gt(means.Man).map({True: PHENOTYPE_FUC_HIGH, False: PHENOTYPE_MAN_HIGH})
    return phenotype.rename("subtype").astype("string")


def annotate_phenotype(wide: pd.DataFrame, cases: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach the assigned phenotype to the freshly extracted tables."""
    phenotype = assign_phenotype(wide[wide.is_tumour.fillna(False)])
    wide = wide.assign(subtype=wide.specimen_id.map(phenotype).astype("string"))
    cases = cases.assign(subtype=cases.specimen_id.map(phenotype).astype("string"))
    return wide, cases


# --------------------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------------------


def read_wide() -> pd.DataFrame:
    return schema.read_wide(PDAC_TMA)


def read_long() -> pd.DataFrame:
    return schema.read_long(PDAC_TMA)


def read_cases() -> pd.DataFrame:
    return schema.read_cases(PDAC_TMA)


def tumour_patches(wide: pd.DataFrame | None = None) -> pd.DataFrame:
    """The annotated tumour patches, which is what every summary here uses."""
    df = read_wide() if wide is None else wide
    return df[df.is_tumour.fillna(False)].dropna(subset=GLYCANS)


def specimen_means(wide: pd.DataFrame | None = None) -> pd.DataFrame:
    """Mean intensity per specimen, with its phenotype."""
    df = tumour_patches(wide)
    means = df.groupby("specimen_id", observed=True)[GLYCANS].mean()
    return means.join(df.groupby("specimen_id", observed=True).subtype.first())


def check(wide: pd.DataFrame, cases: pd.DataFrame) -> list[str]:
    problems = []
    expected = PDAC_TMA.expected_specimens["tumour"]
    if len(cases) != expected:
        problems.append(f"{len(cases)} specimens, expected {expected}")
    if len(wide) != expected * PATCHES_PER_VIEW:
        problems.append(f"{len(wide):,} patches, expected {expected * PATCHES_PER_VIEW:,}")
    if wide.subtype.isna().any():
        problems.append("a specimen was left without a phenotype")
    if not set(wide.subtype.dropna()) <= set(PHENOTYPES):
        problems.append(f"unexpected phenotype values: {set(wide.subtype.dropna())}")
    return problems


# --------------------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------------------


def compare_phenotypes(wide: pd.DataFrame | None = None) -> pd.DataFrame:
    """Mann-Whitney P and rank-biserial r between the two phenotypes, per channel.

    Computed on the patches as selected, before any display-only trimming, so the reported
    effect is not an artefact of what a summary chose to show.
    """
    df = tumour_patches(wide)
    rows = []
    for glycan in GLYCANS:
        a = df.loc[df.subtype == PHENOTYPE_FUC_HIGH, glycan]
        b = df.loc[df.subtype == PHENOTYPE_MAN_HIGH, glycan]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({
            "glycan": glycan,
            "fuc_high_median": a.median(),
            "man_high_median": b.median(),
            "p_value": p,
            # rank-biserial r: how often a patch of one phenotype exceeds one of the other,
            # rescaled so 0 is no difference and +/-1 is complete separation
            "rank_biserial_r": 2 * u / (len(a) * len(b)) - 1,
        })
    return pd.DataFrame(rows).set_index("glycan")


def compare_phenotypes_by_case(wide: pd.DataFrame | None = None) -> pd.DataFrame:
    """The same comparison with the specimen, not the patch, as the unit of observation.

    The patches of one specimen are not independent of each other, so the patch-level test
    above overstates the evidence: it treats 1,628 and 581 patches as separate observations
    although they come from 9 and 6 specimens. Here each specimen contributes one value per
    channel (the mean over its annotated tumour patches), and the two groups of specimens are
    compared with a two-sided exact Mann-Whitney test, which is valid for small groups
    without ties. The rank-biserial r is computed in the same way as at patch level.

    Fucose and mannose are the two channels that DEFINE the phenotype (mean fucose above
    mean mannose), so their separation is partly built in; the other three channels are the
    independent part of the comparison.
    """
    df = tumour_patches(wide)
    case = df.groupby("specimen_id", observed=True)[GLYCANS].mean().join(
        df.groupby("specimen_id", observed=True).subtype.first())
    rows = []
    for glycan in GLYCANS:
        a = case.loc[case.subtype == PHENOTYPE_FUC_HIGH, glycan]
        b = case.loc[case.subtype == PHENOTYPE_MAN_HIGH, glycan]
        u, p = mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rows.append({
            "glycan": glycan,
            "n_fuc_high": len(a),
            "n_man_high": len(b),
            "fuc_high_median": a.median(),
            "man_high_median": b.median(),
            "u_statistic": u,
            "p_value": p,
            "rank_biserial_r": 2 * u / (len(a) * len(b)) - 1,
        })
    return pd.DataFrame(rows).set_index("glycan")


def compare_with_first_cohort() -> pd.DataFrame | None:
    """Would the fucose-against-mannose rule have recovered the first cohort's labels?

    The first cohort was labelled from a principal-component analysis. If the simple rule
    agrees with it, the two cohorts are at least describing the same distinction.
    """
    from .cohorts import HUMAN

    if not HUMAN.table("wide").exists():
        return None
    first = load("phenotype")
    means = first.groupby("specimen_id", observed=True)[["Fuc", "Man"]].mean()
    return pd.DataFrame({
        "labelled": first.groupby("specimen_id", observed=True).subtype.first(),
        "rule_would_say": means.Fuc.gt(means.Man).map({True: "TypeI", False: "TypeII"}),
        "Fuc": means.Fuc.round(1),
        "Man": means.Man.round(1),
    })
