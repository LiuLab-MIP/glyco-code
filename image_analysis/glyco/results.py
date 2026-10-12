"""Result tables: the numbers behind every summary of the analysis.

Nothing is drawn here. Each table is the numeric content that a figure panel of the
manuscript is made from, written as a CSV so it can be checked, re-plotted with any tool, or
deposited as source data:

    cross-validation          loso_summary, metrics_*, confusion_*, per_case_accuracy_*
    metastatic origin         origin_assignment, origin_negative_control
    discriminant projections  lda_scores_*, lda_explained_variance
    standardised profiles     profile_cancer_type, profile_specimen_primary,
                              profile_section_mouse, profile_specimen_pdac_tma, and the
                              dendrogram merges (*_linkage) of the two clustered ones
    pancreatic components     pdac_pca_scores, pdac_pca_variance, pdac_pca_loadings,
                              pdac_pca_specimen_summary
    second pancreatic cohort  pdac_tma_phenotype_stats (per patch), pdac_tma_phenotype_stats_by_case
                              (per specimen), pdac_tma_first_cohort_agreement
    distribution summaries    violin_group_summary

No classifier is fitted in this module: projections use the models serialised by the
modelling step, so a table can never carry a stale or different model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from . import modelling, mouse, pdac_tma
from .config import GLYCANS, RESULTS_DIR
from .selection import load

# --------------------------------------------------------------------------------------
# Fixed choices that shape a table
# --------------------------------------------------------------------------------------

#: Order in which the per-specimen profile of the nine cancer types is listed.
CANCER_PROFILE_ORDER = ("PDAC", "PCa", "GC", "CRC", "HCC", "BC", "CCA", "EC", "LC")
#: Order of the murine sections in their profile table.
ORGAN_PROFILE_ORDER = ("heart", "liver", "spleen", "lung", "kidney")

#: A discriminant axis is defined only up to direction, so which way round each one comes
#: out is an accident of the fit. These are the directions that put the organs where the
#: published coordinates have them; they were fixed by correlating the five class centroids
#: against those coordinates. The human projections are used as the fit returns them.
MOUSE_AXIS_ORIENTATION = {"LD1": 1, "LD2": -1, "LD3": -1}

#: Where, on the second principal component, the two pancreatic glycophenotypes were
#: originally separated by eye. Reported as the share of each patient's patches above it.
PCA_BOUNDARY = 0.15

#: Suffix naming a metastatic lesion of a given origin in the distribution summaries.
LESION_SUFFIX = "-LM"

PATCH_KEYS = ["specimen_id", "view_id", "grid_row", "grid_col", "label"]


# --------------------------------------------------------------------------------------
# Small numeric helpers (display-neutral)
# --------------------------------------------------------------------------------------


def _row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Standardise each row across the five channels."""
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)


def _global_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Standardise against the spread of the whole table, keeping rows comparable."""
    values = df.to_numpy()
    return (df - values.mean()) / values.std()


def _with_zscores(raw: pd.DataFrame, z: pd.DataFrame) -> pd.DataFrame:
    """Raw means and their standardised values side by side (`z_<glycan>` columns)."""
    return raw.join(z.add_prefix("z_"))


def _linkage(z: pd.DataFrame) -> np.ndarray:
    """Ward linkage on Euclidean distance: the one clustering used for every heatmap."""
    return hierarchy.linkage(pdist(z, metric="euclidean"), method="ward")


def _cluster_order(linkage: np.ndarray) -> list[int]:
    """Row order of a dendrogram: the order a clustered heatmap lists its rows in."""
    return hierarchy.dendrogram(linkage, no_plot=True)["leaves"]


def _linkage_table(z: pd.DataFrame, linkage: np.ndarray) -> pd.DataFrame:
    """The merges of a dendrogram, so it can be drawn with any tool.

    Rows are leaves 0..n-1 (in the order of `z`); each merge creates node n, n+1, ... .
    """
    names = list(z.index)
    return pd.DataFrame({
        "node": range(len(names), len(names) + len(linkage)),
        "left": [names[int(i)] if i < len(names) else f"node {int(i)}" for i in linkage[:, 0]],
        "right": [names[int(i)] if i < len(names) else f"node {int(i)}" for i in linkage[:, 1]],
        "height": linkage[:, 2],
        "size": linkage[:, 3].astype(int),
    })


def drop_iqr_outliers(df: pd.DataFrame, value: str, by: list[str]) -> pd.DataFrame:
    """Remove points beyond 1.5 x IQR within each group.

    A display step only: it is applied to the distribution summaries, to show the spread a
    violin panel would draw, and nowhere else. It never touches the classifier or a
    statistical test. Row order is preserved.
    """
    grouped = df.groupby(by, observed=True)[value]
    q1, q3 = grouped.transform("quantile", 0.25), grouped.transform("quantile", 0.75)
    iqr = q3 - q1
    return df[(df[value] >= q1 - 1.5 * iqr) & (df[value] <= q3 + 1.5 * iqr)].copy()


def _region(df: pd.DataFrame) -> pd.Series:
    """"Cancer" or "Adjacent" (tissue next to the tumour), from `sample_type` rather than a label."""
    return df["sample_type"].map({"tumour": "Cancer", "normal": "Adjacent"}).astype("string")


def _lesion_category(df: pd.DataFrame) -> pd.Series:
    """Tumour tissue keeps its cancer type; a metastatic lesion is suffixed.

    Derived from `sample_type`, so it is a way of naming rows, not of choosing them.
    """
    suffix = df["sample_type"].map({"tumour": "", "metastasis": LESION_SUFFIX}).fillna("")
    return (df["label"].astype(str) + suffix).astype("string")


# --------------------------------------------------------------------------------------
# Cross-validation and origin
# --------------------------------------------------------------------------------------


def cross_validation_tables(cvs: dict[str, modelling.CVResult]) -> dict[str, pd.DataFrame]:
    tables = {"loso_summary": pd.DataFrame([cv.summary() for cv in cvs.values()])}
    for name, cv in cvs.items():
        tables[f"metrics_{name}_loso"] = cv.per_class().reset_index()
        tables[f"confusion_{name}_loso"] = cv.confusion().reset_index()
        tables[f"per_case_accuracy_{name}"] = cv.per_case()
    return tables


def origin_tables(origin: dict) -> dict[str, pd.DataFrame]:
    return {
        "origin_assignment": origin["summary"],
        "origin_negative_control": origin["negative_control"],
    }


# --------------------------------------------------------------------------------------
# Discriminant projections
# --------------------------------------------------------------------------------------


def _scores(bundle: dict, df: pd.DataFrame, orientation: dict[str, int] | None = None) -> pd.DataFrame:
    scores = modelling.project(bundle, df)
    for axis, sign in (orientation or {}).items():
        if axis in scores:
            scores[axis] *= sign
    return df[PATCH_KEYS].join(scores)


def projection_tables(human_wide: pd.DataFrame, mouse_wide: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Discriminant coordinates of every patch, from the serialised models."""
    tables, variance = {}, []
    for analysis in modelling.HUMAN_ANALYSES:
        bundle = modelling.load_model(analysis)
        tables[f"lda_scores_{analysis}"] = _scores(bundle, load(analysis, patches=human_wide))
        variance.append((analysis, bundle))
    bundle = modelling.load_model(mouse.ANALYSIS)
    tables["lda_scores_mouse"] = _scores(bundle, mouse.load(mouse_wide), MOUSE_AXIS_ORIENTATION)
    variance.append((mouse.ANALYSIS, bundle))

    tables["lda_explained_variance"] = pd.DataFrame(
        [{"analysis": name, "axis": f"LD{i + 1}", "explained_variance_ratio": ratio}
         for name, b in variance for i, ratio in enumerate(b["explained_variance_ratio"])]
    )
    return tables


# --------------------------------------------------------------------------------------
# Standardised profiles
# --------------------------------------------------------------------------------------


def _grouped_profile(df: pd.DataFrame, order: tuple[str, ...]) -> pd.DataFrame:
    """Mean profile per specimen, ordered by class, with each row's class beside it.

    Matrix and class column are read off the SAME sorted frame. Taking the matrix before the
    sort and the class column after it is how a table ends up carrying one class's name
    against another class's data.
    """
    profile = df.groupby(["label", "specimen_id"], observed=True)[GLYCANS].mean().reset_index()
    profile["_order"] = pd.Categorical(profile["label"], categories=list(order), ordered=True)
    profile = profile.sort_values(["_order", "specimen_id"]).drop(columns="_order")
    return profile.set_index("specimen_id")


def profile_tables(human_wide: pd.DataFrame, mouse_wide: pd.DataFrame,
                   tma_wide: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Mean intensity profiles, raw and standardised, for the heatmap panels."""
    tables = {}

    # Mean per cancer type, standardised per row, in clustered order (Ward linkage,
    # Euclidean distance on the row-standardised profiles).
    raw = load("primary", patches=human_wide).groupby("label", observed=True)[GLYCANS].mean()
    z = _row_zscore(raw)
    linkage = _linkage(z)
    leaves = _cluster_order(linkage)
    table = _with_zscores(raw, z).iloc[leaves].assign(cluster_position=range(len(leaves)))
    tables["profile_cancer_type"] = table.reset_index()
    tables["profile_cancer_type_linkage"] = _linkage_table(z, linkage)

    # Mean per patient, standardised against the whole table.
    patients = _grouped_profile(load("primary", patches=human_wide), CANCER_PROFILE_ORDER)
    tables["profile_specimen_primary"] = _with_zscores(
        patients[GLYCANS], _global_zscore(patients[GLYCANS])).join(patients["label"]).reset_index()

    # Mean per mouse section, standardised per row.
    sections = _grouped_profile(mouse.load(mouse_wide), ORGAN_PROFILE_ORDER)
    tables["profile_section_mouse"] = _with_zscores(
        sections[GLYCANS], _row_zscore(sections[GLYCANS])).join(sections["label"]).reset_index()

    # Mean per second-cohort specimen, standardised per row, in clustered order (Ward).
    means = pdac_tma.specimen_means(tma_wide)
    z2 = _row_zscore(means[GLYCANS])
    linkage2 = _linkage(z2)
    leaves2 = _cluster_order(linkage2)
    table2 = _with_zscores(means[GLYCANS], z2).join(means["subtype"]).iloc[leaves2]
    tables["profile_specimen_pdac_tma"] = table2.assign(
        cluster_position=range(len(leaves2))).reset_index()
    tables["profile_specimen_pdac_tma_linkage"] = _linkage_table(z2, linkage2)
    return tables


# --------------------------------------------------------------------------------------
# Principal components of the first pancreatic cohort
# --------------------------------------------------------------------------------------


def pca_tables(human_wide: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Principal components of the pancreatic patches (standardised channels).

    These are the components the two glycophenotypes of the first pancreatic cohort were read off.
    """
    df = load("phenotype", patches=human_wide)
    X = StandardScaler().fit_transform(df[GLYCANS].to_numpy(dtype="float64"))
    pca = PCA()
    components = pca.fit_transform(X)
    names = [f"PC{i + 1}" for i in range(components.shape[1])]

    scores = df[PATCH_KEYS + ["subtype"]].join(pd.DataFrame(components, columns=names, index=df.index))
    summary = scores.groupby("specimen_id", observed=True).agg(
        subtype=("subtype", "first"),
        n_patches=("PC2", "size"),
        mean_PC2=("PC2", "mean"),
        share_above_boundary=("PC2", lambda s: float((s > PCA_BOUNDARY).mean())),
    ).reset_index()
    return {
        "pdac_pca_scores": scores,
        "pdac_pca_variance": pd.DataFrame(
            {"component": names, "explained_variance_ratio": pca.explained_variance_ratio_}),
        "pdac_pca_loadings": pd.DataFrame(pca.components_, index=names, columns=GLYCANS)
        .rename_axis("component").reset_index(),
        "pdac_pca_specimen_summary": summary,
    }


# --------------------------------------------------------------------------------------
# Distribution summaries
# --------------------------------------------------------------------------------------


def _describe(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    g = df.groupby(by, observed=True)["intensity"]
    return pd.DataFrame({
        "n": g.size(), "mean": g.mean(), "median": g.median(),
        "q1": g.quantile(0.25), "q3": g.quantile(0.75), "min": g.min(), "max": g.max(),
    }).reset_index()


def _summaries(panel: str, df: pd.DataFrame, group: str, trim_by: list[str]) -> list[pd.DataFrame]:
    """Per-group summaries of a panel's intensities, untrimmed and after the 1.5 x IQR trim."""
    out = []
    for trim, frame in (("none", df), ("iqr_1.5", drop_iqr_outliers(df, "intensity", trim_by))):
        table = _describe(frame, [group, "glycan"]).rename(columns={group: "group"})
        out.append(table.assign(panel=panel, trim=trim))
    return out


def violin_summary(human_long: pd.DataFrame, tma_wide: pd.DataFrame) -> pd.DataFrame:
    """Median, quartiles and counts of the intensity distributions the violin panels draw."""
    parts = []

    adjacent = load("adjacent", long=True, patches=human_long)
    adjacent = adjacent.assign(Region=_region(adjacent))
    adjacent = adjacent.assign(group=adjacent["label"].astype(str) + "|" + adjacent["Region"].astype(str))
    parts += _summaries("tumour_vs_adjacent", adjacent, "group", ["Region", "label", "glycan"])

    tumours = load("primary", long=True, patches=human_long)
    lesions = load("origin", long=True, patches=human_long)
    both = pd.concat([tumours, lesions])
    both = both.assign(category=_lesion_category(both))
    origins = sorted(lesions.label.unique())
    wanted = [c for o in origins for c in (o, o + LESION_SUFFIX)] + ["LC"]
    both = both[both.category.isin(wanted)]
    parts += _summaries("metastasis_profile", both, "category", ["glycan"])

    parts += _summaries("pdac_by_patient", load("phenotype", long=True, patches=human_long),
                        "specimen_id", ["glycan"])
    parts += _summaries("breast_subtype", load("subtype", long=True, patches=human_long),
                        "subtype", ["subtype", "glycan"])

    # Second pancreatic cohort: untrimmed, as the phenotype comparison itself is.
    tma = pdac_tma.tumour_patches(tma_wide).melt(
        id_vars=["specimen_id", "subtype"], value_vars=GLYCANS, var_name="glycan", value_name="intensity")
    parts.append(_describe(tma, ["subtype", "glycan"]).rename(columns={"subtype": "group"})
                 .assign(panel="pdac_tma_phenotype", trim="none"))

    table = pd.concat(parts, ignore_index=True)
    first = ["panel", "trim", "group", "glycan"]
    return table[first + [c for c in table.columns if c not in first]]


# --------------------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------------------


def build_all(
    human_wide: pd.DataFrame,
    human_long: pd.DataFrame,
    human: dict,
    mouse_wide: pd.DataFrame,
    mouse_cv: modelling.CVResult,
    tma_wide: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Every result table, as {file stem: DataFrame}."""
    cvs = {**human["cv"], mouse.ANALYSIS: mouse_cv}
    tables = cross_validation_tables(cvs)
    tables |= origin_tables(human["origin"])
    tables |= projection_tables(human_wide, mouse_wide)
    tables |= profile_tables(human_wide, mouse_wide, tma_wide)
    tables |= pca_tables(human_wide)
    tables["pdac_tma_phenotype_stats"] = pdac_tma.compare_phenotypes(tma_wide).reset_index()
    tables["pdac_tma_phenotype_stats_by_case"] = (
        pdac_tma.compare_phenotypes_by_case(tma_wide).reset_index())
    agreement = pdac_tma.compare_with_first_cohort()
    if agreement is not None:
        tables["pdac_tma_first_cohort_agreement"] = agreement.reset_index()
    tables["violin_group_summary"] = violin_summary(human_long, tma_wide)
    return tables


def write_all(tables: dict[str, pd.DataFrame], directory=RESULTS_DIR) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(directory / f"{name}.csv", index=False)
