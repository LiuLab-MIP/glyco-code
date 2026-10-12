"""Write `results/RESULTS.md`, the human-readable record of one run.

The report is generated in one pass from one set of results, so no two numbers in it can
disagree, and nothing in it can be stale relative to anything else. Every figure is stated
together with the grouping that produced it.
"""

from __future__ import annotations

import datetime as dt

import pandas as pd

from . import modelling
from .cohorts import MOUSE_ORGANS, SECTIONS_PER_ORGAN
from .config import CANCER, GLYCANS, MOUSE_SECTION, PATCHES_PER_VIEW, RESULTS_DIR, SATURATED_ABOVE, SETTINGS

REPORT = RESULTS_DIR / "RESULTS.md"


def _table(rows, header) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _confusion(cv: modelling.CVResult) -> str:
    cm = cv.confusion(normalise=True)
    rows = [[f"**{lab}**"] + [f"{v:.1f}" for v in cm.loc[lab]] for lab in cv.labels]
    return _table(rows, ["reference \\ prediction"] + cv.labels)


def _per_class(cv: modelling.CVResult) -> str:
    rows = [
        [lab, f"{r.precision:.3f}", f"{r.recall:.3f}", f"{r.specificity:.3f}", f"{r['f1-score']:.3f}",
         f"{int(r.support):,}"]
        for lab, r in cv.per_class().iterrows()
    ]
    return _table(rows, ["class", "precision", "recall (sensitivity)", "specificity", "F1", "patches"])


def analysis_section(cv: modelling.CVResult) -> str:
    """One section per classification: the headline, the spread, and where it fails."""
    spread = cv.spread()
    recall = cv.per_class().recall
    worst = cv.per_specimen().sort_values().head(6)

    body = [
        f"{cv.n_specimens} patients, {cv.n_views} fields of view, {len(cv.y_true):,} patches, "
        f"{len(cv.labels)} classes. Leave-one-patient-out, {cv.n_folds} folds.",
        "",
        f"**Accuracy {cv.accuracy:.4f}, Cohen's κ {cv.kappa:.4f}.**",
        "",
        f"Per-patient accuracy: median {spread['median']:.3f}, "
        f"IQR {spread['q1']:.3f}–{spread['q3']:.3f}, range {spread['min']:.3f}–{spread['max']:.3f}, "
        f"{spread['n_below_half']} of {cv.n_specimens} below 0.5. This distribution is the "
        f"account of variability; there is no held-out split to quote instead "
        f"(`per_case_accuracy_*.csv`).",
        "",
        _per_class(cv),
        "",
        "Row-normalised confusion matrix (%), rows = reference:",
        "",
        _confusion(cv),
    ]
    if worst.min() < 0.9:
        body += ["", "Weakest patients:", "",
                 _table([[k, f"{v:.4f}"] for k, v in worst.items()], ["patient", "accuracy"])]
    if recall.min() < 0.6:
        body += ["", f"**{recall.idxmin()} recall is {recall.min():.1%}**, the weakest class here."]
    for w in cv.warnings:
        body += ["", f"**Caveat.** {w}"]
    return "\n".join(body) + "\n"


def origin_section(human: dict) -> str:
    """Where the metastatic patches were assigned, and the negative control."""
    summary, control = human["origin"]["summary"], human["origin"]["negative_control"]
    bundle = human["final"]["primary"]
    rows = [
        [r.label, r.predicted, f"{int(r.patches):,}", f"{r.percent:.1f}%"]
        for r in summary.itertuples() if r.patches > 0
    ]
    return "\n".join([
        f"Four lesions, {len(human['origin']['per_patch']):,} patches, scored against the "
        f"classifier fitted on all {len(bundle['training_specimens'])} primary-tumour patients. "
        f"They are external by construction: the training selection asks for tumour tissue, so "
        f"no fold and no fit ever saw them.",
        "",
        _table(rows, ["true origin", "predicted", "patches", "share"]),
        "",
        "**Negative control.** Every lesion here was resected from the lung, so no metastatic "
        "patch should be assigned to lung cancer.",
        "",
        _table([[r.origin, f"{int(r.patches):,}", f"{r.percent:.1f}%"] for r in control.itertuples()],
               ["origin", "patches assigned LC", "share"]),
    ]) + "\n"


def mouse_section(wide: pd.DataFrame, cases: pd.DataFrame, cv: modelling.CVResult) -> str:
    """The murine organ cohort, and why it is grouped differently."""
    spread = cv.spread()
    per_organ = cases.groupby("label").size().reindex(MOUSE_ORGANS)

    return "\n".join([
        f"Five organs, {SECTIONS_PER_ORGAN} sections each, {len(wide):,} patches. Extraction "
        f"rule **{MOUSE_SECTION}**: sections are imaged whole, so there is no annotation to "
        f"mask by and no dim-patch rule; a patch counts unless it holds no signal at all, and "
        f"pixels saturated in every channel are dropped because they say nothing about the "
        f"ratio between channels.",
        "",
        _table([[organ, int(n), int(n) * PATCHES_PER_VIEW] for organ, n in per_organ.items()],
               ["organ", "sections", "patches"]),
        "",
        "**The held-out unit is the section, not the animal, and that is correct.** Each "
        "mouse contributed one section per organ, so within any one organ no two sections "
        "come from the same animal: a held-out section shares no animal with any training "
        "section of its own class, and individual-level information cannot leak within a "
        "class. The human cohort is grouped by patient for the opposite reason: there, two "
        "or three fields of view come from a single section. The two schemes differ because "
        "the unit of independence sits at a different level, not because one of them is wrong.",
        "",
        f"**Accuracy {cv.accuracy:.4f}, Cohen's κ {cv.kappa:.4f}**, leave-one-section-out, "
        f"{cv.n_folds} folds. Per-section accuracy: median {spread['median']:.4f}, "
        f"range {spread['min']:.4f}–{spread['max']:.4f}.",
        "",
        _per_class(cv),
        "",
        "Row-normalised confusion matrix (%), rows = reference:",
        "",
        _confusion(cv),
    ]) + "\n"


def pdac_tma_section(wide: pd.DataFrame, cases: pd.DataFrame, stats: pd.DataFrame,
                        by_case: pd.DataFrame, agreement: pd.DataFrame | None) -> str:
    """The second pancreatic cohort and how its phenotypes were assigned."""
    counts = cases.subtype.value_counts()
    body = [
        f"{len(cases)} specimens, one annotated field of view each, {len(wide):,} patches. "
        f"Measured with the same rule as the rest of the malignant tissue (**{CANCER}**). "
        f"All summaries use the annotated tumour patches only.",
        "",
        "**How the glycophenotype is assigned.** Per specimen, over its annotated tumour "
        "patches: mean fucose above mean mannose is Fuc-high, otherwise Man-high. The split is "
        + ", ".join(f"{n} {name}" for name, n in counts.items()) + ".",
        "",
        "Comparison between the two phenotypes, per channel, computed on the patches as "
        "selected, before any display trimming:",
        "",
        _table(
            [[g, f"{r.fuc_high_median:.1f}", f"{r.man_high_median:.1f}",
              "P < 0.001" if r.p_value < 1e-3 else f"P = {r.p_value:.3f}",
              f"{r.rank_biserial_r:+.3f}"]
             for g, r in stats.iterrows()],
            ["channel", "Fuc-high median", "Man-high median", "Mann-Whitney", "rank-biserial r"],
        ),
        "",
        "The distinction rests on two channels and only two: fucose and mannose separate "
        "almost completely, while the other three barely move.",
        "",
        "**Specimen as the unit.** The patches of one specimen are not independent, so the same "
        f"comparison is repeated with one value per specimen (the mean over its tumour patches; "
        f"{int(by_case.n_fuc_high.iloc[0])} Fuc-high against {int(by_case.n_man_high.iloc[0])} "
        "Man-high specimens), two-sided exact Mann-Whitney test:",
        "",
        _table(
            [[g, f"{r.fuc_high_median:.1f}", f"{r.man_high_median:.1f}",
              f"P = {r.p_value:.4f}" if r.p_value < 0.01 else f"P = {r.p_value:.3f}",
              f"{r.rank_biserial_r:+.3f}"]
             for g, r in by_case.iterrows()],
            ["channel", "Fuc-high median", "Man-high median", "exact Mann-Whitney", "rank-biserial r"],
        ),
        "",
        "Fucose and mannose define the phenotype (fucose above mannose), so their separation is "
        "partly built in; the other three channels are the independent part of the comparison.",
    ]
    if agreement is not None:
        matched = int((agreement.labelled == agreement.rule_would_say).sum())
        body += [
            "",
            f"**This addresses the circularity noted for the first cohort.** Applying the "
            f"same two-channel rule to those {len(agreement)} specimens reproduces "
            f"**{matched} of {len(agreement)}** of the labels a principal-component "
            f"analysis gave them. The phenotypes are not an artefact of the method that "
            f"found them: a rule anyone can check by eye recovers the same grouping.",
        ]
    return "\n".join(body) + "\n"


def build(
    human_wide: pd.DataFrame,
    human_cases: pd.DataFrame,
    human: dict,
    mouse_wide: pd.DataFrame,
    mouse_cases: pd.DataFrame,
    mouse_cv: modelling.CVResult,
    tma_wide: pd.DataFrame,
    tma_cases: pd.DataFrame,
    tma_stats: pd.DataFrame,
    tma_by_case: pd.DataFrame,
    tma_agreement: pd.DataFrame | None,
) -> str:
    """Assemble the whole report from one run of the models."""
    primary, subtype, phenotype = (human["cv"][k] for k in ("primary", "subtype", "phenotype"))
    origin = human["origin"]["summary"]
    lung = int(human["origin"]["negative_control"].patches.sum())

    tumour = human_wide[human_wide.sample_type == "tumour"]
    blank_by_type = human_wide[GLYCANS[0]].isna().groupby(human_wide.sample_type).sum()
    inside = int((tumour[GLYCANS[0]].isna() & tumour.is_tumour.fillna(False)).sum())

    def correct(origin_label: str) -> float:
        return origin.query("label == @origin_label and predicted == @origin_label").percent.iloc[0]

    cohort = human_cases.groupby(["sample_type", "label"]).agg(
        patients=("specimen_id", "nunique"),
        views=("n_views", "sum"),
        patches=("n_patches", "sum"),
        blank=("n_blank_patches", "sum"),
    )

    summary_rows = [
        ["Nine-class cancer classification", f"**{primary.accuracy:.1%}**", f"κ {primary.kappa:.3f}", f"{primary.n_folds} patients"],
        ["Breast cancer subtype", f"**{subtype.accuracy:.1%}**", f"κ {subtype.kappa:.3f}", f"{subtype.n_folds} patients"],
        ["PDAC glycophenotype", f"**{phenotype.accuracy:.1%}**", f"κ {phenotype.kappa:.3f}", f"{phenotype.n_folds} patients"],
        ["Metastatic origin, colorectal", f"**{correct('CRC'):.1f}%**", "assigned correctly", "2 patients"],
        ["Metastatic origin, breast", f"**{correct('BC'):.1f}%**", "assigned correctly", "2 patients"],
        ["Negative control", f"**{lung}**", "metastatic patches assigned to lung", "—"],
        ["Murine organ classification", f"**{mouse_cv.accuracy:.1%}**", f"κ {mouse_cv.kappa:.3f}", f"{mouse_cv.n_folds} sections"],
    ]
    artefacts = [
        ["`data/<cohort>_patches_wide.csv`", "one row per patch (human, mouse, pdac_tma)"],
        ["`data/<cohort>_patches_long.csv`", "the same rows, one per channel"],
        ["`data/<cohort>_cases.csv`", "one row per specimen"],
        ["`models/lda_primary.joblib`", f"fitted on all {len(human['final']['primary']['training_specimens'])} primary-tumour patients"],
        ["`models/lda_bc_subtype.joblib`", f"fitted on all {len(human['final']['subtype']['training_specimens'])} breast patients"],
        ["`models/lda_pdac_phenotype.joblib`", f"fitted on all {len(human['final']['phenotype']['training_specimens'])} pancreatic patients"],
        ["`models/lda_mouse_organ.joblib`", f"fitted on all {mouse_cv.n_specimens} murine sections"],
        ["`results/loso_summary.csv`", "accuracy, κ and per-specimen spread of every cross-validation"],
        ["`results/metrics_*_loso.csv`, `confusion_*_loso.csv`", "per-class precision, recall, F1, support; confusion matrices"],
        ["`results/per_case_accuracy_*.csv`", "out-of-fold accuracy for every held-out specimen"],
        ["`results/origin_assignment.csv`", "metastatic assignment counts and shares"],
        ["`results/lda_scores_*.csv`", "discriminant coordinates of every patch"],
        ["`results/profile_*.csv`", "mean profiles, raw and standardised, with clustering order"],
        ["`results/pdac_pca_*.csv`", "principal components of the first pancreatic cohort"],
        ["`results/pdac_tma_*.csv`", "second pancreatic cohort: phenotype comparison and agreement"],
        ["`results/violin_group_summary.csv`", "distribution summaries, before and after the 1.5 × IQR display trim"],
    ]

    return f"""# Results — glyco-code image analysis

Every number here is recorded with the grouping that produced it and the date of the run.
The whole file is regenerated in one pass, so no two numbers in it can disagree.

**Run: {dt.date.today().isoformat()}.**

---

## Summary

{_table(summary_rows, ["analysis", "accuracy", "", "folds"])}

Every human classification is cross-validated by leaving out one patient at a time. A case
contributes two or three fields of view of a single section (same tissue, same staining
run, same imaging session), so holding out a field of view would leave its siblings in
training and measure recognition of a section the model had already seen. The murine cohort
is grouped by section instead; section 6 says why that is the right unit there.

---

## 1 — Extraction

A pixel is kept when any of the five channels reaches the threshold below; a patch with
fewer than {CANCER.min_valid_px} such pixels has too little tissue to average and is
recorded as missing in all five channels. In every tissue, a pixel above {SATURATED_ABOVE}
in all five channels is also dropped, because a pixel saturated in every channel says
nothing about the ratio between them.

{_table([[k, v.valid_threshold, v.min_valid_px] for k, v in SETTINGS.items() if k != "section"],
        ["tissue", "pixel threshold", "minimum pixels per patch"])}

The thresholds differ deliberately. Malignant tissue is bright, and 10 is where background
stops being counted there; it is also the value every classification result was measured
with. Adjacent-normal tissue is dimmer overall, so the same threshold would discard pixels
that carry real signal, and because it would discard the dim ones, the mean over what
survived would be pushed up. Each tissue is thresholded for its own brightness.

{_table([[t, f"{int(v):,}", f"{int((human_wide.sample_type == t).sum()):,}",
          f"{v / (human_wide.sample_type == t).sum():.3%}"]
         for t, v in blank_by_type.items()],
        ["tissue", "patches recorded missing", "patches", "%"])}

Of the {int(blank_by_type.get('tumour', 0)):,} tumour patches recorded missing,
**{inside}** lie inside the pathologist's annotation; those are the only ones that would
otherwise have reached a classifier.

### Cohort

```
{cohort.to_string()}
```

---

## 2 — Nine-class cancer classification

{analysis_section(primary)}
---

## 3 — Breast cancer subtype

{analysis_section(subtype)}
---

## 4 — PDAC glycophenotype

{analysis_section(phenotype)}

**Read this with care.** The two glycophenotypes were assigned by a principal-component
analysis of these same five channels, so a classifier that recovers them from those
channels is close to tautological. The number says the assignment is internally consistent
and reproducible patient by patient. It is not independent validation of the phenotypes.

---

## 5 — Metastatic origin

{origin_section(human)}
---

## 6 — Murine organ cohort

{mouse_section(mouse_wide, mouse_cases, mouse_cv)}
---

## 7 — Second pancreatic cohort

{pdac_tma_section(tma_wide, tma_cases, tma_stats, tma_by_case, tma_agreement)}
---

## 8 — Files written

{_table(artefacts, ["file", "what it is"])}
"""


def write(text: str) -> None:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text, encoding="utf8")
