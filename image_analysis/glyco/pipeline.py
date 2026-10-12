"""The stages of the pipeline as functions, shared by the command line and the tests."""

from __future__ import annotations

import time

import pandas as pd

from . import extraction, modelling, mouse, pdac_tma, report, results, schema
from .cohorts import COHORTS, HUMAN, MOUSE, PDAC_TMA, Cohort

#: Cohort-specific step applied to freshly extracted tables before they are written.
POSTPROCESS = {PDAC_TMA.name: pdac_tma.annotate_phenotype}
CHECKS = {MOUSE.name: mouse.check, PDAC_TMA.name: pdac_tma.check}


def extract_cohort(cohort: Cohort) -> int:
    """Read every image of a cohort once and write its three tables. Returns an exit code."""
    print(f"\n{'=' * 78}\nextraction: {cohort.name}\n{'=' * 78}")
    wide, cases = extraction.build(cohort)
    if cohort.name in POSTPROCESS:
        wide, cases = POSTPROCESS[cohort.name](wide, cases)
    long = schema.write_tables(cohort, wide, cases)
    extraction.check_round_trip(cohort, wide, long)
    extraction.describe(wide, cases)

    problems = extraction.check_cohort(cohort, wide, cases)
    if cohort.name in CHECKS:
        problems += CHECKS[cohort.name](wide, cases)
    if problems:
        print("\n*** COHORT MISMATCH, stopping ***")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print(f"\n{cohort.name}: matches the expected cohort. "
          f"{len(wide):,} patches, {len(cases)} specimens.")
    return 0


def extract(names: list[str] | None = None) -> int:
    for name in names or list(COHORTS):
        if extract_cohort(COHORTS[name]) != 0:
            return 1
    return 0


def analyse() -> int:
    """Cross-validate, fit and serialise the classifiers, write the result tables and report."""
    print(f"\n{'=' * 78}\nanalysis\n{'=' * 78}")
    human_wide, human_long = schema.read_wide(HUMAN), schema.read_long(HUMAN)
    human_cases = schema.read_cases(HUMAN)
    mouse_wide, mouse_cases = mouse.read_wide(), mouse.read_cases()
    tma_wide, tma_cases = pdac_tma.read_wide(), pdac_tma.read_cases()

    human = modelling.run_human(human_wide)
    mouse_cv = mouse.cross_validate(mouse_wide)
    mouse.fit_final(mouse_wide)

    for name, cv in {**human["cv"], mouse.ANALYSIS: mouse_cv}.items():
        spread = cv.spread()
        print(f"\n{modelling.TITLES[name]}: {cv.n_folds} held-out specimens, "
              f"{len(cv.y_true):,} patches, {len(cv.labels)} classes")
        print(f"  accuracy {cv.accuracy:.4f}   kappa {cv.kappa:.4f}")
        print(f"  per-specimen accuracy: median {spread['median']:.3f}, "
              f"IQR {spread['q1']:.3f}-{spread['q3']:.3f}, min {spread['min']:.3f}, "
              f"{spread['n_below_half']} below 0.5")
        for warning in cv.warnings:
            print(f"  ! {warning}")
    print("\nmetastatic origin:")
    print(human["origin"]["summary"].query("patches > 0").to_string(index=False))
    print("\nnegative control:")
    print(human["origin"]["negative_control"].to_string(index=False))

    tables = results.build_all(human_wide, human_long, human, mouse_wide, mouse_cv, tma_wide)
    results.write_all(tables)
    text = report.build(
        human_wide, human_cases, human, mouse_wide, mouse_cases, mouse_cv,
        tma_wide, tma_cases,
        tables["pdac_tma_phenotype_stats"].set_index("glycan"),
        tables["pdac_tma_phenotype_stats_by_case"].set_index("glycan"),
        pdac_tma.compare_with_first_cohort(),
    )
    report.write(text)
    print(f"\nwrote {len(tables)} result tables and {report.REPORT.name}")
    return 0


def run_all(*, skip_extract: bool = False) -> int:
    began = time.time()
    if not skip_extract and extract() != 0:
        return 1
    status = analyse()
    print(f"\ndone in {time.time() - began:.0f}s")
    return status
