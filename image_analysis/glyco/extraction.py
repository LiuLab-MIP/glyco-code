"""Extraction: raw images to patch tables.

Walks every specimen of a cohort once; cuts each 1024 x 1024 field of view into a 16 x 16
grid of patches; records, per patch, how many pixels carried signal and the mean intensity
of each of the five channels over those pixels; reads the pathologist's annotation where
one exists; and returns three tables.

    wide    one row per patch, five intensity columns
    long    the same rows reshaped to one row per patch per channel (derived, never measured)
    cases   one row per specimen, for reading rather than joining

The images are walked once and only once. Everything else is derived from the wide table,
so the summaries and the classifier cannot drift apart.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from PIL import Image

from . import schema
from .cohorts import (
    EXPECTED_PRIMARY,
    EXPECTED_SUBTYPES,
    HUMAN,
    Cohort,
)
from .config import (
    ANNOTATION_FILE,
    BLOCK,
    COLOUR_CANCER,
    COLOUR_NON_CANCER,
    COLOUR_TOLERANCE,
    FILTER_BLOCK,
    FILTER_OFFSET,
    GLYCANS,
    GRID,
    IMAGE_SIZE,
    LUMA_WEIGHTS,
    PATCHES_PER_VIEW,
    SETTINGS,
    SPECIMEN_MANIFEST,
    ExtractionSettings,
)

# --------------------------------------------------------------------------------------
# One field of view
# --------------------------------------------------------------------------------------


def to_greyscale(path: Path) -> np.ndarray:
    """Read one channel image as an 8-bit greyscale array.

    Colour images are reduced by luma weights and rounded half away from zero; an image
    that is already single-channel passes through untouched.
    """
    arr = np.asarray(Image.open(path))
    if arr.ndim == 2:
        return arr
    r, g, b = LUMA_WEIGHTS
    luma = (
        arr[:, :, 0].astype(np.float64) * r
        + arr[:, :, 1].astype(np.float64) * g
        + arr[:, :, 2].astype(np.float64) * b
    )
    return np.floor(luma + 0.5).clip(0, 255).astype(np.uint8)


def _as_patches(arr: np.ndarray) -> np.ndarray:
    """Reshape (..., 1024, 1024) to (..., 16, 16, 64, 64), indexed [row, col, y, x]."""
    lead = arr.shape[:-2]
    return arr.reshape(*lead, GRID, BLOCK, GRID, BLOCK).swapaxes(-3, -2)


def patch_statistics(
    view_dir: Path, settings: ExtractionSettings
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Signal-pixel count and per-channel mean for every patch of one field of view.

    A pixel counts as signal when ANY of the five channels reaches the threshold, so empty
    background cannot drag a patch mean down. Where the rule asks for it, a pixel saturated
    in every channel is dropped as well. Means are taken over the signal pixels only.

    A patch with too few signal pixels has nothing to average and is recorded as missing in
    all five channels rather than as a small number. Returns a (16, 16) count and one
    (16, 16) array of means per channel.
    """
    stack = np.stack([to_greyscale(view_dir / f"{glycan}.jpg") for glycan in GLYCANS])
    if stack.shape[1:] != (IMAGE_SIZE, IMAGE_SIZE):
        raise ValueError(
            f"every channel must be {IMAGE_SIZE}x{IMAGE_SIZE}; got {stack.shape[1:]} in {view_dir}"
        )

    signal = (stack >= settings.valid_threshold).any(axis=0)
    if settings.saturated_above is not None:
        signal &= ~(stack > settings.saturated_above).all(axis=0)

    signal_patches = _as_patches(signal)
    n_valid = signal_patches.sum(axis=(-2, -1)).astype(np.int32)

    totals = (_as_patches(stack.astype(np.float64)) * signal_patches).sum(axis=(-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        means = totals / n_valid
    means[:, n_valid < settings.min_valid_px] = np.nan

    return n_valid, dict(zip(GLYCANS, means))


def _matches(pixel: np.ndarray, colour: tuple[int, int, int]) -> bool:
    """Is this pixel one of the two annotation colours, within tolerance?"""
    return bool(
        np.all(
            np.abs(pixel.astype(np.int16) - np.array(colour, np.int16))
            <= np.array(COLOUR_TOLERANCE, np.int16)
        )
    )


def read_annotation(path: Path) -> np.ndarray:
    """(16, 16) bool: which patches the pathologist marked as tumour.

    The annotation image carries an 8 px frame around its 256 x 256 grid, so one pixel is
    sampled per block at the block centre. Any colour that is neither of the two expected
    ones stops the run, naming exactly where it was found.
    """
    arr = np.asarray(Image.open(path).convert("RGB"))
    out = np.zeros((GRID, GRID), dtype=bool)
    for row in range(GRID):
        for col in range(GRID):
            y, x = FILTER_OFFSET + FILTER_BLOCK * row, FILTER_OFFSET + FILTER_BLOCK * col
            pixel = arr[y, x, :]
            is_cancer = _matches(pixel, COLOUR_CANCER)
            is_non_cancer = _matches(pixel, COLOUR_NON_CANCER)
            if is_cancer and not is_non_cancer:
                out[row, col] = True
            elif is_non_cancer:
                out[row, col] = False
            else:
                raise ValueError(
                    "Unrecognised colour in annotation\n"
                    f"  file         : {path}\n"
                    f"  block        : row {row + 1}, column {col + 1} of {GRID}\n"
                    f"  pixel (y, x) : ({y}, {x})\n"
                    f"  RGB          : {tuple(int(v) for v in pixel)}\n"
                    f"  expected     : {COLOUR_CANCER} (tumour) or {COLOUR_NON_CANCER} "
                    f"(non-tumour), within +/-{COLOUR_TOLERANCE}"
                )
    return out


# --------------------------------------------------------------------------------------
# Walking the image tree
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ViewRecord:
    """One field of view found in the tree, with everything its path says about it."""

    specimen_id: str
    label: str
    sample_type: str
    view: str
    directory: Path


def _children(path: Path) -> list[Path]:
    """Sorted child directories, ignoring dot-prefixed ones (editor and tooling folders)."""
    return sorted((p for p in path.iterdir() if p.is_dir() and not p.name.startswith(".")),
                  key=lambda p: p.name)


def _natural(path: Path) -> tuple:
    """Sort key putting numeric names in numeric order (`2` before `10`)."""
    return (0, int(path.name), "") if path.name.isdigit() else (1, 0, path.name)


def read_manifest() -> pd.DataFrame:
    """`metadata/specimens.csv`: the sub-classification of every specimen."""
    return pd.read_csv(SPECIMEN_MANIFEST, dtype="string")


def iter_views(cohort: Cohort) -> Iterator[ViewRecord]:
    """Every field of view of a cohort, in a fixed order: sample type, label, specimen, view."""
    root = cohort.image_root
    if not root.is_dir():
        raise FileNotFoundError(
            f"missing image folder {root}\n"
            "Put the raw images there or point GLYCO_RAW_IMAGES at the folder holding grey/."
        )
    for sample_type in cohort.sample_types:
        type_dir = root / sample_type
        labels = {p.name: p for p in _children(type_dir)}
        unknown = set(labels) - set(cohort.label_order)
        if unknown:
            raise ValueError(f"unexpected label folders under {type_dir}: {sorted(unknown)}")
        for label in cohort.label_order:
            if label not in labels:
                continue
            for specimen in _children(labels[label]):
                for view in sorted(_children(specimen), key=_natural):
                    yield ViewRecord(specimen.name, label, sample_type, view.name, view)


def build(cohort: Cohort = HUMAN) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Traverse every image of a cohort once and return (patch table, case summary)."""
    manifest = read_manifest().set_index("specimen_id")
    rows: list[dict] = []
    cases: dict[str, dict] = {}

    for record in iter_views(cohort):
        try:
            entry = manifest.loc[record.specimen_id]
        except KeyError:
            raise ValueError(f"{record.specimen_id} is in the image tree but not in {SPECIMEN_MANIFEST}")
        if (entry.cohort, entry.sample_type, entry.label) != (cohort.name, record.sample_type, record.label):
            raise ValueError(f"{record.specimen_id}: the image path disagrees with the manifest")

        rule = SETTINGS[record.sample_type]
        has_annotation = record.sample_type in cohort.annotated_sample_types
        subtype = None if pd.isna(entry.subtype) else str(entry.subtype)
        view_id = f"{record.specimen_id}-{record.view}"

        n_valid, means = patch_statistics(record.directory, rule)
        annotation = read_annotation(record.directory / ANNOTATION_FILE) if has_annotation else None

        rows.extend(
            {
                "specimen_id": record.specimen_id,
                "label": record.label,
                "sample_type": record.sample_type,
                "subtype": subtype,
                "view_id": view_id,
                "grid_row": row + 1,
                "grid_col": col + 1,
                "n_valid_px": int(n_valid[row, col]),
                "is_tumour": None if annotation is None else bool(annotation[row, col]),
                **{g: means[g][row, col] for g in GLYCANS},
            }
            for row in range(GRID)
            for col in range(GRID)
        )

        case = cases.setdefault(
            record.specimen_id,
            {
                "specimen_id": record.specimen_id,
                "cohort": cohort.name,
                "sample_type": record.sample_type,
                "label": record.label,
                "subtype": subtype,
                "n_views": 0,
                "n_patches": 0,
                "n_tumour_patches": 0 if has_annotation else None,
                "n_blank_patches": 0,
                "has_annotation": has_annotation,
            },
        )
        case["n_views"] += 1
        case["n_patches"] += PATCHES_PER_VIEW
        case["n_blank_patches"] += int((n_valid < rule.min_valid_px).sum())
        if annotation is not None:
            case["n_tumour_patches"] += int(annotation.sum())

    wide = pd.DataFrame(rows, columns=schema.WIDE_COLUMNS).astype(schema.DTYPES)
    case_table = pd.DataFrame(list(cases.values()), columns=schema.CASE_COLUMNS)
    return wide, case_table.astype(schema.CASE_DTYPES)


# --------------------------------------------------------------------------------------
# Checks that run on every build, not only under pytest
# --------------------------------------------------------------------------------------


def check_cohort(cohort: Cohort, wide: pd.DataFrame, cases: pd.DataFrame) -> list[str]:
    """Compare what was extracted against what the cohort should contain.

    An extraction that silently loses a specimen is worse than one that fails, so this runs
    every time and its problems stop the run.
    """
    problems: list[str] = []

    specimens = cases.groupby("sample_type").size().to_dict()
    views = wide.groupby("sample_type")["view_id"].nunique().to_dict()
    for sample_type in cohort.sample_types:
        got, want = specimens.get(sample_type, 0), cohort.expected_specimens[sample_type]
        if got != want:
            problems.append(f"{sample_type}: {got} specimens, expected {want}")
        got, want = views.get(sample_type, 0), cohort.expected_views[sample_type]
        if got != want:
            problems.append(f"{sample_type}: {got} views, expected {want}")

    expected_rows = sum(cohort.expected_views.values()) * PATCHES_PER_VIEW
    if len(wide) != expected_rows:
        problems.append(f"{len(wide):,} rows, expected {expected_rows:,}")
    per_view = wide.groupby("view_id").size()
    if not (per_view == PATCHES_PER_VIEW).all():
        problems.append(f"views without {PATCHES_PER_VIEW} rows: "
                        f"{per_view[per_view != PATCHES_PER_VIEW].to_dict()}")

    if cohort is HUMAN:
        tumour = wide[wide.sample_type == "tumour"]
        for label, (patients, view_count) in EXPECTED_PRIMARY.items():
            sub = tumour[tumour.label == label]
            if (sub.specimen_id.nunique(), sub.view_id.nunique()) != (patients, view_count):
                problems.append(
                    f"{label}: {sub.specimen_id.nunique()} patients / {sub.view_id.nunique()} views, "
                    f"expected {patients} / {view_count}"
                )
        for label, want in EXPECTED_SUBTYPES.items():
            got = tumour[tumour.label == label].specimen_id.nunique()
            if got != want:
                problems.append(f"{label}: {got} specimens carrying a subtype, expected {want}")

    return problems


def check_round_trip(cohort: Cohort, wide: pd.DataFrame, long: pd.DataFrame) -> None:
    """Read the tables back and assert the CSV lost nothing.

    Missing values are the thing at risk: an empty field must come back as missing, not as
    zero or as the string "". Dtypes and per-column missing counts must be identical.
    """
    back_wide, back_long = schema.read_wide(cohort), schema.read_long(cohort)
    for name, before, after in (("wide", wide, back_wide), ("long", long, back_long)):
        assert list(before.columns) == list(after.columns), f"{name}: column order changed"
        assert before.dtypes.equals(after.dtypes), (
            f"{name}: dtypes changed\n{pd.concat([before.dtypes, after.dtypes], axis=1)}"
        )
        assert before.isna().sum().equals(after.isna().sum()), f"{name}: missing counts changed"
        pd.testing.assert_frame_equal(before.reset_index(drop=True), after.reset_index(drop=True))
    pd.testing.assert_frame_equal(
        schema.to_long(back_wide).reset_index(drop=True), back_long.reset_index(drop=True)
    )


def describe(wide: pd.DataFrame, cases: pd.DataFrame) -> None:
    """Print the case table and the cohort breakdown, for reading before anything else."""
    print("\n" + "=" * 100 + "\ncases\n" + "=" * 100)
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(cases.to_string(index=False))

    summary = wide.groupby(["sample_type", "label"]).agg(
        specimens=("specimen_id", "nunique"),
        views=("view_id", "nunique"),
        patches=("view_id", "size"),
    )
    summary["blank"] = wide[GLYCANS[0]].isna().groupby([wide.sample_type, wide.label]).sum()
    print("\n" + summary.to_string())
    print(f"\nblank patches: {int(wide[GLYCANS[0]].isna().sum()):,} of {len(wide):,}")
