"""The three cohorts, where their tables live, and what each should contain.

The raw-image tree encodes everything the pipeline needs to know about a specimen:

    <GREY_IMAGES>/<cohort>/<sample_type>/<label>/<specimen_id>/<view>/<glycan>.jpg
                                                                    /annotation.jpg   (where annotated)

so the code never needs to know a folder name. The only per-specimen information that is
not in the path is the sub-classification (`subtype`), which is read from
`metadata/specimens.csv`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import DATA_DIR, GREY_IMAGES

#: Exactly the kinds of specimen there are.
SAMPLE_TYPES = ("tumour", "metastasis", "normal", "section")


@dataclass(frozen=True)
class Cohort:
    """One experiment, kept in its own tables throughout."""

    name: str
    #: Prefix of the table file names, e.g. `human_patches_wide.csv`.
    table_prefix: str
    #: Sample types present, in the order they are walked.
    sample_types: tuple[str, ...]
    #: Whether specimens of this cohort carry a pathologist's annotation image.
    annotated_sample_types: tuple[str, ...]
    #: Order in which labels are walked within a sample type. Fixes the row order of the
    #: tables, which makes every output reproducible to the last digit.
    label_order: tuple[str, ...]
    expected_specimens: dict[str, int]
    expected_views: dict[str, int]

    @property
    def image_root(self) -> Path:
        return GREY_IMAGES / self.name

    def table(self, kind: str) -> Path:
        """Path of `wide`, `long` or `cases` for this cohort."""
        stem = {"wide": "patches_wide", "long": "patches_long", "cases": "cases"}[kind]
        return DATA_DIR / f"{self.table_prefix}_{stem}.csv"


# --------------------------------------------------------------------------------------
# Human tissue
# --------------------------------------------------------------------------------------

#: Cancer-type codes: hepatocellular carcinoma, colorectal cancer, lung cancer, gastric
#: cancer, cholangiocarcinoma, oesophageal cancer, prostate cancer, breast cancer,
#: pancreatic ductal adenocarcinoma.
CANCER_TYPES = ("HCC", "CRC", "LC", "GC", "CCA", "EC", "PCa", "BC", "PDAC")

#: A metastasis is labelled with its tissue of ORIGIN, not the organ the lesion sits in; an
#: adjacent-normal block with the cancer type it is adjacent to. Membership of an analysis
#: is therefore decided by `sample_type`, never by `label`.
HUMAN = Cohort(
    name="human",
    table_prefix="human",
    sample_types=("tumour", "metastasis", "normal"),
    annotated_sample_types=("tumour", "metastasis"),
    label_order=CANCER_TYPES,
    expected_specimens={"tumour": 72, "metastasis": 4, "normal": 25},
    expected_views={"tumour": 131, "metastasis": 8, "normal": 25},
)

#: label -> (patients, views), primary tumours only.
EXPECTED_PRIMARY = {
    "CRC": (8, 14), "HCC": (8, 17), "GC": (7, 11), "EC": (8, 14), "PDAC": (7, 10),
    "BC": (13, 29), "LC": (7, 11), "CCA": (8, 13), "PCa": (6, 12),
}

#: Specimens carrying a sub-classification. Breast values are clinical (IHC) subtypes taken
#: from the pathology report; pancreatic values are the glycophenotypes this study assigned
#: by principal-component analysis of the same five channels.
EXPECTED_SUBTYPES = {"BC": 13, "PDAC": 7}
BREAST_SUBTYPES = ("HER2", "LuminalA", "LuminalB", "TNBC")
PDAC_PHENOTYPES = ("TypeI", "TypeII")

#: Every metastatic specimen was resected from the lung, so `LC` is the negative control of
#: the origin prediction: no metastatic patch should be assigned to it.
NEGATIVE_CONTROL = "LC"

# --------------------------------------------------------------------------------------
# Mouse organ sections
# --------------------------------------------------------------------------------------

#: Five organs, three sections each, one section per organ per animal. That last point is
#: what makes the section the right unit to hold out: within a class no two sections come
#: from the same animal, so holding one out cannot leak individual-level information into
#: the rest of its class. The human cohort is grouped by patient for the opposite reason:
#: there, several fields of view do share a section.
MOUSE_ORGANS = ("heart", "kidney", "liver", "lung", "spleen")
SECTIONS_PER_ORGAN = 3

MOUSE = Cohort(
    name="mouse",
    table_prefix="mouse",
    sample_types=("section",),
    annotated_sample_types=(),
    label_order=MOUSE_ORGANS,
    expected_specimens={"section": len(MOUSE_ORGANS) * SECTIONS_PER_ORGAN},
    expected_views={"section": len(MOUSE_ORGANS) * SECTIONS_PER_ORGAN},
)

# --------------------------------------------------------------------------------------
# Second pancreatic cohort
# --------------------------------------------------------------------------------------

#: A tissue microarray of annotated pancreatic cases, one field of view each. Acquired apart
#: from the first pancreatic cohort and kept in separate tables, because its glycophenotypes are
#: assigned by a different rule (whichever of fucose and mannose is stronger).
PDAC_TMA = Cohort(
    name="pdac_tma",
    table_prefix="pdac_tma",
    sample_types=("tumour",),
    annotated_sample_types=("tumour",),
    label_order=("PDAC",),
    expected_specimens={"tumour": 15},
    expected_views={"tumour": 15},
)
PHENOTYPE_FUC_HIGH = "FucHigh"
PHENOTYPE_MAN_HIGH = "ManHigh"

COHORTS = {c.name: c for c in (HUMAN, MOUSE, PDAC_TMA)}
