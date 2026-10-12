"""Locations, image geometry and the extraction rules.

Every tunable that affects a number lives here, so a change that moves a result shows up in
a single diff.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------------------
# Locations
# --------------------------------------------------------------------------------------

#: The `image_analysis/` directory (the parent of this package).
PACKAGE_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PACKAGE_ROOT / "data"
MODEL_DIR = PACKAGE_ROOT / "models"
RESULTS_DIR = PACKAGE_ROOT / "results"
SPECIMEN_MANIFEST = PACKAGE_ROOT / "metadata" / "specimens.csv"

#: Raw images are not tracked by git. Point `GLYCO_RAW_IMAGES` at a directory holding the
#: `grey/` and `colour/` trees to read them from somewhere else.
RAW_IMAGES = Path(os.environ.get("GLYCO_RAW_IMAGES", PACKAGE_ROOT / "raw_images"))
GREY_IMAGES = RAW_IMAGES / "grey"
COLOUR_IMAGES = RAW_IMAGES / "colour"

# --------------------------------------------------------------------------------------
# Channels
# --------------------------------------------------------------------------------------

#: Terminal monosaccharides. The order is the channel order of the acquisition and the
#: column order of every table. Each grey image is named after its channel: `<glycan>.jpg`.
GLYCANS = ["GlcNAc", "Gal", "SA", "Fuc", "Man"]

# --------------------------------------------------------------------------------------
# Image geometry
# --------------------------------------------------------------------------------------

IMAGE_SIZE = 1024
GRID = 16
BLOCK = IMAGE_SIZE // GRID  # 64 x 64 px per patch, 256 patches per field of view
PATCHES_PER_VIEW = GRID * GRID

#: The annotation image is 264 x 264: a 256 x 256 grid of 16 x 16 px blocks inside an 8 px
#: frame. Sampling one pixel per block at `FILTER_OFFSET + FILTER_BLOCK * i` lands on the
#: block centre; the offset carries the frame width.
ANNOTATION_FILE = "annotation.jpg"
FILTER_BLOCK = 16
FILTER_OFFSET = 15

#: Annotation colours, with a per-channel tolerance. The tolerance is widest on blue because
#: that is the channel in which the two colours differ least (212 vs 227).
COLOUR_CANCER = (255, 255, 212)
COLOUR_NON_CANCER = (203, 203, 227)
COLOUR_TOLERANCE = (10, 10, 30)

#: ITU-R BT.601 luma weights, applied with round-half-away-from-zero.
LUMA_WEIGHTS = (0.298936021293776, 0.587043074451121, 0.114020904255103)

# --------------------------------------------------------------------------------------
# Extraction rules
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ExtractionSettings:
    """The pixel mask and the blank-patch rule.

    A pixel is kept when ANY channel reaches `valid_threshold` (equivalently, it is dropped
    only when all five channels are below it), so empty background cannot pull a patch mean
    down. A patch with too few surviving pixels carries no usable signal and is recorded as
    missing rather than as a small number.
    """

    name: str
    valid_threshold: int
    #: A pixel is also dropped when EVERY channel is above this: a pixel saturated in all
    #: five channels carries no information about their ratio.
    saturated_above: int | None
    #: Below this many surviving pixels the patch is blank: NaN in all five channels.
    min_valid_px: int

    def __str__(self) -> str:
        sat = "none" if self.saturated_above is None else f">{self.saturated_above} in all"
        return f"{self.name} (T={self.valid_threshold}, saturation {sat}, N={self.min_valid_px})"


#: One saturation rule for every tissue.
SATURATED_ABOVE = 250

#: The pixel threshold differs by tissue, deliberately. Malignant tissue is bright, and 10 is
#: the value at which background stops being counted there; it is also the value every
#: classification result was measured with. Adjacent-normal tissue is dimmer overall, so the
#: same threshold would discard pixels that carry real signal, and because it would discard
#: the dim ones the mean over the survivors would be pushed up. 5 keeps that signal.
CANCER = ExtractionSettings("cancer", valid_threshold=10, saturated_above=SATURATED_ABOVE, min_valid_px=512)
ADJACENT_NORMAL = ExtractionSettings(
    "adjacent-normal", valid_threshold=5, saturated_above=SATURATED_ABOVE, min_valid_px=512
)
#: Mouse sections are imaged whole, with no annotation and no dim-patch rule: a patch counts
#: unless it contains no signal at all.
MOUSE_SECTION = ExtractionSettings(
    "mouse section", valid_threshold=5, saturated_above=SATURATED_ABOVE, min_valid_px=1
)

#: Which rule each kind of specimen is measured with.
SETTINGS = {
    "tumour": CANCER,
    "metastasis": CANCER,
    "normal": ADJACENT_NORMAL,
    "section": MOUSE_SECTION,
}
