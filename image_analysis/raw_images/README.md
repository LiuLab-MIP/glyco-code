# raw_images

The raw fluorescence images (about 3 GB) are **not tracked by git**: everything in this
folder except this file is ignored (see `.gitignore` at the root of the project).

## Download and unpack

The images are deposited on Zenodo: [10.5281/zenodo.23316047](https://doi.org/10.5281/zenodo.23316047), as 19 zip
files of at most 183 MB each (large single files are difficult to upload and download reliably):

| Files | Size | Content | Needed for |
|---|---|---|---|
| `glyco-code_raw_images_grey_01of05.zip` to `_05of05.zip` (5 files) | 788 MB | `grey/`: 1,124 files, the images the code reads | `python -m glyco extract`, the image tests |
| `glyco-code_raw_images_colour_01of14.zip` to `_14of14.zip` (14 files) | 2.3 GB | `colour/`: 1,530 files, pseudo-colour, brightfield, merge and H&E for viewing | looking at the specimens only |

The files are independent zip archives, not pieces of one split archive: each holds whole
specimens, every file is in exactly one of them, and all of them have the same top-level
folder (`grey/` or `colour/`). Unzip **all parts of a set into the same folder**,
`image_analysis/raw_images/`; the folders merge:

```bash
cd image_analysis/raw_images
# after downloading the files into this folder
for f in glyco-code_raw_images_grey_*.zip;   do unzip -q "$f"; done
for f in glyco-code_raw_images_colour_*.zip; do unzip -q "$f"; done   # optional
rm glyco-code_raw_images_*.zip         # the archives are not needed afterwards
```

On Windows, select all the downloaded `.zip` files of a set in File Explorer, right-click and
choose "Extract All" (or extract each one into the same folder). If a download fails, only
that part has to be fetched again. `python -m pytest tests/test_extraction.py` checks the
unpacked tree.

The code needs only `grey/`; the colour images are for people. To keep the images elsewhere,
point `GLYCO_RAW_IMAGES` at the folder that contains `grey/` and `colour/`.

```
raw_images/
  grey/      read by the code: one 8-bit grey image per channel, plus the annotation
  colour/    for people: pseudo-coloured channels, brightfield, merge, H&E
```

The two trees share one directory structure:

```
<tree>/<cohort>/<sample_type>/<label>/<specimen_id>/<view>/<file>.jpg
```

| Level | Values |
|---|---|
| `cohort` | `human`, `mouse`, `pdac_batch2` |
| `sample_type` | `tumour`, `metastasis`, `normal` (human); `section` (mouse); `tumour` (pdac_batch2) |
| `label` | cancer type (`HCC`, `CRC`, `LC`, `GC`, `CCA`, `EC`, `PCa`, `BC`, `PDAC`) or organ |
| `specimen_id` | pseudonym such as `HCC-P01`, `CRC-LM-P01`, `CRC-N-P01`, `PDAC2-P01`, `heart-1` |
| `view` | `1`, `2`, ... (one field of view of the specimen) |

Files: `grey/` has `GlcNAc.jpg Gal.jpg SA.jpg Fuc.jpg Man.jpg` and, for annotated specimens,
`annotation.jpg`; `colour/` has the same five names (pseudo-coloured) plus `brightfield.jpg`,
`merge.jpg` and, where available, `he.jpg`. Totals: 164 human fields of view, 15 mouse
sections, 15 pancreatic specimens of the second cohort.

Contents are described in detail in `../../docs/DATA_DICTIONARY.md`. To check an unpacked
copy, run `python -m pytest tests/test_extraction.py`.
