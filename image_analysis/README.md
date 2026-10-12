# image_analysis

Python package `glyco`: from fluorescence images of five terminal monosaccharides to patch
tables, classifiers and result tables. Installation: `../docs/INSTALL.md`; demos:
`../docs/DEMO.md`; table and file formats: `../docs/DATA_DICTIONARY.md`.

## Commands

Run from this directory, inside the virtual environment.

| Command | What it does | Needs |
|---|---|---|
| `python -m glyco extract [--cohort human\|mouse\|pdac_tma]` | reads the images once and writes `data/<cohort>_patches_wide.csv`, `_patches_long.csv` and `_cases.csv`; checks the cohort against the expected counts | `raw_images/grey/` |
| `python -m glyco analyse` | cross-validates, fits and serialises the classifiers (`models/`), writes the result tables and `results/RESULTS.md` | `data/` |
| `python -m glyco all [--skip-extract]` | extraction followed by analysis; with `--skip-extract` it uses the tables already in `data/` | |
| `python -m pytest` | 63 tests | `data/`, `models/` (image tests need `raw_images/`) |

`GLYCO_RAW_IMAGES=/path/to/folder` makes the code read `grey/` from another location.

## Stages and modules

```
raw_images/grey ──extraction──▶ data/*.csv ──modelling──▶ models/*.joblib ──results──▶ results/*.csv, RESULTS.md
```

| Module | Content |
|---|---|
| `config.py` | paths, image geometry, extraction rules (every tunable that affects a number) |
| `cohorts.py` | the three cohorts: where their tables live, expected counts, label order |
| `schema.py` | table columns and dtypes; the only readers and writers of the tables |
| `extraction.py` | image reading, pixel mask, patch means, annotation parsing, cohort checks |
| `selection.py` | which rows belong to which analysis (`sample_type` only) |
| `modelling.py` | leave-one-specimen-out cross-validation, final fit, projection, origin prediction |
| `mouse.py`, `pdac_tma.py` | cohort-specific rules (section grouping; fucose-versus-mannose phenotype and statistics) |
| `results.py` | the result tables behind the summary figures; fits nothing |
| `report.py` | writes `results/RESULTS.md` |
| `pipeline.py`, `__main__.py` | the stages as functions, and the command line |

## Design rules the tests enforce

* The held-out unit is the specimen (patient or mouse section); no specimen is in both the
  training and the test side of any fold.
* Metastases and adjacent normal tissue never reach a model fit.
* Rows are selected by `sample_type`, never by a label or an identifier string.
* Missing values stay missing: intensities are `float64`, blank patches are `NaN` in all five
  channels, and `is_tumour` is missing (not false) where there is no annotation.
* No random split, random seed or seed search anywhere; the result step fits no model.

## Reference results

`results/RESULTS.md` holds the report of the reference run (2026-10-08); the CSV tables beside
it are described in `../docs/DATA_DICTIONARY.md`.
