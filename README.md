# glyco-code

Analysis code for the article *Visualising and deciphering the glyco-code*, version 1.0.0.

Two terminal-monosaccharide readouts are analysed here:

* **Image analysis** (`image_analysis/`, Python): multicolour fluorescence images of five
  terminal monosaccharides (GlcNAc, Gal, SA, Fuc, Man) are reduced to patch-level intensity
  tables, and cancer type, breast subtype, pancreatic glycophenotype, metastatic origin and
  mouse organ are classified by linear discriminant analysis with leave-one-specimen-out
  cross-validation.
* **Bioinformatics** (`bioinformatics/`, R): the same fucose-high / mannose-high
  glycophenotypes are examined in public transcriptomes (TCGA-PAAD, GSE71729) by ssGSEA
  scoring, differential expression, enrichment and survival analysis.

**Raw images.** The 3 GB of fluorescence images are not in this repository. They are deposited
on Zenodo: [10.5281/zenodo.23316047](https://doi.org/10.5281/zenodo.23316047), as 19 zip files of at most 183 MB
each (`grey`, 5 files, 788 MB, which the code reads; `colour`, 14 files, 2.3 GB, for viewing).
Download them and unzip them all into `image_analysis/raw_images/`; the commands are in
`image_analysis/raw_images/README.md`.
Everything except the image extraction step runs without them.

## Documentation provided

| File | What it contains |
|---|---|
| `README.md` | This overview, the layout, and the index of all documentation |
| `docs/INSTALL.md` | Operating system, languages, dependencies and versions, hardware, typical install time |
| `docs/DEMO.md` | Demos on the bundled example data, expected output, typical run times |
| `docs/METHODS_SOFTWARE.md` | Draft Methods text: key operations, key characteristics, the test data set |
| `docs/DATA_DICTIONARY.md` | Raw-image layout, identifiers, table columns, curation notes |
| `LICENSE` | MIT licence |
| `CHANGELOG.md` | Version history |
| `image_analysis/README.md` | The Python package: commands, modules, tests, outputs |
| `image_analysis/raw_images/README.md` | The raw-image tree and how to obtain it |
| `bioinformatics/README.md` | The R pipeline: steps, inputs, outputs, the Moffitt comparison, tested versions |
| `bioinformatics/data/moffitt/README.md` | Source and licence of the Moffitt gene lists and classifier |

## Layout

```
glyco-code/
  docs/                 installation, demo, methods text, data dictionary
  image_analysis/
    glyco/              Python package (python -m glyco)
    tests/              pytest suite (63 tests)
    metadata/           specimens.csv: sub-classification of every specimen
    data/               patch tables produced by extraction (CSV)
    models/             fitted classifiers (joblib)
    results/            result tables (CSV) and RESULTS.md of the reference run
    raw_images/         grey/ and colour/ image trees (not tracked by git)
  bioinformatics/
    00_install_packages.R ... 07_GEO_validation.R, compare_with_expected.R
    data/               TCGA-PAAD counts and TPM (R objects), the GDC manifest, Moffitt genes and classifier
    expected_results/   reference tables for comparison
    run_all.sh
```

## Quick start

Image analysis, from the bundled patch tables (about 7 seconds, no images needed):

```bash
cd image_analysis
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m glyco all --skip-extract      # models, results/ and results/RESULTS.md
python -m pytest                        # 63 tests
```

Bioinformatics, from the bundled TCGA-PAAD object (about 2.5 minutes, needs internet for GEO and KEGG):

```bash
cd bioinformatics
bash run_all.sh --install               # --install installs the R packages first
```

Details, expected output and timings are in `docs/INSTALL.md` and `docs/DEMO.md`.

## Reference results (leave-one-specimen-out)

| Analysis | Accuracy | Cohen's κ | Held-out unit |
|---|---|---|---|
| Nine-class cancer classification | 89.2 % | 0.877 | 72 patients |
| Breast cancer subtype | 68.9 % | 0.583 | 13 patients |
| PDAC glycophenotype | 99.8 % | 0.995 | 7 patients |
| Murine organ classification | 99.7 % | 0.997 | 15 sections |
| Metastatic origin | 94.0 % (colorectal), 95.2 % (breast) of patches assigned to the true origin; 0 patches assigned to lung (negative control) | | external lesions |

The full report, with per-class metrics and the per-patient spread, is
`image_analysis/results/RESULTS.md`.

## Licence and citation

MIT licence (`LICENSE`). Please cite the article.
