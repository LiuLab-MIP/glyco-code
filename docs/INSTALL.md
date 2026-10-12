# Installation guide

## Systems and requirements

| Item | Image analysis (Python) | Bioinformatics (R) |
|---|---|---|
| Tested operating system | macOS 27.0.1, arm64 (Apple M1, 8 cores, 16 GB RAM) | the same machine |
| Language | Python 3.12.15 (3.11 to 3.13 expected to work) | R 4.6.1, Bioconductor 3.23 |
| Package manager | `pip` (a virtual environment is recommended) | `BiocManager` and CRAN |
| Non-standard hardware | none (no GPU, no cluster) | none |
| Peak memory | 0.40 GB (extraction), 0.28 GB (analysis) | not measured; a few GB is ample |
| Disk | about 65 MB for the code, tables and models; about 3 GB more for the raw images (optional; 788 MB of it, the grey images, is what the code reads) | about 32 MB shipped (counts 10 MB, TPM 22 MB); GEO cache about 35 MB; the optional full GDC download is about 800 MB |
| Internet | only to install packages | to install packages; steps 04 (KEGG) and 07 (GEO) query remote servers |

The Python code uses only pure-Python code and binary wheels, and the R code only CRAN and
Bioconductor packages, so Linux and Windows should work, but **only macOS was tested**.

## Image analysis (Python)

```bash
cd image_analysis
python3.12 -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Pinned versions (the ones that produced the reference results):

| Package | Version | Used for |
|---|---|---|
| numpy | 2.1.3 | arrays, patch statistics |
| pandas | 2.2.3 | tables |
| scikit-learn | 1.5.2 | linear discriminant analysis, cross-validation, PCA |
| scipy | 1.14.1 | hierarchical clustering, Mann-Whitney test |
| pillow | 11.0.0 | reading JPEG images |
| joblib | 1.4.2 | serialising fitted classifiers |
| pytest | 8.3.3 | tests |

Optionally install the package itself so that `glyco` can be imported from anywhere:
`pip install -e .`

**Typical install time:** 2 minutes 7 seconds for creating the environment and installing
all pinned packages from PyPI (measured on a broadband connection; dominated by download).

Check the installation: `python -m pytest` (63 tests, about 13 seconds with the raw images
present; without the images the five image-tree tests are skipped).

## Bioinformatics (R)

```bash
cd bioinformatics
Rscript 00_install_packages.R        # or: bash run_all.sh --install
```

Versions with which the reference comparison in `docs/DEMO.md` was made:

| Package | Version | Package | Version |
|---|---|---|---|
| TCGAbiolinks | 2.40.0 | GEOquery | 2.80.0 |
| SummarizedExperiment | 1.42.0 | survival | 3.8.12 |
| S4Vectors | 0.50.3 | survminer | 0.5.2 |
| GSVA | 2.6.6 (needs 2.0 or later; also used by step 05) | ggplot2 | 4.0.3 |
| DESeq2 | 1.52.0 | pheatmap | 1.0.13 |
| clusterProfiler | 4.20.2 | dplyr | 1.2.1 |
| org.Hs.eg.db | 3.23.1 | ggrepel | 0.9.8 |
| enrichplot | 1.32.0 | ggridges | 0.5.7 |
| DOSE | 4.6.0 | svglite | 2.2.2 |
| msigdbr | 26.1.1 | broom / ggpubr | 1.0.13 / 1.0.0 |

`run_all.sh` writes the exact versions of a run to `results/sessionInfo.txt`.

**Typical install time:** 1 minute 50 seconds for the seven Bioconductor packages and their
dependencies that were missing on the test machine (TCGAbiolinks, DESeq2, clusterProfiler,
org.Hs.eg.db, enrichplot, DOSE, msigdbr), plus 10 seconds to add ggridges; the other
packages were already installed there. On a machine with a bare R installation it will take
longer.
