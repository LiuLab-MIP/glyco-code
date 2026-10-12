# Demo

Three demonstrations, from quickest to most complete. Times were measured on the machine
described in `INSTALL.md` (macOS, Apple M1, 16 GB) and include everything except installation.

## Demo A: regenerate every result from the bundled patch tables (no images needed)

```bash
cd image_analysis
source .venv/bin/activate
python -m glyco all --skip-extract
```

**Typical run time: 6 seconds.**

What it does: loads `data/*_patches_wide.csv`, runs leave-one-specimen-out cross-validation
for each analysis, fits and serialises the classifiers to `models/`, scores the metastatic
lesions, and writes 34 result tables and `results/RESULTS.md`.

Expected output (abridged):

```
Nine-class cancer classification: 72 held-out specimens, 22,667 patches, 9 classes
  accuracy 0.8924   kappa 0.8768
  per-specimen accuracy: median 0.987, IQR 0.854-1.000, min 0.000, 8 below 0.5

Breast cancer subtype: 13 held-out specimens, 4,023 patches, 4 classes
  accuracy 0.6893   kappa 0.5825

PDAC glycophenotype: 7 held-out specimens, 839 patches, 2 classes
  accuracy 0.9976   kappa 0.9952

Murine organ classification: 15 held-out specimens, 3,840 patches, 5 classes
  accuracy 0.9974   kappa 0.9967

metastatic origin:
label predicted  patches    percent
   BC        BC      805  95.153664
  CRC       CRC      833  94.018059
  ...
negative control:   BC -> LC 0 patches, CRC -> LC 0 patches
wrote 34 result tables and RESULTS.md
```

Check the result automatically with `python -m pytest`: it asserts these headline numbers,
that no specimen is ever in both the training and the test side of a fold, and that
metastases and adjacent-normal tissue never reach a model fit.

## Demo B: extract the patch tables from raw images

Needs the grey images: download the five `glyco-code_raw_images_grey_*.zip` files (788 MB together) from the
Zenodo deposit and unzip them all in `image_analysis/raw_images/` (commands in its README).

```bash
python -m glyco extract --cohort mouse      # the smallest cohort
python -m glyco extract                     # all three cohorts
python -m glyco all                         # extraction followed by Demo A
```

**Typical run times:** mouse cohort 7 seconds; the human cohort (164 fields of view) 14 seconds;
all three cohorts 18 seconds; `python -m glyco all` 22 seconds in total.

Extraction checks the cohort against the expected counts and stops with a message if a
specimen is missing. The tables it writes reproduce the ones shipped in `data/` (the test suite
re-extracts the mouse and second pancreatic cohorts and compares them to a relative tolerance
of 1e-12; the numbers come out identical).

## Demo C: bioinformatics pipeline

```bash
cd bioinformatics
bash run_all.sh                 # steps 02 to 07 and the comparison with the reference tables
```

**Typical run time: about 2 minutes 30 seconds** (the Moffitt step 05 adds about 15 seconds to the 131 to 135 seconds of the other steps), of which about 18 seconds is the GEO download
and 45 seconds each DESeq2 and enrichment. It ends with:

```
ok    clustering/ssGSEA_scores_and_subtypes.csv        max relative difference 0
ok    DEG/DESeq2_results.csv                           max relative difference 6.7e-14
ok    Survival/Cox_regression_summary.csv              max relative difference 0
ok    GEO_Validation/GEO_scores_and_subtypes.csv       max relative difference 0
ok    GEO_Validation/GEO_summary_stats.csv             max relative difference 0
ok    Enrichment/GO_Enrichment_Up_FucHigh.csv          max relative difference 0
ok    Enrichment/GO_Enrichment_Down_ManHigh.csv        max relative difference 0
ok    Enrichment/GSEA_Hallmark_results.csv             15 of 15 reference pathways recovered, max |NES difference| 0.020
ok    Enrichment/KEGG_Enrichment_Up_FucHigh.csv        100% of the 21 reference pathways recovered
ok    Enrichment/KEGG_Enrichment_Down_ManHigh.csv      100% of the 4 reference pathways recovered
ok    Moffitt/scores_per_sample.csv                    max relative difference 0
(and the Moffitt group, call and survival tables)
ALL TABLES AGREE WITH THE REFERENCE
```

Notes on what differs from run to run:

* **GSEA** is permutation based and unseeded, so normalised enrichment scores vary by a few
  hundredths between runs (largest difference to the reference 0.04 in the runs made for
  this documentation), and borderline pathways can enter or leave the significant set:
  `HALLMARK_SPERMATOGENESIS` and `HALLMARK_TNFA_SIGNALING_VIA_NFKB` were significant in
  addition to the 15 reference pathways in one run, only the former in another. The
  comparison script therefore checks that all 15 reference pathways are recovered with a
  similar score, not that the table is identical.
* **KEGG** (step 04) reads the annotation from `rest.kegg.jp`. If the server cannot be
  reached the script reports it and continues without the KEGG figures and tables; in one run
  for this documentation the connection failed with an SSL error, and in a later run it
  worked. The comparison script then reports the KEGG tables as skipped rather than failed.
* **GO and KEGG tables.** Step 04 writes the full enrichment tables (`GO_Enrichment_*.csv`,
  `KEGG_Enrichment_*.csv`) beside the dot plots. The GO tables in `expected_results/` were
  produced by this version of the script, not by an independent method; they are deterministic
  given the differential-expression table and the annotation package.
* Step 07 downloads GSE71729 on first use and caches it in `data/GEO_cache/`.

## What the demos cover

Demo A exercises every Python stage downstream of the images, Demo B the image reading,
masking, annotation parsing and table writing, and Demo C every R script except the full GDC
download (`01_data_download.R`, optional, about 800 MB). Together with the bundled
data they reproduce all reported numbers; see `METHODS_SOFTWARE.md` for the description of the
test data set.
