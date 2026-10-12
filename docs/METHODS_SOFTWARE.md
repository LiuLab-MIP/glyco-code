# Software description for the Methods section

*(Draft text covering the three items Nature asks for in the Methods: the key operations of the
code, its key characteristics, and the test data set. Adapt to the journal's length limits.)*

## 1. Key operations performed by the code

**Image analysis (Python).** Each field of view is a set of five 8-bit grey images
(1024 × 1024 pixels), one per terminal monosaccharide (GlcNAc, Gal, SA, Fuc, Man). The code

1. divides every image into a 16 × 16 grid of 64 × 64 pixel patches;
2. keeps a pixel when at least one of the five channels reaches an intensity threshold (10
   for malignant tissue and metastases, 5 for adjacent normal tissue and for mouse sections)
   and drops a pixel that exceeds 250 in all five channels;
3. records, per patch, the mean intensity of each channel over the kept pixels, and marks a
   patch with fewer than 512 kept pixels (1 for mouse sections) as missing in all five channels;
4. reads the pathologist's annotation (cancer or non-cancer, one value per patch) and, for
   annotated specimens, restricts analyses to the annotated tumour patches;
5. classifies patches from the five mean intensities by linear discriminant analysis
   (singular-value solver, class priors from the data, no feature scaling) and estimates
   performance by leave-one-specimen-out cross-validation, pooling the out-of-fold
   predictions to report accuracy, Cohen's κ, per-class precision, recall and F1, and the
   distribution of per-specimen accuracy;
6. fits a final classifier on all primary-tumour specimens and uses it to project the data
   (discriminant coordinates) and to score the metastatic lesions, which were never part of
   any training set; lung is the organ of residence of every lesion, so patches assigned to
   lung form a negative control;
7. assigns the glycophenotype of a second pancreatic cohort (a tissue microarray) by comparing each specimen's mean
   fucose and mannose intensity over its annotated tumour patches, and compares the two
   phenotypes per channel with a two-sided Mann-Whitney test and the rank-biserial
   correlation;
8. writes the tables behind the summary figures: standardised mean profiles with
   hierarchical-clustering order, principal components of the first pancreatic cohort, and
   a specimen-level (exact Mann-Whitney) comparison of the two phenotypes, and per-group quartiles of the intensity distributions (also after a display-only trim of
   points beyond 1.5 × the interquartile range, which is never applied to a classifier or a
   test).

**Bioinformatics (R).** From TCGA-PAAD STAR gene counts (178 primary tumours), the code
computes ssGSEA scores of a 16-gene fucosylation signature, a 7-gene mannose-processing set
and a 9-gene mannose-biosynthesis set; defines a mannose composite (biosynthesis minus
processing) and a glycan shift score (z-score of fucosylation minus z-score of the mannose
composite); splits patients at the median into Fuc-favored and Man-favored groups; runs
DESeq2 between the groups (450 up- and 523 down-regulated genes at |log2 fold change| > 1 and
adjusted P < 0.05), GO, KEGG and Hallmark GSEA enrichment, Kaplan-Meier and multivariable Cox
analysis; repeats the scoring and survival analysis on GSE71729 (357 samples scored, 125
with survival data), and scores every tumour against the 25 basal-like and 25 classical
Moffitt tumour genes (ssGSEA, expressed as the difference of the two z-scored signatures) and
with a 21-gene-pair single-sample classifier, to relate the glycophenotypes to the Moffitt
subtypes.

## 2. Key characteristics

* **Algorithms.** Patch statistics are plain array operations (NumPy). Classification is
  scikit-learn's `LinearDiscriminantAnalysis(solver="svd")` with `LeaveOneGroupOut`. The
  statistics use SciPy (`mannwhitneyu`, hierarchical clustering) and scikit-learn (PCA). The
  R analysis uses GSVA (ssGSEA), DESeq2, clusterProfiler, survival and survminer.
* **Held-out unit.** The unit is the patient for human data (a patient contributes two or
  three fields of view of one section, so holding out a view would leave its siblings in
  training) and the section for the murine data (each animal contributed one section per
  organ, so no two sections of a class share an animal).
* **Safeguards built into the code and enforced by tests.** Membership of an analysis is
  decided by sample type, never by label or identifier string; metastases and adjacent
  normal tissue never reach a model fit; blank patches stay missing (never zero); no random
  split, random seed or seed search exists anywhere; the result step fits no model and reads
  the serialised classifiers; every held-out fold contains exactly one specimen.
* **Determinism.** The Python pipeline uses no random numbers. Re-extracting the images
  reproduces the shipped tables exactly, and re-running the analysis reproduces the result
  tables. GSEA in the R pipeline is permutation based and unseeded, so its enrichment scores
  vary slightly between runs (see `DEMO.md`).
* **Dependencies.** Python 3.12 with NumPy 2.1.3, pandas 2.2.3, scikit-learn 1.5.2, SciPy
  1.14.1, Pillow 11.0.0 and joblib 1.4.2; R 4.6.1 with Bioconductor 3.23. See `INSTALL.md`.
* **Resources.** No special hardware. Extraction of all images takes 18 seconds and the whole
  analysis 6 seconds on a laptop (Apple M1, 16 GB); peak memory 0.4 GB.
* **Availability.** MIT licence; code, data tables, models and results are public; raw
  images are deposited separately on Zenodo (see the root `README.md`).
* **Community input.** The code was not submitted to community commenting.

## 3. Test data set

The data set bundled with the code lets a reader reproduce every reported number without the
raw images:

* `image_analysis/data/`: patch tables of the 101 human specimens (72 primary tumours from
  72 patients and 131 fields of view; 4 lung metastases; 25 adjacent-normal blocks), the 15
  mouse organ sections and the 15 specimens of the second pancreatic cohort, 49,664 patches
  in total, with the pathologist's annotation per patch. Specimens carry pseudonymous
  identifiers.
* `image_analysis/models/` and `results/`: the fitted classifiers and every result table of
  the reference run, to compare against.
* `bioinformatics/data/TCGA_PAAD_processed.rds` (178 TCGA-PAAD primary tumours; counts and
  clinical data) and `expected_results/`; GSE71729 is downloaded from GEO by the script.

Running the code on these inputs (`DEMO.md`) reproduces the leave-one-patient-out accuracy of
89.2 % (κ 0.877) for nine cancer types, 68.9 % (κ 0.583) for breast subtypes, 99.8 %
(κ 0.995) for the pancreatic glycophenotypes, 99.7 % (κ 0.997) for mouse organs, the
assignment of the metastatic lesions to their origin (94.0 % colorectal, 95.2 % breast, none
to lung), and the R tables listed in `bioinformatics/expected_results/`. The raw images
(`image_analysis/raw_images/`, deposited separately) additionally allow the extraction step
to be run. External dependencies are limited to PyPI, CRAN, Bioconductor, GDC and GEO as
specified in `INSTALL.md`.
