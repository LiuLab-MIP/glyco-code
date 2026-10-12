# Changelog

## 1.0.0 (2026-10-08): first public release

Code accompanying "Visualising and deciphering the glyco-code".

* `image_analysis/`: extraction of patch-level terminal-monosaccharide intensities from
  fluorescence images, leave-one-patient-out linear discriminant analysis, metastatic origin
  prediction, murine organ classification, second pancreatic cohort, and the result tables
  behind every summary figure. Plotting code is not part of the release.
* `bioinformatics/`: TCGA-PAAD and GSE71729 analysis in R (ssGSEA glycophenotype scoring,
  differential expression, enrichment, survival, external validation, and the comparison with
  the Moffitt basal-like and classical subtypes).
* Raw images are identified by pseudonymous specimen identifiers; hospital pathology numbers
  are not part of the release.
