# Data dictionary

## Cohorts

| Cohort | Specimens | Fields of view | Content |
|---|---|---|---|
| `human`, sample type `tumour` | 72 patients | 131 | Primary tumours of nine cancer types: HCC (8 patients), CRC (8), LC (7), GC (7), CCA (8), EC (8), PCa (6), BC (13), PDAC (7) |
| `human`, sample type `metastasis` | 4 patients | 8 | Lung metastases: 2 of colorectal and 2 of breast origin. `label` is the tissue of **origin**. |
| `human`, sample type `normal` | 25 blocks | 25 | Adjacent normal tissue, labelled with the cancer type it is adjacent to (3 per type, 2 for HCC and CCA) |
| `mouse`, sample type `section` | 15 sections | 15 | Heart, kidney, liver, lung, spleen; 3 sections each, one section per organ per animal |
| `pdac_tma`, sample type `tumour` | 15 specimens | 15 | A tissue microarray of annotated pancreatic cases (the second pancreatic cohort), one field of view each |

Cancer-type codes: HCC hepatocellular carcinoma, CRC colorectal cancer, LC lung cancer, GC
gastric cancer, CCA cholangiocarcinoma, EC oesophageal cancer, PCa prostate cancer, BC breast
cancer, PDAC pancreatic ductal adenocarcinoma.

**Membership of an analysis is decided by `sample_type` and nothing else.** Metastases carry
the labels `CRC` and `BC`, the same values as two training classes, because that is their
true origin; selecting by label would silently fold them back into training.

## Identifiers

| Pattern | Meaning |
|---|---|
| `HCC-P01` | Tumour specimen number 1 of cancer type HCC |
| `CRC-LM-P01` | Metastasis specimen of colorectal origin |
| `CRC-N-P01` | Adjacent-normal block next to a colorectal cancer |
| `PDAC-TMA-01` | Specimen of the second pancreatic cohort |
| `heart-1` | Mouse heart section 1 |
| `HCC-P01-2` | Field of view 2 of specimen `HCC-P01` (`view_id`) |

The numbers are arbitrary pseudonyms. They carry no information about the patient, and
adjacent-normal blocks are pseudonymised independently of tumour blocks (a normal block and a
tumour block can come from the same patient). Hospital pathology numbers are not part of the
release. The numbering follows the sorted order of the original identifiers, but the original
identifiers are not distributed.

A patient contributes two or three fields of view of a single section; they share tissue,
staining run, imaging session and label, which is why the **patient is the unit of
cross-validation**. Mouse sections come from different animals within a class, so the
section is the unit there.

## Raw images

```
raw_images/grey/<cohort>/<sample_type>/<label>/<specimen_id>/<view>/
    GlcNAc.jpg  Gal.jpg  SA.jpg  Fuc.jpg  Man.jpg        8-bit grey, 1024 x 1024, one per terminal monosaccharide
    annotation.jpg                                        264 x 264 pathologist annotation (tumour and metastasis, PDAC TMA)
raw_images/colour/<cohort>/<sample_type>/<label>/<specimen_id>/<view>/
    GlcNAc.jpg ... Man.jpg   pseudo-coloured channels
    brightfield.jpg          bright-field image of the same field
    merge.jpg                the five channels merged
    he.jpg                   H&E image of the adjacent section, registered to the field (where available)
```

* `grey/` is what the code reads; `colour/` is for viewing and is not read by any code.
* The channel order is fixed: GlcNAc, Gal, SA, Fuc, Man.
* `annotation.jpg` encodes the pathologist's call as a colour per 16 x 16 px block: cancer
  `(255, 255, 212)`, non-cancer `(203, 203, 227)`, with a tolerance of `(10, 10, 30)`. Any
  other colour stops the extraction with the file name and pixel position.
* Normal blocks, mouse sections and the PDAC TMA cases have a single view named `1`.
* Raw images are not tracked by git (`.gitignore`); see `image_analysis/raw_images/README.md`.

**Curation notes.** The tree holds exactly the fields of view that were analysed. Left out at
curation: two view folders that are not fields of view (one colorectal, one gastric specimen),
one adjacent-normal colorectal block with a non-standard name, and two further views per
breast-metastasis specimen that were not analysed; the code therefore needs no exclusion list.
The grey images are the exact files read by the analysis. Five of them (three prostate-tumour
and two adjacent-normal channel images: `PCa-P02/2/GlcNAc`, `PCa-P04/1/GlcNAc`,
`PCa-P04/2/GlcNAc`, `EC-N-P02/1/Fuc`, `EC-N-P02/1/Man`) differ from the corresponding
black-and-white files of the acquisition export; the analysed files were kept so that results
are unchanged. Seven views have no H&E image (`PCa-P06/1`, `BC-P01/1`, `BC-P01/2`,
`BC-P11/1`, `CRC-LM-P02/1`, `CRC-LM-P02/2`, `PDAC-N-P03/1`). The colour renderings of twelve
early specimens come from an earlier acquisition export.

## Metadata

`image_analysis/metadata/specimens.csv`: `specimen_id`, `cohort`, `sample_type`, `label`,
`subtype`. The sub-classification: for BC the clinical (IHC) subtype `HER2`, `LuminalA`,
`LuminalB`, `TNBC`; for first-cohort PDAC the glycophenotype `TypeI` / `TypeII` assigned in
this study by principal-component analysis of the five channels; empty otherwise. For the PDAC
TMA cohort it is empty here and assigned by rule at extraction.

## Patch tables (`image_analysis/data/`)

Per cohort: `<cohort>_patches_wide.csv`, `<cohort>_patches_long.csv`, `<cohort>_cases.csv`.
The long table is always derived from the wide one (`melt`), never measured again.
Intensities are `float64` and missing values are empty CSV fields; they are never zero.

| Column | Meaning |
|---|---|
| `specimen_id`, `label`, `sample_type`, `subtype` | see above |
| `view_id` | `<specimen_id>-<view>` |
| `grid_row`, `grid_col` | position of the patch, 1 to 16 |
| `n_valid_px` | pixels of the 64 x 64 patch that passed the signal rule |
| `is_tumour` | pathologist's annotation of the patch; **missing where there is no annotation** (adjacent normal, mouse), never `False` |
| `GlcNAc`, `Gal`, `SA`, `Fuc`, `Man` | mean 8-bit intensity over the signal pixels; missing for a blank patch (all five at once) |
| `glycan`, `intensity` | (long table) channel name and its value |

Extraction rule: a pixel is kept when any channel is at least `T` (10 for tumour and
metastasis, 5 for adjacent normal and for mouse sections) and not above 250 in all five
channels; a patch with fewer than 512 kept pixels (1 for mouse) is blank. Tumour-only
analyses use the annotated region where there is an annotation and the whole section where
there is none.

`<cohort>_cases.csv` has one row per specimen: `specimen_id`, `cohort`, `sample_type`,
`label`, `subtype`, `n_views`, `n_patches`, `n_tumour_patches`, `n_blank_patches`,
`has_annotation`.

## Result tables (`image_analysis/results/`)

| File | Content |
|---|---|
| `loso_summary.csv` | accuracy, Cohen's κ, and the distribution of per-specimen accuracy for each analysis (`primary`, `subtype`, `phenotype`, `mouse`) |
| `metrics_<analysis>_loso.csv` | per-class precision, recall (= sensitivity), F1, support and one-versus-rest specificity |
| `confusion_<analysis>_loso.csv` | pooled out-of-fold confusion matrix (counts; rows are the reference) |
| `per_case_accuracy_<analysis>.csv` | one row per held-out specimen: `specimen_id`, `label`, `n_patches`, `n_correct`, `accuracy` |
| `origin_assignment.csv`, `origin_negative_control.csv` | assignment of the metastatic patches to each class, and the number assigned to lung |
| `lda_scores_<analysis>.csv`, `lda_explained_variance.csv` | discriminant coordinates `LD1..` of every patch, and their explained variance |
| `profile_cancer_type.csv` | mean profile per cancer type, row-standardised (`z_*`), in clustered order (`cluster_position`; Ward linkage, Euclidean distance); `profile_cancer_type_linkage.csv` holds the dendrogram merges and heights |
| `profile_specimen_primary.csv` | mean profile per patient, standardised against the whole table |
| `profile_section_mouse.csv` | mean profile per mouse section, row-standardised |
| `profile_specimen_pdac_tma.csv` | mean profile per specimen of the second cohort, row-standardised, in Ward clustered order; `profile_specimen_pdac_tma_linkage.csv` holds its merges |
| `pdac_pca_scores.csv`, `_variance.csv`, `_loadings.csv`, `_specimen_summary.csv` | principal components of the first pancreatic cohort and the share of each patient's patches above the PC2 boundary 0.15 |
| `pdac_tma_phenotype_stats.csv` | per channel, patch level: median per phenotype, Mann-Whitney P and rank-biserial r (patches treated as independent) |
| `pdac_tma_phenotype_stats_by_case.csv` | the same comparison with one value per specimen (mean over its tumour patches, 9 against 6 specimens), two-sided exact Mann-Whitney U test; fucose and mannose define the phenotype, so only the other three channels are an independent comparison |
| `pdac_tma_first_cohort_agreement.csv` | whether the Fuc-versus-Man rule reproduces the first cohort's labels |
| `violin_group_summary.csv` | n, mean, median, quartiles, range of the intensity distributions, untrimmed (`none`) and after the display-only 1.5 × IQR trim (`iqr_1.5`) |
| `RESULTS.md` | the report of the run |

## Models (`image_analysis/models/`)

`lda_primary.joblib`, `lda_bc_subtype.joblib`, `lda_pdac_phenotype.joblib`,
`lda_mouse_organ.joblib`: a dictionary with the fitted scikit-learn estimator, the feature
order, the class order, the specimens it was fitted on and the explained variance. Loading a
joblib file executes code; load only files you trust, or regenerate them with
`python -m glyco analyse`. Files are tied to scikit-learn 1.5.2.
