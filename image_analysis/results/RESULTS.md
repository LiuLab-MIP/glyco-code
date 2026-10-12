# Results — glyco-code image analysis

Every number here is recorded with the grouping that produced it and the date of the run.
The whole file is regenerated in one pass, so no two numbers in it can disagree.

**Run: 2026-10-12.**

---

## Summary

| analysis | accuracy |  | folds |
|---|---|---|---|
| Nine-class cancer classification | **89.2%** | κ 0.877 | 72 patients |
| Breast cancer subtype | **68.9%** | κ 0.583 | 13 patients |
| PDAC glycophenotype | **99.8%** | κ 0.995 | 7 patients |
| Metastatic origin, colorectal | **94.0%** | assigned correctly | 2 patients |
| Metastatic origin, breast | **95.2%** | assigned correctly | 2 patients |
| Negative control | **0** | metastatic patches assigned to lung | — |
| Murine organ classification | **99.7%** | κ 0.997 | 15 sections |

Every human classification is cross-validated by leaving out one patient at a time. A case
contributes two or three fields of view of a single section (same tissue, same staining
run, same imaging session), so holding out a field of view would leave its siblings in
training and measure recognition of a section the model had already seen. The murine cohort
is grouped by section instead; section 6 says why that is the right unit there.

---

## 1 — Extraction

A pixel is kept when any of the five channels reaches the threshold below; a patch with
fewer than 512 such pixels has too little tissue to average and is
recorded as missing in all five channels. In every tissue, a pixel above 250
in all five channels is also dropped, because a pixel saturated in every channel says
nothing about the ratio between them.

| tissue | pixel threshold | minimum pixels per patch |
|---|---|---|
| tumour | 10 | 512 |
| metastasis | 10 | 512 |
| normal | 5 | 512 |

The thresholds differ deliberately. Malignant tissue is bright, and 10 is where background
stops being counted there; it is also the value every classification result was measured
with. Adjacent-normal tissue is dimmer overall, so the same threshold would discard pixels
that carry real signal, and because it would discard the dim ones, the mean over what
survived would be pushed up. Each tissue is thresholded for its own brightness.

| tissue | patches recorded missing | patches | % |
|---|---|---|---|
| metastasis | 1 | 2,048 | 0.049% |
| normal | 30 | 6,400 | 0.469% |
| tumour | 56 | 33,536 | 0.167% |

Of the 56 tumour patches recorded missing,
**0** lie inside the pathologist's annotation; those are the only ones that would
otherwise have reached a classifier.

### Cohort

```
                   patients  views  patches  blank
sample_type label                                 
metastasis  BC            2      4     1024      0
            CRC           2      4     1024      1
normal      BC            3      3      768      0
            CCA           2      2      512      0
            CRC           3      3      768      0
            EC            3      3      768      0
            GC            3      3      768      0
            HCC           2      2      512      0
            LC            3      3      768     19
            PCa           3      3      768     11
            PDAC          3      3      768      0
tumour      BC           13     29     7424     37
            CCA           8     13     3328      9
            CRC           8     14     3584      0
            EC            8     14     3584      9
            GC            7     11     2816      1
            HCC           8     17     4352      0
            LC            7     11     2816      0
            PCa           6     12     3072      0
            PDAC          7     10     2560      0
```

---

## 2 — Nine-class cancer classification

72 patients, 131 fields of view, 22,667 patches, 9 classes. Leave-one-patient-out, 72 folds.

**Accuracy 0.8924, Cohen's κ 0.8768.**

Per-patient accuracy: median 0.987, IQR 0.854–1.000, range 0.000–1.000, 8 of 72 below 0.5. This distribution is the account of variability; there is no held-out split to quote instead (`per_case_accuracy_*.csv`).

| class | precision | recall (sensitivity) | specificity | F1 | patches |
|---|---|---|---|---|---|
| BC | 0.884 | 0.932 | 0.974 | 0.907 | 4,023 |
| CCA | 0.774 | 0.856 | 0.973 | 0.813 | 2,206 |
| CRC | 0.983 | 0.949 | 0.998 | 0.966 | 2,564 |
| EC | 0.856 | 0.993 | 0.981 | 0.920 | 2,289 |
| GC | 0.985 | 0.746 | 0.999 | 0.849 | 2,320 |
| HCC | 0.951 | 0.883 | 0.990 | 0.916 | 4,104 |
| LC | 0.954 | 0.998 | 0.996 | 0.976 | 1,858 |
| PCa | 0.908 | 0.954 | 0.988 | 0.930 | 2,464 |
| PDAC | 0.430 | 0.386 | 0.980 | 0.407 | 839 |

Row-normalised confusion matrix (%), rows = reference:

| reference \ prediction | BC | CCA | CRC | EC | GC | HCC | LC | PCa | PDAC |
|---|---|---|---|---|---|---|---|---|---|
| **BC** | 93.2 | 4.2 | 0.0 | 0.2 | 0.0 | 2.3 | 0.0 | 0.0 | 0.0 |
| **CCA** | 5.9 | 85.6 | 0.0 | 0.0 | 0.0 | 3.8 | 1.9 | 2.7 | 0.0 |
| **CRC** | 0.0 | 0.0 | 94.9 | 0.0 | 0.0 | 0.0 | 0.0 | 2.7 | 2.4 |
| **EC** | 0.0 | 0.0 | 0.0 | 99.3 | 0.6 | 0.0 | 0.1 | 0.0 | 0.0 |
| **GC** | 0.0 | 0.0 | 0.0 | 14.4 | 74.6 | 0.0 | 0.0 | 0.0 | 10.9 |
| **HCC** | 8.4 | 2.3 | 0.0 | 0.0 | 0.0 | 88.3 | 1.0 | 0.0 | 0.0 |
| **LC** | 0.0 | 0.2 | 0.0 | 0.0 | 0.0 | 0.0 | 99.8 | 0.0 | 0.0 |
| **PCa** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 95.4 | 4.6 |
| **PDAC** | 2.4 | 33.6 | 4.8 | 4.8 | 1.5 | 1.2 | 0.0 | 13.1 | 38.6 |

Weakest patients:

| patient | accuracy |
|---|---|
| PDAC-P03 | 0.0000 |
| PDAC-P01 | 0.0000 |
| PDAC-P06 | 0.0893 |
| PDAC-P07 | 0.1408 |
| GC-P07 | 0.2480 |
| PDAC-P05 | 0.3585 |

**PDAC recall is 38.6%**, the weakest class here.

---

## 3 — Breast cancer subtype

13 patients, 29 fields of view, 4,023 patches, 4 classes. Leave-one-patient-out, 13 folds.

**Accuracy 0.6893, Cohen's κ 0.5825.**

Per-patient accuracy: median 0.739, IQR 0.474–0.790, range 0.157–1.000, 4 of 13 below 0.5. This distribution is the account of variability; there is no held-out split to quote instead (`per_case_accuracy_*.csv`).

| class | precision | recall (sensitivity) | specificity | F1 | patches |
|---|---|---|---|---|---|
| HER2 | 1.000 | 0.879 | 1.000 | 0.935 | 1,006 |
| LuminalA | 0.597 | 0.680 | 0.816 | 0.636 | 1,150 |
| LuminalB | 0.486 | 0.528 | 0.809 | 0.506 | 1,025 |
| TNBC | 0.791 | 0.672 | 0.953 | 0.727 | 842 |

Row-normalised confusion matrix (%), rows = reference:

| reference \ prediction | HER2 | LuminalA | LuminalB | TNBC |
|---|---|---|---|---|
| **HER2** | 87.9 | 0.0 | 11.9 | 0.2 |
| **LuminalA** | 0.0 | 68.0 | 30.7 | 1.3 |
| **LuminalB** | 0.0 | 34.2 | 52.8 | 13.0 |
| **TNBC** | 0.0 | 21.0 | 11.8 | 67.2 |

Weakest patients:

| patient | accuracy |
|---|---|
| BC-P05 | 0.1569 |
| BC-P10 | 0.2140 |
| BC-P01 | 0.3647 |
| BC-P06 | 0.4745 |
| BC-P04 | 0.5531 |
| BC-P12 | 0.6516 |

**LuminalB recall is 52.8%**, the weakest class here.

---

## 4 — PDAC glycophenotype

7 patients, 10 fields of view, 839 patches, 2 classes. Leave-one-patient-out, 7 folds.

**Accuracy 0.9976, Cohen's κ 0.9952.**

Per-patient accuracy: median 1.000, IQR 1.000–1.000, range 0.972–1.000, 0 of 7 below 0.5. This distribution is the account of variability; there is no held-out split to quote instead (`per_case_accuracy_*.csv`).

| class | precision | recall (sensitivity) | specificity | F1 | patches |
|---|---|---|---|---|---|
| TypeI | 1.000 | 0.995 | 1.000 | 0.997 | 382 |
| TypeII | 0.996 | 1.000 | 0.995 | 0.998 | 457 |

Row-normalised confusion matrix (%), rows = reference:

| reference \ prediction | TypeI | TypeII |
|---|---|---|
| **TypeI** | 99.5 | 0.5 |
| **TypeII** | 0.0 | 100.0 |


**Read this with care.** The two glycophenotypes were assigned by a principal-component
analysis of these same five channels, so a classifier that recovers them from those
channels is close to tautological. The number says the assignment is internally consistent
and reproducible patient by patient. It is not independent validation of the phenotypes.

---

## 5 — Metastatic origin

Four lesions, 1,732 patches, scored against the classifier fitted on all 72 primary-tumour patients. They are external by construction: the training selection asks for tumour tissue, so no fold and no fit ever saw them.

| true origin | predicted | patches | share |
|---|---|---|---|
| BC | BC | 805 | 95.2% |
| BC | CCA | 21 | 2.5% |
| BC | HCC | 20 | 2.4% |
| CRC | CRC | 833 | 94.0% |
| CRC | PCa | 34 | 3.8% |
| CRC | PDAC | 19 | 2.1% |

**Negative control.** Every lesion here was resected from the lung, so no metastatic patch should be assigned to lung cancer.

| origin | patches assigned LC | share |
|---|---|---|
| BC | 0 | 0.0% |
| CRC | 0 | 0.0% |

---

## 6 — Murine organ cohort

Five organs, 3 sections each, 3,840 patches. Extraction rule **mouse section (T=5, saturation >250 in all, N=1)**: sections are imaged whole, so there is no annotation to mask by and no dim-patch rule; a patch counts unless it holds no signal at all, and pixels saturated in every channel are dropped because they say nothing about the ratio between channels.

| organ | sections | patches |
|---|---|---|
| heart | 3 | 768 |
| kidney | 3 | 768 |
| liver | 3 | 768 |
| lung | 3 | 768 |
| spleen | 3 | 768 |

**The held-out unit is the section, not the animal, and that is correct.** Each mouse contributed one section per organ, so within any one organ no two sections come from the same animal: a held-out section shares no animal with any training section of its own class, and individual-level information cannot leak within a class. The human cohort is grouped by patient for the opposite reason: there, two or three fields of view come from a single section. The two schemes differ because the unit of independence sits at a different level, not because one of them is wrong.

**Accuracy 0.9974, Cohen's κ 0.9967**, leave-one-section-out, 15 folds. Per-section accuracy: median 1.0000, range 0.9805–1.0000.

| class | precision | recall (sensitivity) | specificity | F1 | patches |
|---|---|---|---|---|---|
| heart | 1.000 | 0.999 | 1.000 | 0.999 | 768 |
| kidney | 1.000 | 1.000 | 1.000 | 1.000 | 768 |
| liver | 0.999 | 1.000 | 1.000 | 0.999 | 768 |
| lung | 0.999 | 0.988 | 1.000 | 0.993 | 768 |
| spleen | 0.990 | 1.000 | 0.997 | 0.995 | 768 |

Row-normalised confusion matrix (%), rows = reference:

| reference \ prediction | heart | kidney | liver | lung | spleen |
|---|---|---|---|---|---|
| **heart** | 99.9 | 0.0 | 0.0 | 0.1 | 0.0 |
| **kidney** | 0.0 | 100.0 | 0.0 | 0.0 | 0.0 |
| **liver** | 0.0 | 0.0 | 100.0 | 0.0 | 0.0 |
| **lung** | 0.0 | 0.0 | 0.1 | 98.8 | 1.0 |
| **spleen** | 0.0 | 0.0 | 0.0 | 0.0 | 100.0 |

---

## 7 — Second pancreatic cohort

15 specimens, one annotated field of view each, 3,840 patches. Measured with the same rule as the rest of the malignant tissue (**cancer (T=10, saturation >250 in all, N=512)**). All summaries use the annotated tumour patches only.

**How the glycophenotype is assigned.** Per specimen, over its annotated tumour patches: mean fucose above mean mannose is Fuc-high, otherwise Man-high. The split is 9 FucHigh, 6 ManHigh.

Comparison between the two phenotypes, per channel, computed on the patches as selected, before any display trimming:

| channel | Fuc-high median | Man-high median | Mann-Whitney | rank-biserial r |
|---|---|---|---|---|
| GlcNAc | 20.5 | 21.2 | P = 0.022 | -0.064 |
| Gal | 34.6 | 35.9 | P = 0.024 | -0.063 |
| SA | 85.2 | 86.6 | P = 0.014 | -0.068 |
| Fuc | 52.7 | 28.2 | P < 0.001 | +0.884 |
| Man | 32.6 | 55.4 | P < 0.001 | -0.935 |

The distinction rests on two channels and only two: fucose and mannose separate almost completely, while the other three barely move.

**Specimen as the unit.** The patches of one specimen are not independent, so the same comparison is repeated with one value per specimen (the mean over its tumour patches; 9 Fuc-high against 6 Man-high specimens), two-sided exact Mann-Whitney test:

| channel | Fuc-high median | Man-high median | exact Mann-Whitney | rank-biserial r |
|---|---|---|---|---|
| GlcNAc | 20.7 | 21.4 | P = 0.955 | +0.037 |
| Gal | 35.4 | 33.6 | P = 0.607 | +0.185 |
| SA | 80.6 | 81.5 | P = 0.776 | +0.111 |
| Fuc | 53.8 | 27.8 | P = 0.0008 | +0.963 |
| Man | 33.5 | 51.9 | P = 0.0004 | -1.000 |

Fucose and mannose define the phenotype (fucose above mannose), so their separation is partly built in; the other three channels are the independent part of the comparison.

**This addresses the circularity noted for the first cohort.** Applying the same two-channel rule to those 7 specimens reproduces **7 of 7** of the labels a principal-component analysis gave them. The phenotypes are not an artefact of the method that found them: a rule anyone can check by eye recovers the same grouping.

---

## 8 — Files written

| file | what it is |
|---|---|
| `data/<cohort>_patches_wide.csv` | one row per patch (human, mouse, pdac_tma) |
| `data/<cohort>_patches_long.csv` | the same rows, one per channel |
| `data/<cohort>_cases.csv` | one row per specimen |
| `models/lda_primary.joblib` | fitted on all 72 primary-tumour patients |
| `models/lda_bc_subtype.joblib` | fitted on all 13 breast patients |
| `models/lda_pdac_phenotype.joblib` | fitted on all 7 pancreatic patients |
| `models/lda_mouse_organ.joblib` | fitted on all 15 murine sections |
| `results/loso_summary.csv` | accuracy, κ and per-specimen spread of every cross-validation |
| `results/metrics_*_loso.csv`, `confusion_*_loso.csv` | per-class precision, recall, F1, support; confusion matrices |
| `results/per_case_accuracy_*.csv` | out-of-fold accuracy for every held-out specimen |
| `results/origin_assignment.csv` | metastatic assignment counts and shares |
| `results/lda_scores_*.csv` | discriminant coordinates of every patch |
| `results/profile_*.csv` | mean profiles, raw and standardised, with clustering order |
| `results/pdac_pca_*.csv` | principal components of the first pancreatic cohort |
| `results/pdac_tma_*.csv` | second pancreatic cohort: phenotype comparison and agreement |
| `results/violin_group_summary.csv` | distribution summaries, before and after the 1.5 × IQR display trim |
