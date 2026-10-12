# bioinformatics

R pipeline that examines the fucose-high / mannose-low (Fuc-favored) and mannose-high /
fucose-low (Man-favored) glycophenotypes in public pancreatic-cancer transcriptomes:
TCGA-PAAD (discovery, 178 primary tumours) and GSE71729 (validation). Installation and
timings: `../docs/INSTALL.md`, `../docs/DEMO.md`.

## Run

```bash
cd bioinformatics
bash run_all.sh                 # steps 02 to 07, then compare_with_expected.R (about 2.5 min)
bash run_all.sh --install       # install the R packages first (step 00)
bash run_all.sh --download      # rebuild the input from GDC first (step 01, about 800 MB)
```

All scripts must be run from this directory (paths are relative). Outputs go to `results/`
(not tracked by git); the reference tables are in `expected_results/`.

## Steps

| Script | Does | Reads | Writes (in `results/`) |
|---|---|---|---|
| `00_install_packages.R` | installs the packages and fails if one is missing | | |
| `01_data_download.R` | queries GDC for TCGA-PAAD STAR counts, keeps primary tumours, maps Ensembl IDs to symbols (duplicates: highest total expression) | GDC | `data/TCGA_PAAD_processed.rds` |
| `02_ssGSEA_scoring.R` | ssGSEA of three gene sets on log2(counts + 1); mannose composite = biosynthesis − processing; glyco shift score = z(Fuc) − z(Man composite); median split into Fuc-favored / Man-favored | `data/TCGA_PAAD_processed.rds` | `clustering/ssGSEA_scores_and_subtypes.csv`, heatmap; `data/TCGA_PAAD_subtyped.rds` |
| `03_DEG_analysis.R` | DESeq2 between the groups (reference Man-favored; genes with at least 10 counts in at least 10 samples) | subtyped object | `DEG/DESeq2_results.csv`, volcano plot |
| `04_enrichment_analysis.R` | GO biological process and KEGG of up- and down-regulated genes (|log2 FC| > 1, adjusted P < 0.05); Hallmark GSEA on the Wald statistic | DESeq2 table | `Enrichment/GO_Enrichment_*.csv`, `KEGG_Enrichment_*.csv`, `GSEA_Hallmark_results.csv`, dot plots, ridge plot |
| `05_moffitt_scoring.R` | scores each TCGA-PAAD tumour against the 25 basal-like and 25 classical Moffitt genes (ssGSEA, combined as the Moffitt score) and with the 21-gene-pair single-sample classifier, on TPM; compares both with the glycophenotype and with survival; repeats the scoring on raw counts as a sensitivity analysis | `data/TCGA_PAAD_tpm.rds`, `data/TCGA_PAAD_processed.rds`, `data/moffitt/`, step 02 | `Moffitt/scores_per_sample.csv`, `scores_by_glycophenotype.csv`, `calls_by_glycophenotype.csv`, `correlations_with_glyco_scores.csv`, `survival.csv`, `gene_coverage.csv`, `counts_sensitivity.csv` |
| `06_survival_analysis.R` | Kaplan-Meier by group; Cox regression of group, age and stage (stage I and II patients) | subtyped object | `Survival/Cox_regression_summary.csv`, plots |
| `07_GEO_validation.R` | the same scoring on GSE71729, cut-off over all 357 samples, survival in the 125 patients with data | GEO | `GEO_Validation/GEO_scores_and_subtypes.csv`, `GEO_summary_stats.csv`, plot |
| `compare_with_expected.R` | compares the tables with `expected_results/` | | |

Gene sets (identical in 02 and 07): 16 fucosylation genes (FUT1 to FUT11, GMDS, TSTA3, FUK,
FPGT, SLC35C1), 7 mannose-processing genes (MAN1A1, MAN1A2, MAN1C1, MAN2A1, MAN2A2, MGAT1,
MGAT2) and 9 mannose-biosynthesis genes (ALG3, ALG9, ALG12, DPM1, DPAGT1, MPI, PMM2, MRC1,
MRC2).

## Reference results

| | Result |
|---|---|
| TCGA-PAAD groups | 89 Fuc-favored, 89 Man-favored (178 primary tumours) |
| Differential expression | 450 up- and 523 down-regulated genes in Fuc-favored relative to Man-favored |
| GSEA (Hallmark) | 15 significant sets, for example E2F targets, MYC targets, G2M checkpoint (negative NES) and bile-acid and fatty-acid metabolism (positive) |
| TCGA survival (177 patients) | Kaplan-Meier log-rank P = 0.161 (median 21.1 versus 18.9 months, Fuc-favored versus Man-favored) |
| Cox, stage I and II | Fuc-low Man-high: HR 1.86, P = 0.016 (adjusted for age and stage) |
| Moffitt subtype (TPM) | glyco shift score against the Moffitt score: Spearman ρ = −0.57; Man-high tumours are more often basal-like by the classifier (23/89 versus 3/89; odds ratio 9.9, 95% CI 2.8 to 53.6, Fisher P = 2.5e-5); the Moffitt call separates survival (log-rank P = 0.020 by the classifier, 0.0079 by the score), the glycophenotype does not (P = 0.16) |
| GSE71729 validation | 357 samples scored, 125 with survival data; groups 81 versus 44 (cut-off over all samples); log-rank P = 0.101; median survival 19 versus 12 months |

Read the validation result with care: the median cut-off is taken over all 357 samples,
of which only the 145 primary tumours are tumour tissue and only 125 have survival data
(the other samples are normal pancreas, metastases and unlabelled). That is why the groups are
unbalanced (81 versus 44) and the result is not statistically significant.

## The Moffitt comparison

Moffitt et al. (Nat Genet 2015) defined the basal-like and classical subtypes by clustering the
tumour-specific part of the expression. A new cohort is scored against their marker genes
instead of being re-clustered:

* **Moffitt score**: ssGSEA of the 25 basal-like and the 25 classical genes, each z-scored over
  the cohort, and the difference z(basal-like) − z(classical). Higher means more basal-like. A
  call at 0 is arbitrary, so the continuous score is the quantity to report.
* **Classifier**: a penalised logistic model on 21 gene pairs (is the basal-like marker above the
  classical marker in this tumour?). It returns a probability for each tumour on its own.

Gene lists and classifier are the authors' own, from pdacR (`data/moffitt/`, MIT licence). Of the
50 genes, 46 are used: `CTSL2` is now `CTSV`; `ANXA8L2`, `ATAD4` and `LOC400573` are not in the
current annotation. Both readouts agree in direction. The input is TPM: the classifier compares
two different genes inside one tumour, so the result must not depend on gene length, which raw
counts do. The shipped matrix (`data/TCGA_PAAD_tpm.rds`, 22 MB) is the GDC TPM of the same 178
tumours with the same duplicate-symbol rule as the counts, and covers the whole transcriptome
because ssGSEA ranks every gene. Scoring the raw counts instead (`counts_sensitivity.csv`) leaves the association with the
glycophenotype unchanged (Spearman ρ −0.57 with the score, odds ratio 10.2), but 15 of the
classifier calls change and so do the survival P values; the TPM result is the one to quote. TCGA-PAAD tumours differ widely in purity and purity was not
accounted for. The result is an association, not an equivalence: about two thirds of the Man-high
tumours are called classical.

## Things that can differ between runs

* **GSEA** is permutation based and the scripts do not seed it, so scores change by a few
  hundredths and borderline sets can enter or leave the significant list.
* **KEGG** needs `rest.kegg.jp`; if it is unreachable the script says so and continues without the KEGG tables and figures.
* **Package versions.** The reference tables were produced with earlier Bioconductor
  releases; the comparison script was tested with the versions in `../docs/INSTALL.md` and
  every deterministic table agrees to within 1e-6 (relative).

## Changes relative to the exploratory scripts

All comments and messages are in English. No analysis logic was changed. Compatibility
fixes: `04` registers a serial back-end instead of passing `BPPARAM` to `GSEA()` and accepts
both the old and new `msigdbr` argument names; `00` lists exactly the packages that the
scripts load and fails if any is missing; the two versions of the GEO script were merged
into `07`, which now also writes the summary statistics that a separate helper script used to
print.

## Data

* `data/TCGA_PAAD_processed.rds`: counts (59,427 genes × 178 samples) and clinical data of
  TCGA-PAAD primary tumours, prepared by step 01 from the GDC release listed in
  `data/GDC_manifest.txt`.
* GSE71729 is downloaded from GEO by step 07 into `data/GEO_cache/`.
* TCGA and GEO data remain subject to the terms of the Genomic Data Commons and GEO.
