# Moffitt reference data

Taken unchanged from the pdacR package of the Moffitt laboratory
(https://github.com/rmoffitt/pdacR, package version 0.1.2, repository commit `5b044b2f32`, MIT
licence, copyright 2016 Richard A. Moffitt). The three files below were downloaded on
2026-10-12; the classifier file is byte-identical to the one in that commit and the two gene
lists equal `Moffitt.Basal.25` and `Moffitt.Classical.25` there. The classifier file was last
changed in pdacR commit `2b91af4a0f` (2021-07-20) and the gene-list file in commit `302aa2f0a8`
(2022-09-07):

| File | Source in pdacR | Content |
|---|---|---|
| `moffitt_basal25.txt` | `data/gene_lists.rds`, `Moffitt.Basal.25` | the 25 basal-like tumour genes of Moffitt et al. 2015 |
| `moffitt_classical25.txt` | `data/gene_lists.rds`, `Moffitt.Classical.25` | the 25 classical tumour genes |
| `moffitt_classifier_21tsp.RData` | `data/Moffitt_classifier_2019.RData` | single-sample classifier: a penalised logistic model on 21 gene pairs |

Reference: Moffitt, R. A. et al. Virtual microdissection identifies distinct tumor- and
stroma-specific subtypes of pancreatic ductal adenocarcinoma. Nat. Genet. 47, 1168-1178 (2015).

The MIT licence requires its notice to travel with copies: "Permission is hereby granted, free
of charge, to any person obtaining a copy of this software and associated documentation files
..., subject to the following conditions: the above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software. The software is
provided 'as is', without warranty of any kind."

`../TCGA_PAAD_tpm.rds` (one level up) holds the TPM of all genes for the 178 TCGA-PAAD tumours,
taken from the GDC STAR output; `05_moffitt_scoring.R` says how it is rebuilt.
