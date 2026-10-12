##############################################################################
# 01_data_download.R
# Download and pre-process the TCGA-PAAD data
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("data", showWarnings = FALSE)

message("Querying TCGA-PAAD transcriptome profiling data...")
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

message("Downloading data...")
GDCdownload(query, directory = "data/GDCdata")

message("Preparing data into SummarizedExperiment object...")
se <- GDCprepare(query, directory = "data/GDCdata")

# Extract the expression matrix (unstranded raw counts) and the clinical data
counts <- assay(se, "unstranded")
clin <- colData(se)

# Filter samples: keep primary solid tumours (TP) only
table(clin$sample_type)
keep_idx <- clin$sample_type == "Primary Tumor"
counts <- counts[, keep_idx]
clin <- clin[keep_idx, ]

# Convert ENSEMBL identifiers to gene symbols
gene_info <- rowData(se)
row.names(counts) <- gene_info$gene_name

# If a gene symbol is duplicated, keep the row with the highest total expression
sums <- rowSums(counts)
counts <- counts[order(sums, decreasing = TRUE), ]
counts <- counts[!duplicated(rownames(counts)), ]

message(sprintf("Processed data: %d genes across %d PAAD primary tumor samples", 
                nrow(counts), ncol(counts)))

saveRDS(list(counts = counts, clin = clin), file = "data/TCGA_PAAD_processed.rds")
message("=== Data download and preprocessing complete. Saved to data/TCGA_PAAD_processed.rds ===")
