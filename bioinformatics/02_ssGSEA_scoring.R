##############################################################################
# 02_ssGSEA_scoring.R
# Compute ssGSEA scores and stratify the patients by glycophenotype
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.

library(GSVA)
library(pheatmap)
library(dplyr)
library(ggplot2)

message("Loading processed data...")
dat <- readRDS("data/TCGA_PAAD_processed.rds")
counts <- dat$counts
clin <- dat$clin

# TPM conversion or log2(counts + 1)
# For ssGSEA ranking, either log2(CPM + 1) or simply log2(counts + 1) would do
# log2(counts + 1) is used here
log2_counts <- log2(counts + 1)

# Define the gene sets (signatures)
fuc_sig <- c(
  "FUT1", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT7",
  "FUT8", "FUT9", "FUT10", "FUT11", "GMDS", "TSTA3", "FUK", "FPGT", "SLC35C1"
)

man_processing <- c("MAN1A1", "MAN1A2", "MAN1C1", "MAN2A1", "MAN2A2", "MGAT1", "MGAT2")

man_biosynthesis <- c("ALG3", "ALG9", "ALG12", "DPM1", "DPAGT1", "MPI", "PMM2", "MRC1", "MRC2")

# Keep only genes present in the expression matrix
genesets <- list(
  Fuc_Signature = intersect(fuc_sig, rownames(log2_counts)),
  Man_Processing = intersect(man_processing, rownames(log2_counts)),
  Man_Biosynthesis = intersect(man_biosynthesis, rownames(log2_counts))
)

message("Calculating ssGSEA scores...")
param <- ssgseaParam(log2_counts, genesets)
gsva_res <- gsva(param)

# Extract the scores and compute the composite mannose score
scores_df <- as.data.frame(t(gsva_res))
scores_df$SampleID <- rownames(scores_df)

# Composite mannose score based on the biological phenotype
# high-mannose phenotype = high mannose biosynthesis - low mannose processing
scores_df$Man_Composite <- scores_df$Man_Biosynthesis - scores_df$Man_Processing

message("Calculating Glycan Shift Score and stratifying patients...")

# Glyco_Shift_Score: shift towards fucosylation (Fuc) versus mannosylation (Man)
# Z-score both signatures first so that they share one scale
scores_df$Fuc_Z <- scale(scores_df$Fuc_Signature)[, 1]
scores_df$Man_Comp_Z <- scale(scores_df$Man_Composite)[, 1]

scores_df$Glyco_Shift_Score <- scores_df$Fuc_Z - scores_df$Man_Comp_Z

# Split into two halves at the median
median_score <- median(scores_df$Glyco_Shift_Score, na.rm = TRUE)
scores_df$Final_Subtype <- ifelse(scores_df$Glyco_Shift_Score >= median_score, "Fuc-favored", "Man-favored")

dir.create("results/clustering", recursive = TRUE, showWarnings = FALSE)

# Order the samples by score so that the heatmap reads smoothly
scores_df <- scores_df[order(scores_df$Glyco_Shift_Score, decreasing = TRUE), ]
cluster_data <- t(scale(scores_df[, c("Fuc_Signature", "Man_Composite")]))

# Heatmap
annotation_col <- data.frame(glycophenotype = factor(
  scores_df$Final_Subtype,
  levels = c("Fuc-favored", "Man-favored"),
  labels = c("Fuc-high Man-low", "Fuc-low Man-high")
))
rownames(annotation_col) <- scores_df$SampleID
annotation_colours <- list(glycophenotype = c("Fuc-high Man-low" = "#d7301f", "Fuc-low Man-high" = "#0570b0"))

pdf("results/clustering/ssGSEA_heatmap.pdf", width = 7, height = 4)
pheatmap(cluster_data,
  annotation_col = annotation_col,
  annotation_colors = annotation_colours,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  main = ""
)
dev.off()

svglite::svglite("results/clustering/ssGSEA_heatmap.svg", width = 7, height = 4)
pheatmap(cluster_data,
  annotation_col = annotation_col,
  annotation_colors = annotation_colours,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  main = ""
)
dev.off()

# Write the group assignment back to the clinical table and save
clin$Final_Subtype <- scores_df$Final_Subtype[match(clin$barcode, scores_df$SampleID)]

saveRDS(list(counts = counts, clin = clin, scores = scores_df),
  file = "data/TCGA_PAAD_subtyped.rds"
)

write.csv(scores_df, "results/clustering/ssGSEA_scores_and_subtypes.csv", row.names = FALSE)
message("=== Subtyping complete! View results/clustering for heatmap. ===")
