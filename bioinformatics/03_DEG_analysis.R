##############################################################################
# 03_DEG_analysis.R
# Differential expression between the two groups
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

message("Loading subtyped data...")
dat <- readRDS("data/TCGA_PAAD_subtyped.rds")
counts <- dat$counts
clin <- dat$clin
scores <- dat$scores

# Make sure expression columns and clinical rows match
overlap_samples <- intersect(colnames(counts), clin$barcode)
counts <- counts[, overlap_samples]
clin <- clin[match(overlap_samples, clin$barcode), ]

# Integer counts for DESeq2
counts <- round(counts)

# Drop genes with very low expression
keep <- rowSums(counts >= 10) >= 10
counts <- counts[keep, ]

message("Running DESeq2...")
# Fuc-favored versus Man-favored groups
dds <- DESeqDataSetFromMatrix(
     countData = counts,
     colData = clin,
     design = ~Final_Subtype
)

# The reference group is Man-favored, so results are Fuc-favored relative to Man-favored
dds$Final_Subtype <- relevel(factor(dds$Final_Subtype), ref = "Man-favored")

dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# Define up- and down-regulated genes
res_df$Significance <- "ns"
res_df$Significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up"
res_df$Significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down"

# Save the full differential-expression table
dir.create("results/DEG", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/DEG/DESeq2_results.csv", row.names = FALSE)

# Volcano plot
message("Plotting Volcano Plot...")
# Label the top genes
top_up <- res_df %>%
     filter(Significance == "Up") %>%
     arrange(padj) %>%
     head(10)
top_down <- res_df %>%
     filter(Significance == "Down") %>%
     arrange(padj) %>%
     head(10)
label_genes <- rbind(top_up, top_down)

p <- ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj), colour = Significance)) +
     geom_point(alpha = 0.6, size = 1.5) +
     scale_colour_manual(values = c("Up" = "#d7301f", "Down" = "#0570b0", "ns" = "grey80")) +
     geom_vline(xintercept = c(-1, 1), lty = 2, col = "black", lwd = 0.5) +
     geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5) +
     theme_classic() +
     theme(axis.title = element_text(size = 14)) +
     labs(
          title = "",
          x = expression(log[2]~"Fold Change"),
          y = expression(-log[10](p[adj]))
     ) +
     geom_text_repel(
          data = label_genes, aes(label = Gene),
          size = 3,
          box.padding = 0.5,
          point.padding = 0.5,
          segment.colour = "black",
          show.legend = FALSE
     )

pdf("results/DEG/Volcano_Plot.pdf", width = 6, height = 4)
print(p)
dev.off()

svglite::svglite("results/DEG/Volcano_Plot.svg", width = 6, height = 4)
print(p)
dev.off()

message("=== DEG Analysis complete ===")
