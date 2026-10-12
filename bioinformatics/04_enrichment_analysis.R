##############################################################################
# 04_enrichment_analysis.R
# GO, KEGG and GSEA enrichment of the differentially expressed genes
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(DOSE)

message("Loading DEG results...")
res_df <- read.csv("results/DEG/DESeq2_results.csv")
res_df <- res_df[!is.na(res_df$padj), ]

# Extract the up-regulated (Fuc-favored) and down-regulated (Man-favored) gene lists
up_genes <- res_df$Gene[res_df$log2FoldChange > 1 & res_df$padj < 0.05]
down_genes <- res_df$Gene[res_df$log2FoldChange < -1 & res_df$padj < 0.05]

message(sprintf("Found %d Up-regulated (Fuc-favored) and %d Down-regulated (Man-favored) DEGs.", length(up_genes), length(down_genes)))

# Convert gene symbols to Entrez identifiers
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

dir.create("results/Enrichment", recursive = TRUE, showWarnings = FALSE)

# 1. GO analysis (biological process)
message("Running GO Enrichment...")
go_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

go_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Save the full enrichment tables (source data of the dot plots)
write.csv(as.data.frame(go_up), "results/Enrichment/GO_Enrichment_Up_FucHigh.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "results/Enrichment/GO_Enrichment_Down_ManHigh.csv", row.names = FALSE)

pdf("results/Enrichment/GO_Enrichment_Up_FucHigh.pdf", width = 8, height = 6)
print(dotplot(go_up, showCategory = 15, title = expression("GO BP Enrichment: " * plain(Fuc)^plain(high) ~ plain(Man)^plain(low) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
dev.off()

svglite::svglite("results/Enrichment/GO_Enrichment_Up_FucHigh.svg", width = 8, height = 6)
print(dotplot(go_up, showCategory = 15, title = expression("GO BP Enrichment: " * plain(Fuc)^plain(high) ~ plain(Man)^plain(low) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
dev.off()

pdf("results/Enrichment/GO_Enrichment_Down_ManHigh.pdf", width = 8, height = 6)
print(dotplot(go_down, showCategory = 15, title = expression("GO BP Enrichment: " * plain(Fuc)^plain(low) ~ plain(Man)^plain(high) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
dev.off()

svglite::svglite("results/Enrichment/GO_Enrichment_Down_ManHigh.svg", width = 8, height = 6)
print(dotplot(go_down, showCategory = 15, title = expression("GO BP Enrichment: " * plain(Fuc)^plain(low) ~ plain(Man)^plain(high) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
dev.off()

# 2. KEGG analysis
message("Running KEGG Enrichment...")
tryCatch(
  {
    kegg_up <- enrichKEGG(gene = up_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
    kegg_down <- enrichKEGG(gene = down_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

    write.csv(as.data.frame(kegg_up), "results/Enrichment/KEGG_Enrichment_Up_FucHigh.csv", row.names = FALSE)
    write.csv(as.data.frame(kegg_down), "results/Enrichment/KEGG_Enrichment_Down_ManHigh.csv", row.names = FALSE)

    pdf("results/Enrichment/KEGG_Enrichment_Up_FucHigh.pdf", width = 8, height = 6)
    print(dotplot(kegg_up, showCategory = 15, title = expression("KEGG Enrichment: " * plain(Fuc)^plain(high) ~ plain(Man)^plain(low) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
    dev.off()

    svglite::svglite("results/Enrichment/KEGG_Enrichment_Up_FucHigh.svg", width = 8, height = 6)
    print(dotplot(kegg_up, showCategory = 15, title = expression("KEGG Enrichment: " * plain(Fuc)^plain(high) ~ plain(Man)^plain(low) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
    dev.off()

    pdf("results/Enrichment/KEGG_Enrichment_Down_ManHigh.pdf", width = 8, height = 6)
    print(dotplot(kegg_down, showCategory = 15, title = expression("KEGG Enrichment: " * plain(Fuc)^plain(low) ~ plain(Man)^plain(high) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
    dev.off()

    svglite::svglite("results/Enrichment/KEGG_Enrichment_Down_ManHigh.svg", width = 8, height = 6)
    print(dotplot(kegg_down, showCategory = 15, title = expression("KEGG Enrichment: " * plain(Fuc)^plain(low) ~ plain(Man)^plain(high) * " glycophenotype")) + theme(axis.title = element_text(size = 14)))
    dev.off()
  },
  error = function(e) {
    message("KEGG enrichment failed, likely due to a network connection timeout to rest.kegg.jp.")
    message("Caught error: ", e$message)
  }
)

# 3. GSEA with the Hallmark gene sets
message("Running GSEA Analysis...")
library(msigdbr)
# msigdbr >= 10 renamed the `category` argument to `collection`
if (packageVersion("msigdbr") >= "10.0.0") {
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
} else {
  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
}
m_t2g <- hallmark %>%
  dplyr::select(gs_name, gene_symbol)

# Ranked list (descending Wald statistic, more robust than log2 fold change alone)
gene_list <- res_df$stat
names(gene_list) <- res_df$Gene
gene_list <- sort(gene_list, decreasing = TRUE)

# Serial back-end, registered globally: not every clusterProfiler version accepts a BPPARAM argument
BiocParallel::register(BiocParallel::SerialParam())
gsea_res <- GSEA(gene_list, TERM2GENE = m_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH")

write.csv(as.data.frame(gsea_res), "results/Enrichment/GSEA_Hallmark_results.csv")

pdf("results/Enrichment/GSEA_Ridgeplot.pdf", width = 12, height = 10)
print(ridgeplot(gsea_res) + labs(title = "") + theme(axis.title = element_text(size = 14)))
dev.off()

svglite::svglite("results/Enrichment/GSEA_Ridgeplot.svg", width = 12, height = 10)
print(ridgeplot(gsea_res) + labs(title = "") + theme(axis.title = element_text(size = 14)))
dev.off()

message("=== Enrichment Analysis complete ===")
