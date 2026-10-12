##############################################################################
# compare_with_expected.R
# Compare the tables written by run_all.sh with the reference tables in
# expected_results/
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.
#
# Deterministic tables must agree to within a small numerical tolerance. The GSEA table
# comes from permutation-based enrichment, so it is compared by pathway and by normalised
# enrichment score (NES) instead.

tolerance <- 1e-6

compare_table <- function(file, key = NULL) {
  new <- read.csv(file.path("results", file), check.names = FALSE)
  ref <- read.csv(file.path("expected_results", file), check.names = FALSE)
  if (is.null(key)) key <- names(ref)[1] # default: the first column identifies the rows
  if (!identical(dim(new), dim(ref)) || !setequal(new[[key]], ref[[key]])) {
    cat(sprintf("FAIL  %-48s different rows or columns\n", file))
    return(FALSE)
  }
  new <- new[match(ref[[key]], new[[key]]), ]
  worst <- 0
  same_text <- TRUE
  for (column in names(ref)) {
    if (is.numeric(ref[[column]])) {
      worst <- max(worst, abs(new[[column]] - ref[[column]]) / pmax(1, abs(ref[[column]])), na.rm = TRUE)
    } else if (column == "geneID") {
      # the genes of an enriched term, "/"-separated: their order can differ, the set cannot
      same_genes <- mapply(function(x, y) setequal(strsplit(x, "/")[[1]], strsplit(y, "/")[[1]]),
        as.character(new[[column]]), as.character(ref[[column]]))
      same_text <- same_text && all(same_genes)
    } else {
      same_text <- same_text && identical(as.character(new[[column]]), as.character(ref[[column]]))
    }
  }
  ok <- same_text && worst <= tolerance
  cat(sprintf("%s  %-48s max relative difference %.2g\n", if (ok) "ok  " else "FAIL", file, worst))
  ok
}

results <- c(
  compare_table("clustering/ssGSEA_scores_and_subtypes.csv", "SampleID"),
  compare_table("DEG/DESeq2_results.csv", "Gene"),
  compare_table("Survival/Cox_regression_summary.csv"),
  compare_table("GEO_Validation/GEO_scores_and_subtypes.csv", "SampleID"),
  compare_table("GEO_Validation/GEO_summary_stats.csv", "version"),
  compare_table("Enrichment/GO_Enrichment_Up_FucHigh.csv", "ID"),
  compare_table("Enrichment/GO_Enrichment_Down_ManHigh.csv", "ID"),
  compare_table("Moffitt/scores_per_sample.csv", "SampleID"),
  compare_table("Moffitt/scores_by_glycophenotype.csv", "variable"),
  compare_table("Moffitt/calls_by_glycophenotype.csv", "call"),
  compare_table("Moffitt/survival.csv", "comparison"),
  compare_table("Moffitt/counts_sensitivity.csv", "input")
)

# GSEA: every reference pathway must be recovered with a similar NES
gsea_new <- read.csv("results/Enrichment/GSEA_Hallmark_results.csv")
gsea_ref <- read.csv("expected_results/Enrichment/GSEA_Hallmark_results.csv")
shared <- merge(gsea_new[, c("ID", "NES")], gsea_ref[, c("ID", "NES")], by = "ID", suffixes = c(".new", ".ref"))
recovered <- all(gsea_ref$ID %in% gsea_new$ID)
worst_nes <- max(abs(shared$NES.new - shared$NES.ref))
extra <- setdiff(gsea_new$ID, gsea_ref$ID)
ok_gsea <- recovered && worst_nes < 0.1
cat(sprintf("%s  %-48s %d of %d reference pathways recovered, max |NES difference| %.3f%s\n",
            if (ok_gsea) "ok  " else "FAIL", "Enrichment/GSEA_Hallmark_results.csv",
            sum(gsea_ref$ID %in% gsea_new$ID), nrow(gsea_ref), worst_nes,
            if (length(extra)) paste0("; additional borderline sets: ", paste(extra, collapse = ", ")) else ""))
results <- c(results, ok_gsea)

# KEGG reads its annotation from a live web service, so only the recovery of the reference
# pathways is checked, and a missing table (server unreachable) is reported but not fatal.
for (file in c("Enrichment/KEGG_Enrichment_Up_FucHigh.csv", "Enrichment/KEGG_Enrichment_Down_ManHigh.csv")) {
  ref <- read.csv(file.path("expected_results", file))
  if (!file.exists(file.path("results", file))) {
    cat(sprintf("skip  %-48s not produced (KEGG server unreachable?)\n", file))
    next
  }
  new <- read.csv(file.path("results", file))
  recovered <- mean(ref$ID %in% new$ID)
  cat(sprintf("%s  %-48s %.0f%% of the %d reference pathways recovered\n",
              if (recovered == 1) "ok  " else "info", file, 100 * recovered, nrow(ref)))
}

cat(if (all(results)) "\nALL TABLES AGREE WITH THE REFERENCE\n" else "\nSOME TABLES DIFFER, see FAIL lines above\n")
quit(status = if (all(results)) 0 else 1)
