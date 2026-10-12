##############################################################################
# 05_moffitt_scoring.R
# Relate the fucose / mannose glycophenotype to the Moffitt basal-like and classical subtypes
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it. Needs the output of
# 02_ssGSEA_scoring.R.
#
# Moffitt et al. (Nat Genet 2015) defined the two subtypes by clustering the tumour-specific
# part of the expression. A new cohort is not re-clustered; each tumour is scored against the
# 25 basal-like and 25 classical tumour genes in two ways:
#   * Moffitt score: ssGSEA scores of the two signatures, each z-scored over the cohort, and the
#     difference z(basal-like) - z(classical). Higher means more basal-like.
#   * Classifier probability: a penalised logistic model on 21 gene pairs ("is the basal-like
#     marker expressed above the classical marker in this tumour?"). It uses only the ranks of
#     genes within one tumour, so it needs no cut-off taken from the cohort.
#
# The input is TPM. The classifier compares the expression of two different genes inside one
# tumour; raw counts depend on gene length, so a pair of genes of very different length could be
# ordered by length instead of by expression. TPM corrects for that. The counts-based result is
# kept as a sensitivity analysis (results/Moffitt/counts_sensitivity.csv).
#
# The gene lists and the classifier are the authors' own, from the pdacR package (data/moffitt/).

library(GSVA)
library(survival)

dir.create("results/Moffitt", recursive = TRUE, showWarnings = FALSE)
tpm_file <- "data/TCGA_PAAD_tpm.rds"
glyco_file <- "results/clustering/ssGSEA_scores_and_subtypes.csv"
if (!file.exists(glyco_file)) stop("Run 02_ssGSEA_scoring.R first: ", glyco_file, " is missing.")

# ---- functions --------------------------------------------------------------------------------------

# Symbols of the 2015 lists that have been renamed since: CTSL2 is now CTSV.
renamed_genes <- c(CTSL2 = "CTSV")

# The two signatures, mapped onto the symbols of an expression matrix, and a table of what
# happened to every gene.
load_signatures <- function(genes_available, dir = "data/moffitt") {
  lists <- list(
    "basal-like" = readLines(file.path(dir, "moffitt_basal25.txt")),
    "classical" = readLines(file.path(dir, "moffitt_classical25.txt"))
  )
  coverage <- do.call(rbind, lapply(names(lists), function(s) data.frame(signature = s, gene = lists[[s]])))
  new_symbol <- ifelse(coverage$gene %in% names(renamed_genes), renamed_genes[coverage$gene], coverage$gene)
  coverage$used_as <- ifelse(new_symbol %in% genes_available, new_symbol, NA)
  coverage$status <- ifelse(is.na(coverage$used_as), "not in annotation",
    ifelse(coverage$used_as == coverage$gene, "present", "present under new symbol")
  )
  list(
    sets = list(
      Moffitt_Basal = na.omit(coverage$used_as[coverage$signature == "basal-like"]),
      Moffitt_Classical = na.omit(coverage$used_as[coverage$signature == "classical"])
    ),
    coverage = coverage
  )
}

# ssGSEA scores of the two signatures and the Moffitt score, z(basal-like) - z(classical).
moffitt_scores <- function(log_expression, sets) {
  s <- as.data.frame(t(gsva(ssgseaParam(log_expression, sets))))
  s$SampleID <- rownames(s)
  s$moffitt_score <- as.numeric(scale(s$Moffitt_Basal)) - as.numeric(scale(s$Moffitt_Classical))
  s
}

# Probability of the basal-like subtype from the 21-gene-pair classifier.
classifier_probability <- function(log_expression, file = "data/moffitt/moffitt_classifier_21tsp.RData") {
  e <- new.env()
  load(file, envir = e)
  cl <- e$Moffitt_classifier_2019
  stopifnot(all(unique(c(cl$TSPs)) %in% rownames(log_expression)))
  # One indicator per pair: 1 when the first gene (a basal-like marker) is above the second
  # (a classical marker) in the same tumour.
  indicator <- sapply(seq_len(nrow(cl$TSPs)), function(i) {
    as.numeric(log_expression[cl$TSPs[i, 1], ] > log_expression[cl$TSPs[i, 2], ])
  })
  # The model was fitted on standardised indicators: apply the training centre and scale.
  x_train <- cl$fit$X
  z <- sweep(sweep(indicator, 2, attr(x_train, "center")), 2, attr(x_train, "scale"), "/")
  as.numeric(1 / (1 + exp(-(cl$fit$beta[1, 1] + z %*% cl$fit$beta[-1, 1]))))
}

# Score one expression matrix: returns one row per tumour with the glyco scores attached.
score_cohort <- function(log_expression, glyco) {
  sig <- load_signatures(rownames(log_expression))
  s <- moffitt_scores(log_expression, sig$sets)
  s$classifier_basal_prob <- classifier_probability(log_expression)
  m <- merge(s, glyco, by = "SampleID")
  m$glycophenotype <- factor(m$Final_Subtype,
    levels = c("Fuc-favored", "Man-favored"), labels = c("FucHigh_ManLow", "FucLow_ManHigh")
  )
  m$call_score <- ifelse(m$moffitt_score > 0, "basal-like", "classical")
  m$call_classifier <- ifelse(m$classifier_basal_prob > 0.5, "basal-like", "classical")
  list(data = m, coverage = sig$coverage)
}

# Comparison with the glycophenotype scores of 02_ssGSEA_scoring.R.
compare_with_glyco <- function(m) {
  readouts <- c("moffitt_score", "Moffitt_Basal", "Moffitt_Classical", "classifier_basal_prob")
  correlations <- do.call(rbind, lapply(c("Glyco_Shift_Score", "Fuc_Z", "Man_Comp_Z"), function(v) {
    do.call(rbind, lapply(readouts, function(w) {
      r <- cor.test(m[[v]], m[[w]], method = "spearman", exact = FALSE)
      data.frame(glyco_variable = v, moffitt_variable = w, spearman_rho = unname(r$estimate), p_value = r$p.value)
    }))
  }))
  by_group <- do.call(rbind, lapply(readouts, function(w) {
    a <- m[[w]][m$glycophenotype == "FucHigh_ManLow"]
    b <- m[[w]][m$glycophenotype == "FucLow_ManHigh"]
    wt <- wilcox.test(a, b, exact = FALSE) # W counts the pairs in which a Fuc-high value exceeds a Man-high value
    data.frame(
      variable = w, n_FucHigh = length(a), n_ManHigh = length(b),
      median_FucHigh = median(a), median_ManHigh = median(b),
      # probability that a random Man-high tumour scores above a random Fuc-high tumour
      auc_ManHigh_above_FucHigh = mean(outer(b, a, ">")) + 0.5 * mean(outer(b, a, "==")),
      wilcoxon_W = unname(wt$statistic), wilcoxon_p = wt$p.value
    )
  }))
  calls <- do.call(rbind, lapply(c("call_score", "call_classifier"), function(cn) {
    tab <- table(glycophenotype = m$glycophenotype, call = m[[cn]])
    ft <- fisher.test(tab)
    data.frame(
      call = cn,
      FucHigh_basal = tab["FucHigh_ManLow", "basal-like"], FucHigh_classical = tab["FucHigh_ManLow", "classical"],
      ManHigh_basal = tab["FucLow_ManHigh", "basal-like"], ManHigh_classical = tab["FucLow_ManHigh", "classical"],
      share_basal_FucHigh = tab["FucHigh_ManLow", "basal-like"] / sum(tab["FucHigh_ManLow", ]),
      share_basal_ManHigh = tab["FucLow_ManHigh", "basal-like"] / sum(tab["FucLow_ManHigh", ]),
      # odds of a basal-like call in Man-high relative to Fuc-high
      odds_ratio_ManHigh_vs_FucHigh = 1 / unname(ft$estimate),
      or_ci_low = 1 / ft$conf.int[2], or_ci_high = 1 / ft$conf.int[1], fisher_p = ft$p.value
    )
  }))
  list(correlations = correlations, by_group = by_group, calls = calls)
}

# Overall survival by Moffitt call and by glycophenotype (log-rank), and the hazard ratio of the
# Moffitt score per standard deviation.
survival_summary <- function(m, clin) {
  clin$OS.time <- ifelse(is.na(clin$days_to_death), clin$days_to_last_follow_up, clin$days_to_death)
  clin$OS <- ifelse(clin$vital_status == "Dead", 1, 0)
  sv <- merge(m, clin[, c("barcode", "OS.time", "OS")], by.x = "SampleID", by.y = "barcode")
  sv <- sv[!is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS), ]
  sv$months <- sv$OS.time / 30
  # the hazard ratio is that of the second level: basal-like versus classical, Fuc-low Man-high versus Fuc-high Man-low
  sv$call_score <- factor(sv$call_score, levels = c("classical", "basal-like"))
  sv$call_classifier <- factor(sv$call_classifier, levels = c("classical", "basal-like"))
  logrank <- function(label, f) {
    lr <- survdiff(f, data = sv)
    hr <- summary(coxph(f, data = sv))$conf.int
    data.frame(
      comparison = label, n = nrow(sv), events = sum(lr$obs),
      logrank_chisq = lr$chisq, logrank_df = length(lr$n) - 1, logrank_p = 1 - pchisq(lr$chisq, length(lr$n) - 1),
      hr = hr[1, "exp(coef)"], hr_ci_low = hr[1, "lower .95"], hr_ci_high = hr[1, "upper .95"]
    )
  }
  out <- rbind(
    logrank("Moffitt call, from the score", Surv(months, OS) ~ call_score),
    logrank("Moffitt call, from the classifier", Surv(months, OS) ~ call_classifier),
    logrank("glycophenotype", Surv(months, OS) ~ glycophenotype)
  )
  cox <- summary(coxph(Surv(months, OS) ~ scale(moffitt_score), data = sv))$coefficients
  out$hr_per_sd_of_score <- c(exp(cox[1, "coef"]), NA, NA)
  out$hr_p <- c(cox[1, "Pr(>|z|)"], NA, NA)
  out
}

# ---- data -------------------------------------------------------------------------------------------
dat <- readRDS("data/TCGA_PAAD_processed.rds")
clin <- as.data.frame(dat$clin)
glyco <- read.csv(glyco_file)
barcodes <- colnames(dat$counts)
stopifnot(setequal(glyco$SampleID, barcodes))

# The TPM matrix of the 178 tumours covers the whole transcriptome (ssGSEA ranks every gene). It
# is shipped; if it is missing it is rebuilt from the GDC STAR files that 01_data_download.R
# downloaded (data/GDCdata, or GDC_DIR). The rebuild needs TCGAbiolinks and one query to the GDC,
# which only asks which file belongs to which sample.
if (!file.exists(tpm_file)) {
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  gdc_dir <- Sys.getenv("GDC_DIR", "data/GDCdata")
  if (!dir.exists(gdc_dir)) stop("Neither ", tpm_file, " nor the GDC files in ", gdc_dir, " are available.")
  se <- GDCprepare(GDCquery(
    project = "TCGA-PAAD", data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts"
  ), directory = gdc_dir)
  tpm <- assay(se, "tpm_unstrand")[, barcodes]
  rownames(tpm) <- rowData(se)$gene_name
  # same rule as for the counts: for a duplicated symbol keep the row with the highest total
  tpm <- tpm[order(rowSums(tpm), decreasing = TRUE), ]
  tpm <- tpm[!duplicated(rownames(tpm)), ]
  saveRDS(tpm, tpm_file)
}
tpm <- readRDS(tpm_file)
stopifnot(setequal(colnames(tpm), barcodes))

# ---- scoring on TPM (the result to report) ----------------------------------------------------------------
tpm_result <- score_cohort(log2(tpm[, barcodes] + 1), glyco)
m <- tpm_result$data
write.csv(tpm_result$coverage, "results/Moffitt/gene_coverage.csv", row.names = FALSE)
print(table(tpm_result$coverage$signature, tpm_result$coverage$status))
cat(sprintf(
  "\nClassifier probability versus Moffitt score: Spearman rho = %.3f (positive: the orientation is right)\n",
  cor(m$classifier_basal_prob, m$moffitt_score, method = "spearman")
))

res <- compare_with_glyco(m)
surv <- survival_summary(m, clin)
cat("\nSpearman correlation with the glyco shift score:\n")
print(subset(res$correlations, glyco_variable == "Glyco_Shift_Score"), digits = 3, row.names = FALSE)
cat("\nScores by glycophenotype:\n")
print(res$by_group, digits = 3, row.names = FALSE)
cat("\nCalls by glycophenotype:\n")
print(res$calls, digits = 3, row.names = FALSE)
cat("\nSurvival:\n")
print(surv, digits = 3, row.names = FALSE)

write.csv(res$correlations, "results/Moffitt/correlations_with_glyco_scores.csv", row.names = FALSE)
write.csv(res$by_group, "results/Moffitt/scores_by_glycophenotype.csv", row.names = FALSE)
write.csv(res$calls, "results/Moffitt/calls_by_glycophenotype.csv", row.names = FALSE)
write.csv(surv, "results/Moffitt/survival.csv", row.names = FALSE)
write.csv(m[, c(
  "SampleID", "Moffitt_Basal", "Moffitt_Classical", "moffitt_score", "classifier_basal_prob",
  "call_score", "call_classifier", "glycophenotype", "Glyco_Shift_Score"
)], "results/Moffitt/scores_per_sample.csv", row.names = FALSE)

# ---- sensitivity: the same on raw counts --------------------------------------------------------------------
counts_result <- score_cohort(log2(dat$counts[, barcodes] + 1), glyco)
mc <- counts_result$data
rc <- compare_with_glyco(mc)
sc <- survival_summary(mc, clin)
summarise <- function(input, m, r, s) {
  shift <- subset(r$correlations, glyco_variable == "Glyco_Shift_Score")
  call <- subset(r$calls, call == "call_classifier")
  data.frame(
    input = input,
    basal_calls_by_classifier = sum(m$call_classifier == "basal-like"),
    spearman_shift_vs_score = shift$spearman_rho[shift$moffitt_variable == "moffitt_score"],
    spearman_shift_vs_probability = shift$spearman_rho[shift$moffitt_variable == "classifier_basal_prob"],
    odds_ratio_classifier_call = call$odds_ratio_ManHigh_vs_FucHigh,
    fisher_p_classifier_call = call$fisher_p,
    logrank_p_score_call = s$logrank_p[s$comparison == "Moffitt call, from the score"],
    logrank_p_classifier_call = s$logrank_p[s$comparison == "Moffitt call, from the classifier"]
  )
}
sensitivity <- rbind(summarise("TPM (reported)", m, res, surv), summarise("raw counts", mc, rc, sc))
both <- merge(m, mc[, c("SampleID", "moffitt_score", "classifier_basal_prob", "call_classifier")],
  by = "SampleID", suffixes = c("_tpm", "_counts")
)
sensitivity$spearman_score_tpm_vs_counts <- c(cor(both$moffitt_score_tpm, both$moffitt_score_counts, method = "spearman"), NA)
sensitivity$spearman_probability_tpm_vs_counts <- c(cor(both$classifier_basal_prob_tpm, both$classifier_basal_prob_counts, method = "spearman"), NA)
sensitivity$classifier_calls_that_differ <- c(sum(both$call_classifier_tpm != both$call_classifier_counts), NA)
write.csv(sensitivity, "results/Moffitt/counts_sensitivity.csv", row.names = FALSE)
cat("\nSensitivity, TPM versus raw counts:\n")
print(t(sensitivity), quote = FALSE)
message("=== Moffitt scoring complete ===")
