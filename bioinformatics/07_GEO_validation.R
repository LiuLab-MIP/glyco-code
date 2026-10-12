##############################################################################
# 07_GEO_validation.R
# Independent validation on a GEO data set (GSE71729)
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.
#
# All 357 samples of GSE71729 are scored, and the median cut-off of the
# Glyco_Shift_Score is taken over all of them (145 primary tumours, 134 normal
# pancreas, 61 metastases, 17 unlabelled). Survival data exist for 125 primary
# tumours only, so the Kaplan-Meier analysis uses those 125 patients.

library(GEOquery)
library(GSVA)
library(survival)
library(survminer)
library(dplyr)
library(pheatmap)

message("Downloading GSE71729 from GEO...")
options(timeout = 3000)
dir.create("data/GEO_cache", recursive = TRUE, showWarnings = FALSE)
gset <- getGEO("GSE71729", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/GEO_cache")
if (length(gset) > 1) {
  idx <- grep("GPL", attr(gset, "names"))[1]
  gset <- gset[[idx]]
} else {
  gset <- gset[[1]]
}

# Expression matrix and platform annotation
exprs_mat <- exprs(gset)
clin_geo <- pData(gset)

# No filtering by tissue type: every sample is retained for scoring
message("Full GSE71729 composition (all samples retained for scoring):")
print(table(clin_geo[["tissue type:ch2"]], useNA = "ifany"))

# The GPL annotation maps probes to gene symbols
fData_gset <- fData(gset)
probe_to_gene <- fData_gset$ORF

# A probe may map to several genes: take the first
probe_to_gene <- sapply(strsplit(as.character(probe_to_gene), "///"), `[`, 1)

# Use gene symbols as the row names of the expression matrix
rownames(exprs_mat) <- probe_to_gene

# Drop probes without a gene symbol
exprs_mat <- exprs_mat[!is.na(rownames(exprs_mat)) & rownames(exprs_mat) != "", ]

# De-duplicate probes: for each gene keep the probe with the highest total expression
sums <- rowSums(exprs_mat)
exprs_mat <- exprs_mat[order(sums, decreasing = TRUE), ]
exprs_mat <- exprs_mat[!duplicated(rownames(exprs_mat)), ]

message(sprintf(
  "Processing GSE71729 (all samples): %d genes across %d samples",
  nrow(exprs_mat), ncol(exprs_mat)
))

# Clinical data of GSE71729: survival time (months) and status
clin_geo$OS.time <- as.numeric(clin_geo$`survival_months:ch2`)
clin_geo$OS <- as.numeric(clin_geo$`death_event_1death_0censor:ch2`)

# Gene sets (identical to 02_ssGSEA_scoring.R)
fuc_sig <- c(
  "FUT1", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT7",
  "FUT8", "FUT9", "FUT10", "FUT11", "GMDS", "TSTA3", "FUK", "FPGT", "SLC35C1"
)
man_processing <- c("MAN1A1", "MAN1A2", "MAN1C1", "MAN2A1", "MAN2A2", "MGAT1", "MGAT2")
man_biosynthesis <- c("ALG3", "ALG9", "ALG12", "DPM1", "DPAGT1", "MPI", "PMM2", "MRC1", "MRC2")

genesets_geo <- list(
  Fuc_Signature = intersect(fuc_sig, rownames(exprs_mat)),
  Man_Processing = intersect(man_processing, rownames(exprs_mat)),
  Man_Biosynthesis = intersect(man_biosynthesis, rownames(exprs_mat))
)

message("Calculating ssGSEA scores on the validation cohort (all samples)...")
param_geo <- ssgseaParam(exprs_mat, genesets_geo)
gsva_res_geo <- gsva(param_geo)
scores_geo <- as.data.frame(t(gsva_res_geo))
scores_geo$SampleID <- rownames(scores_geo)
scores_geo$Man_Composite <- scores_geo$Man_Biosynthesis - scores_geo$Man_Processing

# Median split of the Glyco_Shift_Score (over all 357 samples), in place of unsupervised clustering
scores_geo$Fuc_Z <- scale(scores_geo$Fuc_Signature)[, 1]
scores_geo$Man_Comp_Z <- scale(scores_geo$Man_Composite)[, 1]
scores_geo$Glyco_Shift_Score <- scores_geo$Fuc_Z - scores_geo$Man_Comp_Z

median_score_geo <- median(scores_geo$Glyco_Shift_Score, na.rm = TRUE)
scores_geo$Final_Subtype <- ifelse(scores_geo$Glyco_Shift_Score >= median_score_geo, "Fuc-favored", "Man-favored")

dir.create("results/GEO_Validation", recursive = TRUE, showWarnings = FALSE)
write.csv(scores_geo, "results/GEO_Validation/GEO_scores_and_subtypes.csv", row.names = FALSE)

# Survival analysis
clin_geo$Final_Subtype <- scores_geo$Final_Subtype[match(rownames(clin_geo), scores_geo$SampleID)]
clin_geo_surv <- clin_geo %>% filter(!is.na(OS.time) & !is.na(OS))
message(sprintf("Survival cohort: %d patients", nrow(clin_geo_surv)))

fit_geo <- survfit(Surv(OS.time, OS) ~ Final_Subtype, data = clin_geo_surv)

p_geo <- ggsurvplot(fit_geo,
  data = clin_geo_surv,
  pval = TRUE, conf.int = TRUE, risk.table = FALSE,
  title = "",
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  font.x = 14,
  font.y = 14,
  font.legend = 12,
  font.tickslab = 12,
  palette = c("#d7301f", "#0570b0"),
  legend.title = "Glycophenotype",
  legend.labs = c("FuchighManlow", "FuclowManhigh"),
  short.panel.labs = TRUE
)

# Custom legend labels with plotmath (upright superscripts for high / low)
labels_expr <- c(
  expression(plain(Fuc)^plain(high) ~ plain(Man)^plain(low)),
  expression(plain(Fuc)^plain(low) ~ plain(Man)^plain(high))
)
p_geo$plot <- p_geo$plot +
  scale_colour_manual(values = c("#d7301f", "#0570b0"), labels = labels_expr) +
  scale_fill_manual(values = c("#d7301f", "#0570b0"), labels = labels_expr)

pdf("results/GEO_Validation/GEO_KM_Survival.pdf", width = 6, height = 4, onefile = FALSE)
print(p_geo)
dev.off()

svglite::svglite("results/GEO_Validation/GEO_KM_Survival.svg", width = 6, height = 4)
print(p_geo)
dev.off()

# Summary statistics: log-rank test and median survival per group
sd_geo <- survdiff(Surv(OS.time, OS) ~ Final_Subtype, data = clin_geo_surv)
pval <- 1 - pchisq(sd_geo$chisq, length(sd_geo$n) - 1)
cox_geo <- summary(coxph(Surv(OS.time, OS) ~ Final_Subtype, data = clin_geo_surv))$conf.int # hazard ratio of the second level (Man-favored)
med <- summary(fit_geo)$table[, "median"]
stats <- data.frame(
  version = "all_samples",
  n_scored = nrow(scores_geo),
  n_surv = nrow(clin_geo_surv),
  events = sum(sd_geo$obs),
  median_cut = median_score_geo,
  logrank_chisq = sd_geo$chisq,
  logrank_df = length(sd_geo$n) - 1,
  logrank_p = pval,
  hr_man_favored = cox_geo[1, "exp(coef)"],
  hr_ci_low = cox_geo[1, "lower .95"],
  hr_ci_high = cox_geo[1, "upper .95"],
  t(as.data.frame(med))
)
write.csv(stats, "results/GEO_Validation/GEO_summary_stats.csv", row.names = FALSE)
print(stats)

message("=== GEO Validation complete ===")
