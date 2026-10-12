##############################################################################
# 06_survival_analysis.R
# Survival analysis: Kaplan-Meier curves and Cox regression
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

message("Loading clinical and subtyping data...")
dat <- readRDS("data/TCGA_PAAD_subtyped.rds")
clin <- as.data.frame(dat$clin)

# Pre-process the clinical data: overall survival time (OS) and status
# In TCGA, OS.time = days_to_death or days_to_last_follow_up
clin$OS.time <- ifelse(is.na(clin$days_to_death), clin$days_to_last_follow_up, clin$days_to_death)
clin$OS <- ifelse(clin$vital_status == "Dead", 1, 0)

# Drop samples with missing survival data or a survival time of zero
clin_surv <- clin %>%
  filter(!is.na(OS.time) & OS.time > 0 & !is.na(OS) & !is.na(Final_Subtype))

message(sprintf("Running Survival Analysis on %d samples...", nrow(clin_surv)))

# 1. Kaplan-Meier survival curves
dir.create("results/Survival", recursive = TRUE, showWarnings = FALSE)

fit <- survfit(Surv(OS.time / 30, OS) ~ Final_Subtype, data = clin_surv) # time converted to months (OS.time / 30)

# The log-rank statistics of this comparison (chi-square, degrees of freedom, hazard ratio) are in Moffitt/survival.csv, row "glycophenotype".

# Plot with ggsurvplot
p <- ggsurvplot(fit,
  data = clin_surv,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = FALSE,
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
  ggtheme = theme_classic()
)

# Custom legend labels with plotmath (upright superscripts for high / low)
labels_expr <- c(
  expression(plain(Fuc)^plain(high) ~ plain(Man)^plain(low)),
  expression(plain(Fuc)^plain(low) ~ plain(Man)^plain(high))
)
p$plot <- p$plot +
  scale_colour_manual(values = c("#d7301f", "#0570b0"), labels = labels_expr) +
  scale_fill_manual(values = c("#d7301f", "#0570b0"), labels = labels_expr)

pdf("results/Survival/KM_Survival_Plot.pdf", width = 6, height = 4, onefile = FALSE)
print(p)
dev.off()

svglite::svglite("results/Survival/KM_Survival_Plot.svg", width = 6, height = 4)
print(p)
dev.off()

# 2. Cox proportional-hazards regression
# Add the other clinical factors (stage and age) to the model
# Pre-process the clinical variables
clin_surv$Age <- as.numeric(clin_surv$age_at_index)

# Stage: remove NA and "not reported"
clin_surv$Stage <- as.character(clin_surv$ajcc_pathologic_stage)
clin_surv$Stage[grepl("stage i$|stage ia|stage ib", clin_surv$Stage, ignore.case = TRUE)] <- "Stage I"
clin_surv$Stage[grepl("stage ii$|stage iia|stage iib", clin_surv$Stage, ignore.case = TRUE)] <- "Stage II"
clin_surv$Stage[grepl("stage iii", clin_surv$Stage, ignore.case = TRUE)] <- "Stage III"
clin_surv$Stage[grepl("stage iv", clin_surv$Stage, ignore.case = TRUE)] <- "Stage IV"
clin_surv$Stage[clin_surv$Stage %in% c("not reported", "discrepancy", "Stage 0", "Stage 0a")] <- NA

clin_surv_cox <- clin_surv %>% filter(Stage %in% c("Stage I", "Stage II") & !is.na(Age))
clin_surv_cox$Glycophenotype <- factor(
  clin_surv_cox$Final_Subtype,
  levels = c("Fuc-favored", "Man-favored"),
  labels = c("Fuchigh Manlow", "Fuclow Manhigh")
)
ggforest_custom <- function(
  model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7, refLabel = "Reference", noDigits = 2
) {
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    } else if (terms[i] == "numeric") {
      data.frame(
        var = var, Var1 = "", Freq = nrow(data),
        pos = 1
      )
    } else {
      vars <- grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(
        var = vars, Var1 = "", Freq = nrow(data),
        pos = seq_along(vars)
      )
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  # Custom: exact match, avoid partial row-name matching (Stage I -> Stage II).
  # Guard the mapping: every coefficient must be used exactly once, and each
  # factor/character term must contribute exactly one reference (NA) row.
  # Catches both the old duplication bug and a silent name-mismatch that would
  # render a real estimate as "Reference".
  matched <- match(inds, rownames(coef))
  used <- matched[!is.na(matched)]
  if (!setequal(used, seq_len(nrow(coef))) || anyDuplicated(used) ||
    sum(is.na(matched)) != sum(terms %in% c("factor", "character"))) {
    stop(
      "ggforest_custom: model coefficients do not map 1:1 onto plotted rows.\n",
      "  coefficients: ", paste(rownames(coef), collapse = ", "), "\n",
      "  plotted rows: ", paste(inds, collapse = ", ")
    )
  }
  toShow <- cbind(allTermsDF, coef[matched, ])[, c(
    "var", "level",
    "N", "p.value", "estimate", "conf.low", "conf.high",
    "pos"
  )]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp[] <- lapply(toShowExp, function(x) sprintf("%.2f", exp(x)))
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[
    ,
    4
  ], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(
    toShowExpClean$p.value,
    noDigits + 1
  ), " ", ifelse(toShowExpClean$p.value < 0.05,
    "*", ""
  ), ifelse(toShowExpClean$p.value < 0.01, "*",
    ""
  ), ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0(
    "(", toShowExpClean[, "conf.low.1"],
    " - ", toShowExpClean[, "conf.high.1"], ")"
  )
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] <- refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] <- "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] <- 0
  toShowExpClean$var <- as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] <- ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high,
    na.rm = TRUE
  )
  breaks <- axisTicks(rangeb / 2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - 1.5 * diff(rangeb) # Custom: expanded left text margin
  rangeplot[2] <- rangeplot[2] + 0.22 * diff(rangeb) # Custom: expanded right margin for p-values
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(grid::unit(
    theme_get()$text$size,
    "pt"
  ), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    geom_rect(aes(
      xmin = seq_along(var) - 0.5, xmax = seq_along(var) +
        0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
      fill = ordered(seq_along(var) %% 2 + 1)
    )) +
    scale_fill_manual(values = c(
      "#FFFFFF33",
      "#00000033"
    ), guide = "none") +
    geom_point(
      pch = 15,
      size = 4
    ) +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)),
      width = 0.15
    ) +
    geom_hline(yintercept = 1, linetype = 3) +
    coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    scale_y_log10(
      name = "",
      labels = sprintf("%g", breaks), expand = c(0.02, 0.02),
      breaks = breaks
    ) +
    scale_x_continuous(limits = c(0.4, nrow(toShowExpClean) + 0.65), expand = c(0, 0)) +
    theme_light() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
      legend.position = "none", panel.border = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6, "pt")
    ) +
    xlab("") +
    annotate(
      geom = "text", x = x_annotate, y = exp(y_variable),
      label = toShowExpClean$var, fontface = "bold", hjust = 0,
      size = annot_size_mm
    ) +
    annotate(
      geom = "text", x = x_annotate,
      y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level,
      vjust = -0.3, size = annot_size_mm
    ) +
    annotate(
      geom = "text",
      x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N,
      fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level ==
        "", 0.5, 1.5), size = annot_size_mm
    ) +
    annotate(
      geom = "text",
      x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1,
      size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 ==
        "Reference", 0.5, -0.1)
    ) +
    annotate(
      geom = "text",
      x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci,
      size = annot_size_mm, vjust = 1.5, fontface = "italic"
    ) +
    annotate(
      geom = "text", x = x_annotate, y = exp(y_stars),
      label = toShowExpClean$stars, size = annot_size_mm,
      hjust = -0.2, fontface = "italic"
    )
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

clin_surv_cox$Stage <- as.factor(clin_surv_cox$Stage)

message("Running Multivariate Cox Regression...")
tryCatch(
  {
    # Time in months, matching the KM plot above. Cox partial likelihood is
    # rank-based, so rescaling time leaves all HRs unchanged; this is purely
    # so time units are reported consistently across figures.
    cox_model <- coxph(Surv(OS.time / 30, OS) ~ Glycophenotype + Age + Stage, data = clin_surv_cox)
    summary_cox <- summary(cox_model)

    # Cox forest plot
    pdf("results/Survival/Cox_Forest_Plot.pdf", width = 7.5, height = 4, onefile = FALSE)
    p <- ggforest_custom(cox_model, data = clin_surv_cox, main = "", fontsize = 1.0, cpositions = c(0.02, 0.23, 0.47))
    print(p)
    dev.off()

    svglite::svglite("results/Survival/Cox_Forest_Plot.svg", width = 7.5, height = 4)
    p <- ggforest_custom(cox_model, data = clin_surv_cox, main = "", fontsize = 1.0, cpositions = c(0.02, 0.23, 0.47))
    print(p)
    dev.off()

    # coefficients plus the 95% confidence interval of each hazard ratio
    write.csv(
      cbind(as.data.frame(summary_cox$coefficients), as.data.frame(summary_cox$conf.int[, c("lower .95", "upper .95")])),
      "results/Survival/Cox_regression_summary.csv"
    )
  },
  error = function(e) {
    message("Cox regression failed, likely due to perfect separation (zero events in a subgroup).")
    message("Caught error: ", e$message)
  }
)

message("=== Survival Analysis complete ===")
