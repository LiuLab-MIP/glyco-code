##############################################################################
# 00_install_packages.R
# Install every R package the pipeline needs
##############################################################################
# Run from the bioinformatics/ directory: all paths are relative to it.
# Packages are installed into the first library in .libPaths(); set R_LIBS_USER
# to install into a separate, project-specific library.

# Install BiocManager if it is not installed yet
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_pkgs <- c(
  "TCGAbiolinks", # 01: GDC query, download and preparation
  "SummarizedExperiment", # 01: counts and clinical data
  "S4Vectors", # needed to read the TCGA clinical DataFrame
  "GSVA", # 02, 06: ssGSEA (>= 2.0, uses ssgseaParam())
  "DESeq2", # 03: differential expression
  "clusterProfiler", # 04: GO, KEGG, GSEA
  "org.Hs.eg.db", # 04: gene identifier mapping
  "enrichplot", # 04: enrichment plots
  "DOSE", # 04
  "msigdbr", # 04: Hallmark gene sets
  "GEOquery" # 06: GSE71729
)

# CRAN packages
cran_pkgs <- c(
  "survival", # 05, 06
  "survminer", # 05, 06: Kaplan-Meier plots
  "ggplot2",
  "pheatmap", # 02, 06
  "dplyr",
  "ggrepel", # 03: volcano labels
  "ggridges", # 04: needed by enrichplot::ridgeplot()
  "svglite", # every script writes SVG figures
  "broom", # 05: ggforest_custom(): tidy() / glance()
  "ggpubr" # 05: ggforest_custom(): as_ggplot()
)

# Install the Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing Bioconductor package: ", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(paste0("Already installed: ", pkg))
  }
}

# Install the CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing CRAN package: ", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message(paste0("Already installed: ", pkg))
  }
}

# Fail loudly if any package is still missing (an install can fail without stopping the loop)
missing <- Filter(function(pkg) !requireNamespace(pkg, quietly = TRUE), c(bioc_pkgs, cran_pkgs))
if (length(missing) > 0) {
  stop("These packages could not be installed or loaded: ", paste(missing, collapse = ", "))
}

message("=== All packages installed successfully! ===")
