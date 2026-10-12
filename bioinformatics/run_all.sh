#!/usr/bin/env bash
# Run the bioinformatics pipeline from the shipped TCGA-PAAD input.
#
#   bash run_all.sh                run steps 02-07 (about 2.5 minutes)
#   bash run_all.sh --install      install the R packages first (step 00)
#   bash run_all.sh --download     rebuild data/TCGA_PAAD_processed.rds from GDC first (step 01, ~800 MB)
#
# Steps 04 (KEGG) and 07 (GEO) need internet access. Results are written to results/.
set -euo pipefail
cd "$(dirname "$0")"

if [[ " $* " == *" --install "* ]]; then
  Rscript 00_install_packages.R
fi
if [[ " $* " == *" --download "* ]]; then
  Rscript 01_data_download.R
fi

mkdir -p results
for script in 02_ssGSEA_scoring.R 03_DEG_analysis.R 04_enrichment_analysis.R \
              05_moffitt_scoring.R 06_survival_analysis.R 07_GEO_validation.R; do
  echo "=== ${script}"
  Rscript "${script}"
done

# Record the exact R and package versions of this run
Rscript -e 'sink("results/sessionInfo.txt"); print(sessionInfo()); sink()'
# Compare the tables with the reference tables
Rscript compare_with_expected.R
echo "Done. Tables and figures are in results/; reference tables are in expected_results/."
