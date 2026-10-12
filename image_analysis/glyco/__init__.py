"""glyco: quantitative analysis of terminal-monosaccharide multicolour fluorescence images.

The pipeline has three stages, each reading only what the one before wrote:

    extraction   raw images          ->  data/<cohort>_patches_*.csv, <cohort>_cases.csv
    modelling    patch tables        ->  models/*.joblib
    results      tables and models   ->  results/*.csv, results/RESULTS.md

Run from the `image_analysis/` directory:  python -m glyco --help
"""

__version__ = "1.0.0"
