#!/usr/bin/env bash
# Reproduce every statistic, table and figure in the manuscript from the raw task files.
#   DATA_DIR  folder with the raw data (default ./data)
#   OUT_DIR   output folder (default ./outputs)
# Usage:  DATA_DIR=/path/to/data OUT_DIR=/path/to/outputs ./run_all.sh
# Total run time is a few minutes (most of it the PRL model fits and the Bayes factors).
set -euo pipefail
cd "$(dirname "$0")/code"
export DATA_DIR="${DATA_DIR:-$(cd .. && pwd)/data}"
export OUT_DIR="${OUT_DIR:-$(cd .. && pwd)/outputs}"
export LC_ALL=C.UTF-8 2>/dev/null || true

echo "== 1. Scoring (Python): raw task files -> one row per participant x session"
python3 01_allocation.py
python3 02_score_prl.py
python3 03_score_rule_switching.py
python3 04_score_demand_selection.py
python3 05_score_msit.py
python3 06_score_ospan.py
python3 07_assemble_measures.py

echo "== 2. Statistics (R)"
Rscript 10_drug_effects_and_coupling.R
Rscript 11_sample_and_clinical.R
Rscript 12_prediction_and_wmc.R
Rscript 13_demand_avoidance.R
Rscript 14_bayes_factors.R
Rscript 15_sensitivity_rule_switching.R
Rscript 16_test_statistics.R
Rscript 17_supplementary_analyses.R

echo "== 3. Task validity and construct validation (Python)"
python3 20_task_validity.py
python3 21_construct_validation.py

echo "== 4. Figures (R)"
Rscript 30_figures.R
echo "Done. Results in $OUT_DIR/results, figures in $OUT_DIR/figures."
