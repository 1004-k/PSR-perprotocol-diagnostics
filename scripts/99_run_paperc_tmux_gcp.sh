#!/usr/bin/env bash
# scripts/99_run_paperc_tmux_gcp.sh
# ------------------------------------------------------------
# GCP/tmux launcher for Paper C (main + methods + sensitivity).
# Usage (on the VM):
#   cd ~/spdt_final_gcp_github
#   bash scripts/99_run_paperc_tmux_gcp.sh
#
# Notes:
# - Avoid expressions like log(0.8) in env vars; use numeric strings.
# - Prevent thread oversubscription (OMP/BLAS/MKL threads set to 1).
# ------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

tmux kill-session -t pc_main 2>/dev/null || true
tmux kill-session -t pc_methods 2>/dev/null || true
tmux kill-session -t pc_sens 2>/dev/null || true

# ---- common knobs ----
export T_MAX="${T_MAX:-5}"
export DT="${DT:-0.25}"
export T_CUT="${T_CUT:-5}"
export MIS_SPEC_LEVELS="${MIS_SPEC_LEVELS:-0,1}"
export TRUTH_LEVELS="${TRUTH_LEVELS:-null,non_null}"
export BETA_TRUE_PP="${BETA_TRUE_PP:--0.2231436}"

export DELTA_DGM_LEVELS="${DELTA_DGM_LEVELS:-D1_zero,D2_constant,D3_late_only}"
export DELTA_MAG="${DELTA_MAG:-0.6}"
export T0_DELTA="${T0_DELTA:-3.5}"
export CORR_MEAS_UNM="${CORR_MEAS_UNM:-0.0}"

export BREAK_DT="${BREAK_DT:-1.0}"
export DELTA_SHAPES="${DELTA_SHAPES:-constant,late_only}"
export DELTA_MAX="${DELTA_MAX:-2.0}"
export DELTA_STEP="${DELTA_STEP:-0.1}"
export T0_SENS="${T0_SENS:-3.5}"
export DECISION_MODE="${DECISION_MODE:-cross_zero}"

export MC_TRUTH="${MC_TRUTH:-1}"
export N_TRUTH="${N_TRUTH:-200000}"

# anti oversubscription
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# ---- pc_main: full grid, IPCW + DR(glm), delta* summary ----
tmux new-session -d -s pc_main "
cd '${REPO_ROOT}' && set -euo pipefail
OUT_DIR='output_paperc_main'
mkdir -p \"\${OUT_DIR}\"
export OUT_DIR
export N_CORES='${N_CORES_MAIN:-10}'
export B='${B_MAIN:-200}'
export N='${N_MAIN:-2000}'
export RUN_DR=1
export RUN_ML=0
export RUN_TMLE=0
export SAVE_SENS_CURVES='${SAVE_SENS_CURVES_MAIN:-0}'
export TIME_RESOLVED='${TIME_RESOLVED_MAIN:-0}'

Rscript scripts/21_run_paperc_simulations.R > \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/22_make_paperc_tables.R     >> \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/23_make_paperc_figures.R    >> \"\${OUT_DIR}/run.log\" 2>&1
echo \"DONE \$(date)\" >> \"\${OUT_DIR}/run.log\"
"

# ---- pc_methods: subset scenarios, adds ML/TMLE, saves sensitivity curves for fragility figure ----
tmux new-session -d -s pc_methods "
cd '${REPO_ROOT}' && set -euo pipefail
OUT_DIR='output_paperc_methods'
mkdir -p \"\${OUT_DIR}\"
export OUT_DIR
export N_CORES='${N_CORES_METHODS:-6}'
export B='${B_METHODS:-50}'
export N='${N_METHODS:-2000}'
export SCENARIO_SUBSET='${SCENARIO_SUBSET:-S01,S09,S18}'
export RUN_DR=1
export RUN_ML=1
export ML_METHOD='${ML_METHOD:-glmnet}'
export CROSSFIT='${CROSSFIT:-1}'
export CF_FOLDS='${CF_FOLDS:-2}'
export RUN_TMLE=1
export SAVE_SENS_CURVES=1
export TIME_RESOLVED='${TIME_RESOLVED_METHODS:-0}'

Rscript scripts/21_run_paperc_simulations.R > \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/22_make_paperc_tables.R     >> \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/23_make_paperc_figures.R    >> \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/24_make_paperc_method_compare_figures.R >> \"\${OUT_DIR}/run.log\" 2>&1
echo \"DONE \$(date)\" >> \"\${OUT_DIR}/run.log\"
"

# ---- pc_sens: optional heavier sensitivity check (sign-change D4 + time-resolved delta*(t)) ----
tmux new-session -d -s pc_sens "
cd '${REPO_ROOT}' && set -euo pipefail
OUT_DIR='output_paperc_sensitivity'
mkdir -p \"\${OUT_DIR}\"
export OUT_DIR
export N_CORES='${N_CORES_SENS:-6}'
export B='${B_SENS:-50}'
export N='${N_SENS:-2000}'
export SCENARIO_SUBSET='${SCENARIO_SUBSET_SENS:-S01,S09,S18}'
export DELTA_DGM_LEVELS='${DELTA_DGM_LEVELS_SENS:-D4_sign_change}'
export RUN_DR=1
export RUN_ML=0
export RUN_TMLE=0
export SAVE_SENS_CURVES=1
export TIME_RESOLVED=1
export DELTA_SHAPES='${DELTA_SHAPES_SENS:-late_only,sign_change}'

Rscript scripts/21_run_paperc_simulations.R > \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/22_make_paperc_tables.R     >> \"\${OUT_DIR}/run.log\" 2>&1
Rscript scripts/23_make_paperc_figures.R    >> \"\${OUT_DIR}/run.log\" 2>&1
echo \"DONE \$(date)\" >> \"\${OUT_DIR}/run.log\"
"

tmux ls
echo ""
echo "To monitor:"
echo "  tail -f output_paperc_main/run.log"
echo "  tail -f output_paperc_methods/run.log"
echo "  tail -f output_paperc_sensitivity/run.log"
