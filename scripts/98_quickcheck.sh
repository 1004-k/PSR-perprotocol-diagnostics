#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
OUT_DIR=${OUT_DIR:-quickcheck_paperc}
rm -rf "$OUT_DIR" || true
mkdir -p "$OUT_DIR"

export OUT_DIR
export N_CORES=${N_CORES:-2}
export B=${B:-10}
export N=${N:-500}
export T_MAX=${T_MAX:-5}
export DT=${DT:-0.25}
export T_CUT=${T_CUT:-5}
export MIS_SPEC_LEVELS=${MIS_SPEC_LEVELS:-"0"}
export TRUTH_LEVELS=${TRUTH_LEVELS:-"null"}
export DELTA_DGM_LEVELS=${DELTA_DGM_LEVELS:-"D2_constant"}
export DELTA_SHAPES=${DELTA_SHAPES:-"constant"}
export DELTA_MAX=${DELTA_MAX:-2}
export DELTA_STEP=${DELTA_STEP:-0.25}
export BETA_TRUE_PP=${BETA_TRUE_PP:-"-0.2231436"}
export RUN_DR=${RUN_DR:-1}
export RUN_ML=${RUN_ML:-0}
export RUN_TMLE=${RUN_TMLE:-0}
export SAVE_SENS_CURVES=${SAVE_SENS_CURVES:-1}
export MC_TRUTH=${MC_TRUTH:-0}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript scripts/21_run_paperc_simulations.R > "$OUT_DIR/run.log" 2>&1
Rscript scripts/22_make_paperc_tables.R     >> "$OUT_DIR/run.log" 2>&1
Rscript scripts/23_make_paperc_figures.R
Rscript scripts/24_make_paperc_method_compare_figures.R    >> "$OUT_DIR/run.log" 2>&1


for f in \
  raw/replicate_results_paperc.csv \
  raw/spd_curves_paperc.csv \
  raw/weight_diagnostics_paperc.csv \
  raw/delta_star_paperc.csv \
  raw/sensitivity_curves_paperc.csv \
  tables/perf_summary_paperc_ipcw.csv \
  tables/perf_summary_paperc_dr.csv \
  tables/perf_summary_paperc_delta_star.csv \
  figures/Figure1.pdf \
  figures/Figure2.pdf \
  figures/Figure3.pdf
do
  if [ ! -f "$OUT_DIR/$f" ]; then
    echo "Missing expected quickcheck output: $OUT_DIR/$f" | tee -a "$OUT_DIR/run.log"
    exit 1
  fi
done

echo "DONE $(date)" | tee -a "$OUT_DIR/run.log"
