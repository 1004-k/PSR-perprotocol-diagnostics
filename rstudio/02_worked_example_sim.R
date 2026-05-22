# RStudio wrapper: run simulation-based worked example (single scenario/replicate)
# ------------------------------------------------------------
# Produces: OUT_DIR/raw and OUT_DIR/figures for one illustrative run.
# ------------------------------------------------------------

OUT_DIR <- "output_paperc_worked_sim"  # <-- change if needed
SCENARIO_ID <- "S09"
MIS_SPEC <- "0"
TRUTH <- "null"
DELTA_DGM <- "D2_constant"

# Core simulation size for demo
N <- "2000"
T_MAX <- "5"
DT <- "0.25"
T_CUT <- "5"

Sys.setenv(
  OUT_DIR = OUT_DIR,
  SCENARIO_ID = SCENARIO_ID,
  MIS_SPEC = MIS_SPEC,
  TRUTH = TRUTH,
  DELTA_DGM = DELTA_DGM,
  N = N,
  T_MAX = T_MAX,
  DT = DT,
  T_CUT = T_CUT
)

source("scripts/25_run_paperc_worked_example.R")

message("Done. Worked example saved under: ", OUT_DIR)
