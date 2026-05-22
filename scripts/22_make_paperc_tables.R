#!/usr/bin/env Rscript
# scripts/22_make_paperc_tables.R
# ------------------------------------------------------------
# Paper C tables:
# - Scenario table (A/B/C) expanded with mis_spec and delta_dgm
# - Copies key summaries into tables/
# ------------------------------------------------------------

source("R/00_utils.R")
cfg <- init_project(n_cores = 1L, seed = 2026L, out_dir = Sys.getenv("OUT_DIR", "output_paperc"))
cfg$require_pkgs(c("data.table"))

source("R/01_scenarios.R")

tab_dir <- file.path(cfg$out_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

grid <- make_scenario_grid()
mis <- as.integer(trimws(strsplit(Sys.getenv("MIS_SPEC_LEVELS","0,1"), ",", fixed = TRUE)[[1]]))
mis <- mis[mis %in% c(0L,1L)]
if (length(mis) == 0) mis <- 0:1

dd <- trimws(strsplit(Sys.getenv("DELTA_DGM_LEVELS","D1_zero,D2_constant,D3_late_only"), ",", fixed = TRUE)[[1]])
dd <- dd[dd %in% c("D1_zero","D2_constant","D3_late_only","D4_sign_change")]
if (length(dd) == 0) dd <- c("D1_zero","D2_constant","D3_late_only")

scn <- data.table::CJ(scenario_id = grid$scenario_id, mis_spec = mis, delta_dgm = dd)
scn <- merge(scn, grid, by = "scenario_id", all.x = TRUE)
data.table::setorder(scn, scenario_num, mis_spec, delta_dgm)

data.table::fwrite(
  scn[, .(scenario_id, mis_spec, delta_dgm, axisA, axisB, rho_meas, panel_title)],
  file.path(tab_dir, "TableC1_scenarios.csv")
)

# Copy summaries if present
for (nm in c("perf_summary_paperc_ipcw.csv",
             "perf_summary_paperc_dr.csv",
             "perf_summary_paperc_tmle.csv",
             "perf_summary_paperc_delta_star.csv")) {
  p <- file.path(cfg$out_dir, nm)
  if (file.exists(p)) data.table::fwrite(data.table::fread(p, showProgress = FALSE), file.path(tab_dir, nm))
}

# Copy MC truth if present
tp <- file.path(cfg$out_dir, "raw", "mc_truth_risk_paperc.csv")
if (file.exists(tp)) data.table::fwrite(data.table::fread(tp, showProgress = FALSE), file.path(tab_dir, "mc_truth_risk_paperc.csv"))

message("Saved tables under: ", tab_dir)
