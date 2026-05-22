# scripts/98_quickcheck.R
# RStudio/Windows-friendly quick check for the Paper C pipeline.
# Run from RStudio Console with: source("scripts/98_quickcheck.R")

get_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(paste0("^", file_arg), cmd, value = TRUE)
  if (length(hit) > 0) {
    return(normalizePath(file.path(dirname(sub(file_arg, "", hit[1])), ".."), winslash = "/", mustWork = TRUE))
  }
  ofile <- NULL
  for (i in rev(seq_along(sys.frames()))) {
    if (!is.null(sys.frames()[[i]]$ofile)) {
      ofile <- sys.frames()[[i]]$ofile
      break
    }
  }
  if (!is.null(ofile)) {
    return(normalizePath(file.path(dirname(ofile), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- get_repo_root()
setwd(repo_root)

out_dir <- Sys.getenv("OUT_DIR", "quickcheck_paperc")
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OUT_DIR = out_dir,
  N_CORES = Sys.getenv("N_CORES", "2"),
  B = Sys.getenv("B", "10"),
  N = Sys.getenv("N", "500"),
  T_MAX = Sys.getenv("T_MAX", "5"),
  DT = Sys.getenv("DT", "0.25"),
  T_CUT = Sys.getenv("T_CUT", "5"),
  MIS_SPEC_LEVELS = Sys.getenv("MIS_SPEC_LEVELS", "0"),
  TRUTH_LEVELS = Sys.getenv("TRUTH_LEVELS", "null"),
  DELTA_DGM_LEVELS = Sys.getenv("DELTA_DGM_LEVELS", "D2_constant"),
  DELTA_SHAPES = Sys.getenv("DELTA_SHAPES", "constant"),
  DELTA_MAX = Sys.getenv("DELTA_MAX", "2"),
  DELTA_STEP = Sys.getenv("DELTA_STEP", "0.25"),
  BETA_TRUE_PP = Sys.getenv("BETA_TRUE_PP", "-0.2231436"),
  RUN_DR = Sys.getenv("RUN_DR", "1"),
  RUN_ML = Sys.getenv("RUN_ML", "0"),
  RUN_TMLE = Sys.getenv("RUN_TMLE", "0"),
  SAVE_SENS_CURVES = Sys.getenv("SAVE_SENS_CURVES", "1"),
  MC_TRUTH = Sys.getenv("MC_TRUTH", "0"),
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_log <- file.path(out_dir, "run.log")

run_step <- function(script) {
  cat("Running", script, "\n")
  cat("\n## Running ", script, "\n", file = run_log, append = TRUE, sep = "")
  out <- system2(rscript, script, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (length(out) > 0) cat(out, sep = "\n", file = run_log, append = TRUE)
  if (!identical(as.integer(status), 0L)) {
    cat("\nQuickcheck failed while running:", script, "\n")
    cat("See log:", normalizePath(run_log, winslash = "/", mustWork = FALSE), "\n")
    stop("quickcheck failed", call. = FALSE)
  }
}

run_step("scripts/21_run_paperc_simulations.R")
run_step("scripts/22_make_paperc_tables.R")
run_step("scripts/23_make_paperc_figures.R")
run_step("scripts/24_make_paperc_method_compare_figures.R")

# Minimal smoke-test validation. These checks do not validate the published
# numerical results; they only confirm that the quickcheck pipeline produced
# the expected output classes and that the residual-sensitivity engine returned
# finite curves.
validate_quickcheck <- function(out_dir) {
  required <- c(
    "raw/replicate_results_paperc.csv",
    "raw/spd_curves_paperc.csv",
    "raw/weight_diagnostics_paperc.csv",
    "raw/delta_star_paperc.csv",
    "raw/sensitivity_curves_paperc.csv",
    "tables/perf_summary_paperc_ipcw.csv",
    "tables/perf_summary_paperc_dr.csv",
    "tables/perf_summary_paperc_delta_star.csv",
    "figures/Figure1.pdf",
    "figures/Figure2.pdf",
    "figures/Figure3.pdf",
    "figures/Figure4.pdf"
  )
  missing <- required[!file.exists(file.path(out_dir, required))]
  if (length(missing) > 0) {
    stop(paste("Quickcheck did not create expected output(s):", paste(missing, collapse = ", ")), call. = FALSE)
  }

  sc <- data.table::fread(file.path(out_dir, "raw/sensitivity_curves_paperc.csv"), showProgress = FALSE)
  if (!("logrr" %in% names(sc)) || !any(is.finite(sc$logrr))) {
    stop("Quickcheck sensitivity curves contain no finite logRR values.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_quickcheck(out_dir)

cat("DONE ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = run_log, append = TRUE)
cat("Quickcheck completed. Outputs under:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
