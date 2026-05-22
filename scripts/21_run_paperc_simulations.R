#!/usr/bin/env Rscript
# scripts/21_run_paperc_simulations.R
# ------------------------------------------------------------
# Paper C runner: PSR bundle + residual non-ignorability tipping summaries
#
# Outputs (OUT_DIR):
#   raw/replicate_results_paperc.csv
#   raw/spd_curves_paperc.csv
#   raw/weight_diagnostics_paperc.csv
#   raw/sensitivity_curves_paperc.csv          (optional; SAVE_SENS_CURVES=1)
#   raw/mc_truth_risk_paperc.csv               (cached)
#   perf_summary_paperc_ipcw.csv
#   perf_summary_paperc_dr.csv                 (optional)
#   perf_summary_paperc_delta_star.csv
# ------------------------------------------------------------

# ---- load modules robustly ----
source("R/00_utils.R")
cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "4")),
  seed    = 2026L,
  out_dir = Sys.getenv("OUT_DIR", "output_paperc")
)
cfg$require_pkgs(c("data.table", "survival", "future.apply", "future"))

# Core modules reused from Paper B
source("R/01_scenarios.R")
source("R/02_dgp_simulate.R")      # for make_gamma_path()
source("R/03_spd_curve.R")
source("R/04_pp_ipcw.R")
source("R/05_tipping.R")
source("R/06_performance.R")
source("R/09_dr_aipw_ml.R")
source("R/07_plotting.R")          # for open_pdf/close_device (used by downstream scripts)

# Paper C modules
source("R/10_dgp_paperc.R")
source("R/11_paperc_sensitivity.R")

# ---- knobs ----
B     <- as.integer(Sys.getenv("B", "200"))
N     <- as.integer(Sys.getenv("N", "2000"))
t_max <- as.numeric(Sys.getenv("T_MAX", "5"))
dt    <- as.numeric(Sys.getenv("DT", "0.25"))
T_CUT <- as.numeric(Sys.getenv("T_CUT", as.character(t_max)))

# Truth levels (null/non_null)
parse_truth_levels <- function() {
  raw <- trimws(Sys.getenv("TRUTH_LEVELS", unset = "null"))
  if (!nzchar(raw) || tolower(raw) %in% c("null", "none")) return(c("null"))
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out <- out[out %in% c("null","non_null")]
  if (length(out) == 0) out <- c("null")
  out
}
truth_levels <- parse_truth_levels()

beta_true_pp <- as.numeric(Sys.getenv("BETA_TRUE_PP", "0"))

# Mis-spec levels (affects nuisance models for IPCW/DR)
MIS_SPEC_LEVELS <- strsplit(Sys.getenv("MIS_SPEC_LEVELS", "0,1"), ",", fixed = TRUE)[[1]]
MIS_SPEC_LEVELS <- as.integer(trimws(MIS_SPEC_LEVELS))
MIS_SPEC_LEVELS <- MIS_SPEC_LEVELS[MIS_SPEC_LEVELS %in% c(0L, 1L)]
if (length(MIS_SPEC_LEVELS) == 0) MIS_SPEC_LEVELS <- c(0L, 1L)

# Paper C: residual MNAR DGP axis
DELTA_DGM_LEVELS <- trimws(strsplit(Sys.getenv("DELTA_DGM_LEVELS", "D1_zero,D2_constant,D3_late_only"), ",", fixed = TRUE)[[1]])
DELTA_DGM_LEVELS <- DELTA_DGM_LEVELS[DELTA_DGM_LEVELS %in% c("D1_zero","D2_constant","D3_late_only","D4_sign_change")]
if (length(DELTA_DGM_LEVELS) == 0) DELTA_DGM_LEVELS <- c("D1_zero","D2_constant","D3_late_only")

DELTA_MAG <- as.numeric(Sys.getenv("DELTA_MAG", "0.6"))       # magnitude for the DGP axis D (unmeasured effect)
T0_DELTA  <- as.numeric(Sys.getenv("T0_DELTA", "3.5"))        # default late-only switch time (approx 0.7*5)

# DGP correlation between measured and unmeasured components
CORR_MEAS_UNM <- as.numeric(Sys.getenv("CORR_MEAS_UNM", "0.0"))

# SPD(t) breaks (align delta(t) piecewise to these breaks)
BREAK_DT <- as.numeric(Sys.getenv("BREAK_DT", "1.0"))
breaks <- seq(0, t_max, by = BREAK_DT)
if (tail(breaks, 1) < t_max) breaks <- c(breaks, t_max)

# tipping thresholds for weight instability
c_ess  <- as.numeric(Sys.getenv("C_ESS", "0.25"))
c_tail <- as.numeric(Sys.getenv("C_TAIL", "0.10"))

# Comparators
RUN_DR    <- as.integer(Sys.getenv("RUN_DR", "1")) == 1L
RUN_ML    <- as.integer(Sys.getenv("RUN_ML", "0")) == 1L
ML_METHOD <- Sys.getenv("ML_METHOD", "glmnet")
CROSSFIT  <- as.integer(Sys.getenv("CROSSFIT", "0")) == 1L
CF_FOLDS  <- as.integer(Sys.getenv("CF_FOLDS", "2"))
RUN_TMLE  <- as.integer(Sys.getenv("RUN_TMLE", "0")) == 1L

# Sensitivity grid (delta for the residual MNAR model used in Paper C)
DELTA_MAX  <- as.numeric(Sys.getenv("DELTA_MAX", "2.0"))
DELTA_STEP <- as.numeric(Sys.getenv("DELTA_STEP", "0.1"))
delta_grid <- seq(0, DELTA_MAX, by = DELTA_STEP)

DELTA_SHAPES <- trimws(strsplit(Sys.getenv("DELTA_SHAPES", "constant,late_only"), ",", fixed = TRUE)[[1]])
DELTA_SHAPES <- DELTA_SHAPES[DELTA_SHAPES %in% c("constant","late_only","sign_change")]
if (length(DELTA_SHAPES) == 0) DELTA_SHAPES <- c("constant","late_only")
T0_SENS <- as.numeric(Sys.getenv("T0_SENS", as.character(0.7 * t_max)))

DECISION_MODE <- Sys.getenv("DECISION_MODE", "cross_zero")  # cross_zero or sign
SAVE_SENS_CURVES <- as.integer(Sys.getenv("SAVE_SENS_CURVES", "0")) == 1L
TIME_RESOLVED <- as.integer(Sys.getenv("TIME_RESOLVED", "0")) == 1L  # compute delta*(t) curve (expensive)

# Monte Carlo truth (per-protocol, no deviation)
MC_TRUTH <- as.integer(Sys.getenv("MC_TRUTH", "1")) == 1L
N_TRUTH  <- as.integer(Sys.getenv("N_TRUTH", "200000"))
SEED_TRUTH <- as.integer(Sys.getenv("SEED_TRUTH", "20261"))

# scenario grid (Axis A/B/C from Paper B)
grid <- make_scenario_grid()
scn_grid <- data.table::CJ(
  scenario_id = grid$scenario_id,
  mis_spec = MIS_SPEC_LEVELS,
  delta_dgm = DELTA_DGM_LEVELS,
  truth = truth_levels,
  unique = TRUE
)
scn_grid <- merge(scn_grid, grid, by = "scenario_id", all.x = TRUE)
data.table::setorder(scn_grid, scenario_num, mis_spec, delta_dgm, truth)

# optional scenario subset
subset_ids <- trimws(strsplit(Sys.getenv("SCENARIO_SUBSET", ""), ",", fixed = TRUE)[[1]])
subset_ids <- subset_ids[nzchar(subset_ids)]
if (length(subset_ids) > 0) {
  scn_grid <- scn_grid[scenario_id %in% subset_ids]
  data.table::setorder(scn_grid, scenario_num, mis_spec, delta_dgm, truth)
}

# ---- output dirs ----
raw_dir <- file.path(cfg$out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- MC truth cache ----
mc_truth <- NULL
truth_path <- file.path(raw_dir, "mc_truth_risk_paperc.csv")
if (MC_TRUTH && file.exists(truth_path)) {
  tmp <- tryCatch(data.table::fread(truth_path, showProgress = FALSE), error = function(e) NULL)
  if (!is.null(tmp) && all(c("delta_dgm","delta_mag","t_cut","beta_true","logrr","rd") %in% names(tmp))) {
    ok <- isTRUE(all(tmp$t_cut == T_CUT)) && isTRUE(all(tmp$beta_true == beta_true_pp)) &&
      isTRUE(all(tmp$delta_mag == DELTA_MAG))
    if (ok) {
      mc_truth <- tmp
      message("Loaded cached MC truth: ", truth_path)
    }
  }
}

if (MC_TRUTH && is.null(mc_truth)) {
  message("Computing MC truth under adherence (no deviation) ...")
  truth_list <- lapply(DELTA_DGM_LEVELS, function(dd) {
    mc_truth_risk_paperc(
      N_truth = N_TRUTH, t_max = t_max, dt = dt, t_cut = T_CUT,
      beta_true = beta_true_pp,
      theta_meas = 0.50,
      delta_dgm = dd, delta_mag = DELTA_MAG, t0_delta = T0_DELTA,
      corr_meas_unm = CORR_MEAS_UNM,
      seed = SEED_TRUTH + sum(utf8ToInt(dd))
    )
  })
  mc_truth <- data.table::rbindlist(truth_list, fill = TRUE)
  data.table::fwrite(mc_truth, truth_path)
  message("Saved MC truth: ", truth_path)
}

# ---- parallel plan ----
future::plan(future::multisession, workers = cfg$n_cores)

TAG <- sprintf("B%d_N%d", B, N)
log_file <- file.path(cfg$log_dir, sprintf("paperc_progress_%s.log", TAG))
if (file.exists(log_file)) file.remove(log_file)
cfg$log_line(log_file, sprintf("[%s] Paper C started | OUT_DIR=%s | rows=%d | B=%d | N=%d\n",
                               format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cfg$out_dir, nrow(scn_grid), B, N))

cat("Paper C grid rows:", nrow(scn_grid), " | B=", B, " | N=", N, "\n")
cat("OUT_DIR:", cfg$out_dir, "\n")

# ---- one replicate ----
run_one <- function(scn_row, rep_id) {
  truth <- as.character(scn_row$truth)
  beta <- if (truth == "non_null") beta_true_pp else 0

  seed <- seed_for_job(cfg$seed * 200000L, scn_row$scenario_id,
                       truth = paste0(truth, "_", scn_row$delta_dgm),
                       rep_id = rep_id)

  sim <- simulate_one_dataset_paperc(
    N = N, t_max = t_max, dt = dt,
    beta_true = beta,
    axisA = scn_row$axisA,
    axisB = scn_row$axisB,
    rho_meas = scn_row$rho_meas,
    delta_dgm = scn_row$delta_dgm,
    delta_mag = DELTA_MAG,
    t0_delta = T0_DELTA,
    corr_meas_unm = CORR_MEAS_UNM,
    seed = seed
  )

  long_dt <- sim$long_dt
  pp_dt   <- sim$pp_dt

  # SPD(t) based on measured prognostic score
  spd_fit <- estimate_spd_piecewise(
    long_dt, breaks = breaks,
    z_col = "z_obs", dev_col = "dev", id_col = "id", robust = TRUE
  )
  spd_dt <- compute_cum_pressure(spd_fit$spd)

  # IPCW censoring model (mis-spec knob)
  cens_fml <- if (scn_row$mis_spec == 1L) {
    survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi
  } else {
    survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall
  }
  # IPCW risk at horizon (aligned with DR/TMLE estimands)
  fit_cens <- survival::coxph(cens_fml, data = pp_dt, model = TRUE, x = TRUE)
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
  lp_cens <- stats::predict(fit_cens, type = "lp")

  ipcw_risk <- fit_ipcw_risk(
    pp_dt = pp_dt, t_cut = T_CUT,
    fit_cens = fit_cens, lp = lp_cens,
    w_floor = 1e-3
  )

  # weight diagnostics on time grid
  t_eval <- sort(unique(long_dt$tstart))
  weights_dt <- data.table::CJ(id = pp_dt$id, t = t_eval)
  lp_dt <- data.table::data.table(id = pp_dt$id, lp = lp_cens)
  weights_dt <- merge(weights_dt, lp_dt, by = "id", all.x = TRUE)
  weights_dt[, S := get_S_from_cox(fit_cens, time_vec = t, lp_vec = lp)]
  weights_dt[, w := 1 / pmax(S, 1e-3)]

  diag_dt <- compute_weight_diagnostics(long_dt, weights_dt,
                                        time_col = "tstart", id_col = "id", w_col = "w", q_tail = 0.99)
  tip <- detect_tipping(diag_dt, c_ess = c_ess, c_tail = c_tail)

  # DR risk at horizon
  dr <- NULL
  if (RUN_DR) {
    dr <- fit_dr_risk_rescue(
      pp_dt = pp_dt,
      t_cut = T_CUT,
      mis_spec = scn_row$mis_spec,
      use_ml_Q = RUN_ML,
      ml_method = ML_METHOD,
      crossfit = CROSSFIT,
      cf_folds = CF_FOLDS,
      seed = seed + 19L
    )
  }

  # TMLE (optional)
  tm <- NULL
  if (RUN_TMLE) {
    tm <- fit_tmle_risk_rescue(pp_dt = pp_dt, t_cut = T_CUT, mis_spec = scn_row$mis_spec)
  }

  # Sensitivity (anchored)
  sens_rows <- list()
  delta_star_rows <- list()

  # Anchor to IPCW and DR separately
  base_methods <- list(
    IPCW = list(logrr = ipcw_risk$logrr, rd = ipcw_risk$rd)
  )
  if (!is.null(dr)) base_methods$DR <- list(logrr = dr$logrr, rd = dr$rd)
  if (!is.null(tm)) base_methods$TMLE <- list(logrr = tm$logrr, rd = tm$rd)

  for (mname in names(base_methods)) {
    th0 <- base_methods[[mname]]

    for (shape in DELTA_SHAPES) {
      curve <- anchored_theta_curve(
        pp_dt = pp_dt, t_cut = T_CUT, breaks = breaks,
        delta_grid = delta_grid, shape = shape, t0 = T0_SENS,
        theta0_logrr = th0$logrr, theta0_rd = th0$rd
      )

      ds_logrr <- compute_delta_star(delta_grid, curve, measure = "logrr", decision_mode = DECISION_MODE)
      ds_rd    <- compute_delta_star(delta_grid, curve, measure = "rd",    decision_mode = DECISION_MODE)

      delta_star_rows[[length(delta_star_rows)+1L]] <- data.table::data.table(
        scenario_id = scn_row$scenario_id,
        mis_spec = scn_row$mis_spec,
        delta_dgm = scn_row$delta_dgm,
        truth = truth,
        replicate = rep_id,
        method = mname,
        delta_shape = shape,
        delta_star_logrr = ds_logrr,
        delta_star_rd = ds_rd
      )

      if (SAVE_SENS_CURVES) {
        curve2 <- data.table::copy(curve)
        curve2[, `:=`(
          scenario_id = scn_row$scenario_id,
          mis_spec = scn_row$mis_spec,
          delta_dgm = scn_row$delta_dgm,
          truth = truth,
          replicate = rep_id,
          method = mname,
          delta_shape = shape
        )]
        sens_rows[[length(sens_rows)+1L]] <- curve2
      }

      if (TIME_RESOLVED && shape == "late_only") {
        # compute delta*(t) curve only for late-only (expensive)
        curve_t <- delta_star_curve_by_start(
          pp_dt = pp_dt, t_cut = T_CUT, breaks = breaks, delta_grid = delta_grid,
          theta0_logrr = th0$logrr, theta0_rd = th0$rd, decision_mode = DECISION_MODE
        )
        curve_t[, `:=`(
          scenario_id = scn_row$scenario_id,
          mis_spec = scn_row$mis_spec,
          delta_dgm = scn_row$delta_dgm,
          truth = truth,
          replicate = rep_id,
          method = mname,
          delta_shape = "late_by_start"
        )]
        # store as sensitivity rows (long)
        if (SAVE_SENS_CURVES) sens_rows[[length(sens_rows)+1L]] <- curve_t
      }
    }
  }

  list(
    rep = data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = scn_row$axisA, axisB = scn_row$axisB, rho_meas = scn_row$rho_meas,
      mis_spec = scn_row$mis_spec,
      delta_dgm = scn_row$delta_dgm,
      truth = truth,
      replicate = rep_id,
      # summaries for operating map
      max_spd = max(abs(spd_dt$gamma_hat), na.rm = TRUE),
      Gamma_end = max(spd_dt$Gamma_hat, na.rm = TRUE),
      tip_time = tip$t_star,
      min_rESS = min(tip$diag$rESS, na.rm = TRUE),
      max_tail = max(tip$diag$tail_share, na.rm = TRUE),
      # baseline effect estimates for reporting
      ipcw_logrr = ipcw_risk$logrr,
      ipcw_se_logrr = ipcw_risk$se_logrr,
      ipcw_rd = ipcw_risk$rd,
      ipcw_se_rd = ipcw_risk$se_rd,
      dr_logrr = if (!is.null(dr)) dr$logrr else NA_real_,
      dr_se_logrr = if (!is.null(dr)) dr$se_logrr else NA_real_,
      dr_rd = if (!is.null(dr)) dr$rd else NA_real_,
      dr_se_rd = if (!is.null(dr)) dr$se_rd else NA_real_,
      tmle_logrr = if (!is.null(tm)) tm$logrr else NA_real_,
      tmle_se_logrr = if (!is.null(tm)) tm$se_logrr else NA_real_,
      tmle_rd = if (!is.null(tm)) tm$rd else NA_real_,
      tmle_se_rd = if (!is.null(tm)) tm$se_rd else NA_real_
    ),
    spd = spd_dt[, .(
      scenario_id = scn_row$scenario_id,
      mis_spec = scn_row$mis_spec,
      delta_dgm = scn_row$delta_dgm,
      truth = truth,
      replicate = rep_id,
      interval, tstart, tstop, t_mid,
      gamma_hat, se, HR_1SD, CI_low, CI_high,
      Gamma_hat
    )],
    diag = tip$diag[, .(
      scenario_id = scn_row$scenario_id,
      mis_spec = scn_row$mis_spec,
      delta_dgm = scn_row$delta_dgm,
      truth = truth,
      replicate = rep_id,
      t, N, ESS, rESS, tail_share, tipped
    )],
    delta_star = data.table::rbindlist(delta_star_rows, fill = TRUE),
    sens_curve = if (SAVE_SENS_CURVES) data.table::rbindlist(sens_rows, fill = TRUE) else NULL
  )
}

# ---- run all ----
jobs <- data.table::CJ(job_id = seq_len(nrow(scn_grid)), replicate = seq_len(B), unique = TRUE)

# Use plain data.frames inside multisession workers. On some Windows/RStudio
# setups, data.table methods are not attached in newly spawned workers before
# globals are restored, so one-column row extraction like DT[i] can be
# interpreted as data.frame column extraction. The explicit data.frame row
# indexing below avoids this platform-specific failure while preserving all
# analysis inputs.
jobs_df <- as.data.frame(jobs)
scn_grid_df <- as.data.frame(scn_grid)

# run in chunks; each job maps to one scn_row + rep
res_list <- future.apply::future_lapply(
  seq_len(nrow(jobs_df)),
  function(i) {
    j <- jobs_df[i, , drop = FALSE]
    scn_row <- scn_grid_df[as.integer(j$job_id), , drop = FALSE]
    out <- run_one(scn_row, rep_id = as.integer(j$replicate))
    if (i %% 50 == 0) cfg$log_line(log_file, sprintf("[%s] progress %d / %d\n", format(Sys.time(), "%H:%M:%S"), i, nrow(jobs_df)))
    out
  },
  future.packages = c("data.table", "survival")
)

rep_dt <- data.table::rbindlist(lapply(res_list, `[[`, "rep"), fill = TRUE)
spd_dt <- data.table::rbindlist(lapply(res_list, `[[`, "spd"), fill = TRUE)
diag_dt<- data.table::rbindlist(lapply(res_list, `[[`, "diag"), fill = TRUE)
ds_dt  <- data.table::rbindlist(lapply(res_list, `[[`, "delta_star"), fill = TRUE)

write_dt_csv(rep_dt, file.path(raw_dir, "replicate_results_paperc.csv"))
write_dt_csv(spd_dt, file.path(raw_dir, "spd_curves_paperc.csv"))
write_dt_csv(diag_dt, file.path(raw_dir, "weight_diagnostics_paperc.csv"))
write_dt_csv(ds_dt, file.path(raw_dir, "delta_star_paperc.csv"))

if (SAVE_SENS_CURVES) {
  sc_dt <- data.table::rbindlist(lapply(res_list, `[[`, "sens_curve"), fill = TRUE)
  if (!is.null(sc_dt) && nrow(sc_dt) > 0) {
    write_dt_csv(sc_dt, file.path(raw_dir, "sensitivity_curves_paperc.csv"))
  }
}

# ---- performance summaries for baseline estimators (risk scale) ----
# join truth by delta_dgm (truth is identical across scenario axes, but we keep keys explicit)
if (!is.null(mc_truth) && nrow(mc_truth) > 0) {
  rep2 <- merge(rep_dt, mc_truth[, .(delta_dgm, logrr_true = logrr, rd_true = rd)], by = "delta_dgm", all.x = TRUE)
} else {
  rep2 <- rep_dt
  rep2[, `:=`(logrr_true = ifelse(truth == "non_null", beta_true_pp, 0), rd_true = NA_real_)]
}

summ_base <- function(dt, est_col, se_col, true_col) {
  out <- dt[, {
    est <- get(est_col); se <- get(se_col); tr <- get(true_col)
    ok <- is.finite(est) & is.finite(se) & is.finite(tr)
    if (!any(ok)) return(.(
      n_ok = 0L, bias = NA_real_, rmse = NA_real_, cover = NA_real_, sign_error = NA_real_, ci_excl0 = NA_real_
    ))
    z  <- stats::qnorm(0.975)
    lo <- est - z * se
    hi <- est + z * se
    .(
      n_ok = sum(ok),
      bias = mean(est[ok] - tr[ok]),
      rmse = sqrt(mean((est[ok] - tr[ok])^2)),
      cover = mean(lo[ok] <= tr[ok] & tr[ok] <= hi[ok]),
      sign_error = mean(sign(est[ok]) != sign(tr[ok]) & sign(tr[ok]) != 0),
      ci_excl0 = mean(!(lo[ok] <= 0 & 0 <= hi[ok]))
    )
  }, by = .(scenario_id, mis_spec, delta_dgm, truth, axisA, axisB, rho_meas)]
  out
}

ipcw_sum <- summ_base(rep2, "ipcw_logrr", "ipcw_se_logrr", "logrr_true")
write_dt_csv(ipcw_sum, file.path(cfg$out_dir, "perf_summary_paperc_ipcw.csv"))

if (RUN_DR) {
  dr_sum <- summ_base(rep2[is.finite(dr_logrr)], "dr_logrr", "dr_se_logrr", "logrr_true")
  # annotate method_Q if ML used
  dr_sum[, method_Q := if (RUN_ML) paste0("ML_", ML_METHOD) else "glm"]
  write_dt_csv(dr_sum, file.path(cfg$out_dir, "perf_summary_paperc_dr.csv"))
}
if (RUN_TMLE) {
  tm_sum <- summ_base(rep2[is.finite(tmle_logrr)], "tmle_logrr", "tmle_se_logrr", "logrr_true")
  write_dt_csv(tm_sum, file.path(cfg$out_dir, "perf_summary_paperc_tmle.csv"))
}

# ---- delta-star summary (median/IQR) ----
ds_sum <- ds_dt[, .(
  n = .N,
  n_finite = sum(is.finite(delta_star_logrr)),
  med_delta_star_logrr = stats::median(delta_star_logrr[is.finite(delta_star_logrr)], na.rm = TRUE),
  q25_delta_star_logrr = stats::quantile(delta_star_logrr[is.finite(delta_star_logrr)], 0.25, na.rm = TRUE, type = 7),
  q75_delta_star_logrr = stats::quantile(delta_star_logrr[is.finite(delta_star_logrr)], 0.75, na.rm = TRUE, type = 7),
  p_inf = mean(is.infinite(delta_star_logrr), na.rm = TRUE)
), by = .(scenario_id, mis_spec, delta_dgm, truth, method, delta_shape)]

write_dt_csv(ds_sum, file.path(cfg$out_dir, "perf_summary_paperc_delta_star.csv"))

cfg$write_session_info(prefix = paste0("paperc_", TAG))
cfg$log_line(log_file, sprintf("[%s] DONE\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Done. Outputs under: ", cfg$out_dir, "\n", sep = "")
