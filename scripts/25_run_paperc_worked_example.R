#!/usr/bin/env Rscript
# scripts/25_run_paperc_worked_example.R
# ------------------------------------------------------------
# Worked example generator for Paper C:
# - Runs ONE dataset for a chosen scenario and delta_dgm
# - Saves:
#   * SPD(t) curve
#   * rESS(t), tail-share(t) curves
#   * Fragility curve logRR(delta) for IPCW and DR (if enabled)
# ------------------------------------------------------------

source("R/00_utils.R")
cfg <- init_project(n_cores = 1L, seed = 2026L, out_dir = Sys.getenv("OUT_DIR", "output_paperc_worked"))
cfg$require_pkgs(c("data.table", "survival"))

source("R/01_scenarios.R")
source("R/02_dgp_simulate.R")
source("R/03_spd_curve.R")
source("R/04_pp_ipcw.R")
source("R/05_tipping.R")
source("R/09_dr_aipw_ml.R")
source("R/10_dgp_paperc.R")
source("R/11_paperc_sensitivity.R")
source("R/07_plotting.R")

sid <- Sys.getenv("SCENARIO_ID", "S09")
mis_spec <- as.integer(Sys.getenv("MIS_SPEC", "0"))
delta_dgm <- Sys.getenv("DELTA_DGM", "D2_constant")
truth <- Sys.getenv("TRUTH", "non_null")
beta_true_pp <- as.numeric(Sys.getenv("BETA_TRUE_PP", "-0.2231436"))

N <- as.integer(Sys.getenv("N", "2000"))
t_max <- as.numeric(Sys.getenv("T_MAX","5"))
dt <- as.numeric(Sys.getenv("DT","0.25"))
T_CUT <- as.numeric(Sys.getenv("T_CUT", as.character(t_max)))

DELTA_MAG <- as.numeric(Sys.getenv("DELTA_MAG", "0.6"))
T0_DELTA <- as.numeric(Sys.getenv("T0_DELTA", "3.5"))
CORR_MEAS_UNM <- as.numeric(Sys.getenv("CORR_MEAS_UNM","0.0"))

BREAK_DT <- as.numeric(Sys.getenv("BREAK_DT","1.0"))
breaks <- seq(0, t_max, by = BREAK_DT); if (tail(breaks,1) < t_max) breaks <- c(breaks, t_max)

DELTA_MAX <- as.numeric(Sys.getenv("DELTA_MAX","2.0"))
DELTA_STEP <- as.numeric(Sys.getenv("DELTA_STEP","0.1"))
delta_grid <- seq(0, DELTA_MAX, by = DELTA_STEP)

grid <- make_scenario_grid()
scn <- grid[scenario_id == sid][1]
if (nrow(scn) == 0) stop("Unknown SCENARIO_ID: ", sid)

beta <- if (truth == "non_null") beta_true_pp else 0

sim <- simulate_one_dataset_paperc(
  N = N, t_max = t_max, dt = dt,
  beta_true = beta,
  axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
  delta_dgm = delta_dgm, delta_mag = DELTA_MAG, t0_delta = T0_DELTA,
  corr_meas_unm = CORR_MEAS_UNM,
  seed = 2026L
)

long_dt <- sim$long_dt
pp_dt <- sim$pp_dt

# SPD(t)
spd_fit <- estimate_spd_piecewise(long_dt, breaks = breaks, z_col = "z_obs", dev_col = "dev", id_col = "id", robust = TRUE)
spd_dt <- compute_cum_pressure(spd_fit$spd)

# IPCW
cens_fml <- if (mis_spec == 1L) {
  survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi
} else {
  survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall
}
fit_cens <- survival::coxph(cens_fml, data = pp_dt, x = TRUE, model = TRUE)
fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
lp <- stats::predict(fit_cens, type = "lp")

ipcw_risk <- fit_ipcw_risk(pp_dt, t_cut = T_CUT, fit_cens = fit_cens, lp = lp, w_floor = 1e-3)

# weight diagnostics
t_eval <- sort(unique(long_dt$tstart))
weights_dt <- data.table::CJ(id = pp_dt$id, t = t_eval)
weights_dt <- merge(weights_dt, data.table::data.table(id = pp_dt$id, lp = lp), by = "id", all.x = TRUE)
weights_dt[, S := get_S_from_cox(fit_cens, time_vec = t, lp_vec = lp)]
weights_dt[, w := 1 / pmax(S, 1e-3)]
diag_dt <- compute_weight_diagnostics(long_dt, weights_dt, time_col = "tstart", id_col = "id", w_col = "w", q_tail = 0.99)

# DR (glm only for worked example)
dr <- fit_dr_risk_rescue(pp_dt, t_cut = T_CUT, mis_spec = mis_spec, use_ml_Q = FALSE, ml_method = "glmnet", crossfit = FALSE, cf_folds = 2L, seed = 999)

# Sens curves
curve_ipcw <- anchored_theta_curve(pp_dt, t_cut = T_CUT, breaks = breaks,
                                   delta_grid = delta_grid, shape = "constant", t0 = 0.7*t_max,
                                   theta0_logrr = ipcw_risk$logrr, theta0_rd = ipcw_risk$rd)
curve_dr <- anchored_theta_curve(pp_dt, t_cut = T_CUT, breaks = breaks,
                                 delta_grid = delta_grid, shape = "constant", t0 = 0.7*t_max,
                                 theta0_logrr = dr$logrr, theta0_rd = dr$rd)

# Save raw artifacts
raw_dir <- file.path(cfg$out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(spd_dt, file.path(raw_dir, "worked_spd_curve.csv"))
data.table::fwrite(diag_dt, file.path(raw_dir, "worked_weight_diag.csv"))
data.table::fwrite(curve_ipcw, file.path(raw_dir, "worked_fragility_ipcw.csv"))
data.table::fwrite(curve_dr, file.path(raw_dir, "worked_fragility_dr.csv"))

# Figures
fig_dir <- file.path(cfg$out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Bundle
fig1 <- file.path(fig_dir, "Worked_bundle.pdf")
open_pdf(fig1, width = 10.8, height = 4.2, pointsize = 12)
par(mfrow = c(1,3), mar = c(4.2,4.4,2.4,1.0), oma = c(0.6,0.6,1.0,0.4))
plot(spd_dt$t_mid, spd_dt$gamma_hat, type = "l", lwd = 3, xlab = "Time", ylab = "SPD(t)", main = "SPD(t)")
abline(h = 0, lty = 3, col = "grey60")
plot(diag_dt$t, diag_dt$rESS, type = "l", lwd = 3, xlab = "Time", ylab = "rESS(t)", main = "rESS(t)")
abline(h = 0.25, lty = 3, col = "grey60")
plot(diag_dt$t, diag_dt$tail_share, type = "l", lwd = 3, xlab = "Time", ylab = "Tail-share(t)", main = "Tail-share(t)")
abline(h = 0.10, lty = 3, col = "grey60")
# Figure title is provided in the manuscript legend, not in the graphic file.
close_device()

# Fragility curve
fig2 <- file.path(fig_dir, "Worked_fragility.pdf")
open_pdf(fig2, width = 8.2, height = 4.2, pointsize = 12)
par(mfrow = c(1,1), mar = c(4.2,4.4,2.4,1.0))
plot(curve_ipcw$delta, curve_ipcw$logrr, type = "l", lwd = 3,
     xlab = expression(delta), ylab = "Anchored logRR(delta)",
     main = "")
lines(curve_dr$delta, curve_dr$logrr, lwd = 3, lty = 2)
abline(h = 0, lty = 3, col = "grey60")
legend("topright", legend = c("IPCW", "DR"), lty = c(1,2), lwd = 3, bty = "n")
close_device()

cat("Worked example saved under: ", cfg$out_dir, "\n", sep = "")
