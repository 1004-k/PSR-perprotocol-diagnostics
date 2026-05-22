# R/10_dgp_paperc.R
# ------------------------------------------------------------
# Paper C DGP extension: adds an unmeasured prognostic component to induce
# residual non-ignorability after conditioning on the observed prognostic score.
#
# Design philosophy (acceptance-risk minimized):
# - Keep the original Paper B DGP structure and axes A/B/C intact.
# - Add a latent "unmeasured" component z_unm(t) that affects deviation AND outcome,
#   but is not captured by the observed prognostic score z_obs(t).
# - Axis D controls the time-pattern of the unmeasured component's effect on outcome hazard.
#
# Output matches simulate_one_dataset(): base, long_dt, pp_dt, grid
# ------------------------------------------------------------

# A helper that returns the coefficient of the unmeasured component on the outcome hazard
# (piecewise constant on the dt grid).
.make_theta_unm <- function(delta_dgm,
                            tstart_vec,
                            delta_mag = 0.0,
                            t0 = NULL) {
  delta_dgm <- match.arg(delta_dgm, c("D1_zero", "D2_constant", "D3_late_only", "D4_sign_change"))
  if (is.null(t0)) t0 <- 0.7 * max(tstart_vec)

  if (delta_dgm == "D1_zero") {
    return(rep(0, length(tstart_vec)))
  }
  if (delta_dgm == "D2_constant") {
    return(rep(delta_mag, length(tstart_vec)))
  }
  if (delta_dgm == "D3_late_only") {
    return(ifelse(tstart_vec < t0, 0, delta_mag))
  }
  # D4: sign change (supplement)
  ifelse(tstart_vec < t0, -delta_mag, delta_mag)
}

simulate_one_dataset_paperc <- function(N = 1000,
                                        t_max = 5,
                                        dt = 0.25,
                                        beta_true = -0.2231436,  # numeric string recommended in bash
                                        # hazards
                                        lambda_event_base = 0.10,
                                        lambda_dev_base   = 0.25,
                                        # measured prognosis effect on outcome
                                        theta_meas = 0.50,
                                        # Axis A/B/C (same semantics as Paper B)
                                        axisA = c("flat","increasing","late_surge"),
                                        gamma_max = 1.0,
                                        t0_gamma = NULL,
                                        axisB = c("none","threshold_jump"),
                                        thresh_c = 1.0,
                                        thresh_mult = 3.0,
                                        rho_meas = 1.0,
                                        # latent dynamics (AR(1))
                                        rho_z = 0.80,
                                        # NEW: Axis D (residual MNAR pattern) + magnitude
                                        delta_dgm = c("D1_zero","D2_constant","D3_late_only","D4_sign_change"),
                                        delta_mag = 0.0,
                                        t0_delta = NULL,
                                        # correlation between measured and unmeasured components (0 = independent)
                                        corr_meas_unm = 0.0,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  axisA <- match.arg(axisA, c("flat","increasing","late_surge"))
  axisB <- match.arg(axisB, c("none","threshold_jump"))
  delta_dgm <- match.arg(delta_dgm, c("D1_zero","D2_constant","D3_late_only","D4_sign_change"))

  # --- baseline covariates (copy from Paper B for continuity) ---
  base <- data.table::data.table(
    id  = 1:N,
    age = pmin(pmax(stats::rnorm(N, 62, 10), 40), 90),
    sex = stats::rbinom(N, 1, 0.45),
    bmi = pmin(pmax(stats::rnorm(N, 31, 6), 18), 60),
    egfr= pmin(pmax(stats::rnorm(N, 75, 18), 15), 120),
    util= stats::rgamma(N, shape = 2, rate = 0.5),
    gall= stats::rbinom(N, 1, stats::plogis(-2 + 0.02 * (stats::rnorm(N, 62, 10) - 60)))
  )

  # treatment assignment (baseline)
  lin <- with(base,
              -1.0 + 0.02*(age-60) + 0.08*(bmi-30) - 0.015*(egfr-75) + 0.25*util + 0.7*gall)
  ps <- stats::plogis(lin)
  base[, A0 := stats::rbinom(.N, 1, ps)]

  # baseline prognosis drivers (measured + unmeasured)
  lp_prog0 <- with(base,
                   0.02*(age-60) + 0.10*(bmi-30) - 0.01*(egfr-75) + 0.25*gall + 0.15*util)
  base[, lp_prog0 := lp_prog0]
  base[, z0_meas := as.numeric(scale(lp_prog0))]

  # baseline unmeasured component (optionally correlated with measured)
  if (abs(corr_meas_unm) < 1e-8) {
    base[, z0_unm := stats::rnorm(.N)]
  } else {
    # create correlated normal: z_unm = r*z_meas + sqrt(1-r^2)*eps
    r <- max(min(corr_meas_unm, 0.99), -0.99)
    base[, z0_unm := r * z0_meas + sqrt(1 - r^2) * stats::rnorm(.N)]
  }

  # time grid
  t_grid <- seq(0, t_max, by = dt)
  if (tail(t_grid, 1) < t_max) t_grid <- c(t_grid, t_max)
  tstart_vec <- t_grid[-length(t_grid)]
  tstop_vec  <- t_grid[-1L]
  K <- length(tstart_vec)

  gamma_vec <- make_gamma_path(axisA, tstart_vec, gamma_max = gamma_max, t0 = t0_gamma)
  theta_unm_vec <- .make_theta_unm(delta_dgm, tstart_vec, delta_mag = delta_mag, t0 = t0_delta)

  # allocate containers
  long_list <- vector("list", N)

  for (i in seq_len(N)) {
    id <- base$id[i]
    A0 <- base$A0[i]

    # latent processes
    z_meas <- numeric(K)
    z_unm  <- numeric(K)
    z_meas[1] <- base$z0_meas[i]
    z_unm[1]  <- base$z0_unm[i]

    if (K >= 2) {
      for (k in 2:K) {
        z_meas[k] <- rho_z * z_meas[k-1] + sqrt(1 - rho_z^2) * stats::rnorm(1)
        z_unm[k]  <- rho_z * z_unm[k-1]  + sqrt(1 - rho_z^2) * stats::rnorm(1)
      }
    }

    # observed measured score only (measurement error on z_meas)
    z_obs <- rho_meas * z_meas + sqrt(1 - rho_meas^2) * stats::rnorm(K)

    # deviation hazard uses total latent prognosis
    z_tot <- z_meas + z_unm

    alive <- TRUE
    deviated <- FALSE
    rows <- list()

    for (k in seq_len(K)) {
      if (!alive || deviated) break

      t0k <- tstart_vec[k]
      t1k <- tstop_vec[k]

      # event hazard depends on measured + unmeasured with time-patterned coefficient
      haz_y <- lambda_event_base * exp(theta_meas * z_meas[k] + theta_unm_vec[k] * z_unm[k] + beta_true * A0)

      # deviation hazard depends on total prognosis and selection pressure path
      haz_d <- lambda_dev_base * exp(gamma_vec[k] * z_tot[k])
      if (axisB == "threshold_jump" && z_tot[k] > thresh_c) {
        haz_d <- haz_d * thresh_mult
      }

      eps_dt <- 1e-6
      ty <- max(stats::rexp(1, rate = haz_y), eps_dt)
      td <- max(stats::rexp(1, rate = haz_d), eps_dt)
      t_event <- min(ty, td)

      if (t_event < dt) {
        t_end <- t0k + t_event
        if ((t_end - t0k) <= eps_dt) t_end <- t0k + eps_dt
        dev <- as.integer(td < ty)
        ev  <- as.integer(ty <= td)
        rows[[length(rows)+1L]] <- data.table::data.table(
          id = id,
          tstart = t0k,
          tstop  = t_end,
          interval = k,
          A0 = A0,
          z_meas = z_meas[k],
          z_unm  = z_unm[k],
          z_tot  = z_tot[k],
          z_obs  = z_obs[k],
          gamma_true = gamma_vec[k],
          theta_unm_true = theta_unm_vec[k],
          dev   = dev,
          event = ev
        )
        if (dev == 1) deviated <- TRUE
        if (ev  == 1) alive    <- FALSE
        break
      } else {
        rows[[length(rows)+1L]] <- data.table::data.table(
          id = id,
          tstart = t0k,
          tstop  = t1k,
          interval = k,
          A0 = A0,
          z_meas = z_meas[k],
          z_unm  = z_unm[k],
          z_tot  = z_tot[k],
          z_obs  = z_obs[k],
          gamma_true = gamma_vec[k],
          theta_unm_true = theta_unm_vec[k],
          dev   = 0L,
          event = 0L
        )
      }
    }

    if (length(rows) == 0) {
      rows[[1]] <- data.table::data.table(
        id = id, tstart = 0, tstop = dt, interval = 1L,
        A0 = A0,
        z_meas = base$z0_meas[i],
        z_unm  = base$z0_unm[i],
        z_tot  = base$z0_meas[i] + base$z0_unm[i],
        z_obs  = base$z0_meas[i],
        gamma_true = gamma_vec[1],
        theta_unm_true = theta_unm_vec[1],
        dev = 0L, event = 0L
      )
    }

    long_list[[i]] <- data.table::rbindlist(rows)
  }

  long_dt <- data.table::rbindlist(long_list)

  # per-protocol summary per id (censor at deviation)
  dev_time <- long_dt[dev == 1, .(dev_time = min(tstop)), by = id]
  ev_time  <- long_dt[event == 1, .(event_time = min(tstop)), by = id]

  pp_dt <- data.table::copy(base)[, .(id, A0, age, sex, bmi, egfr, util, gall)]

  pp_dt <- merge(pp_dt, dev_time, by = "id", all.x = TRUE)
  pp_dt <- merge(pp_dt, ev_time,  by = "id", all.x = TRUE)

  pp_dt[is.na(dev_time), dev_time := Inf]
  pp_dt[is.na(event_time), event_time := Inf]

  pp_dt[, time_pp := pmin(event_time, dev_time, t_max)]
  pp_dt[, delta_pp := as.integer(event_time <= dev_time & event_time <= t_max)]
  pp_dt[, dev_ind  := as.integer(dev_time < pmin(event_time, t_max))]

  # baseline observed z for convenience (risk-set history is in long_dt)
  z0_obs <- long_dt[interval == 1, .(id, z0_obs = z_obs)]
  pp_dt <- merge(pp_dt, z0_obs, by = "id", all.x = TRUE)

  list(
    base   = base,
    long_dt= long_dt,
    pp_dt  = pp_dt,
    grid   = data.table::data.table(
      k = seq_len(K),
      tstart = tstart_vec,
      tstop  = tstop_vec,
      gamma_true = gamma_vec,
      theta_unm_true = theta_unm_vec
    )
  )
}

# ------------------------------------------------------------
# Monte Carlo truth for risk-scale estimand under per-protocol adherence
# (no deviation/censoring), given the same covariate + latent process distribution.
#
# We compute risk deterministically given piecewise-constant hazards:
#   H = sum_k lambda_k * dt_k,  risk = 1 - exp(-H)
# ------------------------------------------------------------

mc_truth_risk_paperc <- function(N_truth = 200000,
                                 t_max = 5,
                                 dt = 0.25,
                                 t_cut = 5,
                                 beta_true = -0.2231436,
                                 theta_meas = 0.50,
                                 delta_dgm = c("D1_zero","D2_constant","D3_late_only","D4_sign_change"),
                                 delta_mag = 0.0,
                                 t0_delta = NULL,
                                 rho_z = 0.80,
                                 corr_meas_unm = 0.0,
                                 lambda_event_base = 0.10,
                                 seed = 1L) {
  delta_dgm <- match.arg(delta_dgm, c("D1_zero","D2_constant","D3_late_only","D4_sign_change"))
  if (!is.null(seed)) set.seed(seed)

  # Generate baseline covariates (same as simulate_one_dataset_paperc)
  base <- data.table::data.table(
    id  = 1:N_truth,
    age = pmin(pmax(stats::rnorm(N_truth, 62, 10), 40), 90),
    sex = stats::rbinom(N_truth, 1, 0.45),
    bmi = pmin(pmax(stats::rnorm(N_truth, 31, 6), 18), 60),
    egfr= pmin(pmax(stats::rnorm(N_truth, 75, 18), 15), 120),
    util= stats::rgamma(N_truth, shape = 2, rate = 0.5),
    gall= stats::rbinom(N_truth, 1, stats::plogis(-2 + 0.02 * (stats::rnorm(N_truth, 62, 10) - 60)))
  )

  lp_prog0 <- with(base,
                   0.02*(age-60) + 0.10*(bmi-30) - 0.01*(egfr-75) + 0.25*gall + 0.15*util)
  base[, z0_meas := as.numeric(scale(lp_prog0))]

  if (abs(corr_meas_unm) < 1e-8) {
    base[, z0_unm := stats::rnorm(.N)]
  } else {
    r <- max(min(corr_meas_unm, 0.99), -0.99)
    base[, z0_unm := r * z0_meas + sqrt(1 - r^2) * stats::rnorm(.N)]
  }

  # time grid
  t_grid <- seq(0, t_max, by = dt)
  if (tail(t_grid, 1) < t_max) t_grid <- c(t_grid, t_max)
  tstart_vec <- t_grid[-length(t_grid)]
  tstop_vec  <- t_grid[-1L]
  dt_vec <- pmin(tstop_vec, t_cut) - pmin(tstart_vec, t_cut)
  keep <- dt_vec > 0
  tstart_vec <- tstart_vec[keep]
  dt_vec <- dt_vec[keep]

  theta_unm_vec <- .make_theta_unm(delta_dgm, tstart_vec, delta_mag = delta_mag, t0 = t0_delta)

  K <- length(tstart_vec)

  # simulate latent processes for all subjects (AR(1))
  z_meas <- matrix(0, nrow = N_truth, ncol = K)
  z_unm  <- matrix(0, nrow = N_truth, ncol = K)
  z_meas[,1] <- base$z0_meas
  z_unm[,1]  <- base$z0_unm
  if (K >= 2) {
    for (k in 2:K) {
      z_meas[,k] <- rho_z * z_meas[,k-1] + sqrt(1 - rho_z^2) * stats::rnorm(N_truth)
      z_unm[,k]  <- rho_z * z_unm[,k-1]  + sqrt(1 - rho_z^2) * stats::rnorm(N_truth)
    }
  }

  # cumulative hazard under A=0 and A=1
  # lambda_k(a) = lambda0 * exp(theta_meas*z_meas + theta_unm*z_unm + beta_true * a)
  lin0 <- theta_meas * z_meas + (matrix(theta_unm_vec, nrow = N_truth, ncol = K, byrow = TRUE) * z_unm)
  H0 <- rowSums(lambda_event_base * exp(lin0 + 0) * matrix(dt_vec, nrow = N_truth, ncol = K, byrow = TRUE))
  H1 <- rowSums(lambda_event_base * exp(lin0 + beta_true) * matrix(dt_vec, nrow = N_truth, ncol = K, byrow = TRUE))

  risk0 <- mean(1 - exp(-H0))
  risk1 <- mean(1 - exp(-H1))

  rd <- risk1 - risk0
  logrr <- log(risk1 / pmax(risk0, 1e-9))

  data.table::data.table(
    t_cut = t_cut,
    N_truth = N_truth,
    beta_true = beta_true,
    delta_dgm = delta_dgm,
    delta_mag = delta_mag,
    risk1 = risk1,
    risk0 = risk0,
    rd = rd,
    logrr = logrr
  )
}
