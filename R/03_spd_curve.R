# R/03_spd_curve.R
# ------------------------------------------------------------
# SPD(t) estimation (piecewise slope per time interval)
# SPD(t) is the log-HR for deviation per 1 SD higher prognostic score at time t.
#
# We standardize z within each interval (risk-set proxy) to achieve scale-invariance.
# ------------------------------------------------------------

standardize_by_interval <- function(dt, z_col = "z_obs", interval_col = "interval") {
  z <- dt[[z_col]]
  # standardize within interval (only finite values)
  dt[, z_std := {
    zz <- get(z_col)
    m <- mean(zz, na.rm = TRUE)
    s <- stats::sd(zz, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, .N) else (zz - m) / s
  }, by = interval_col]
  dt[]
}

estimate_spd_piecewise <- function(long_dt,
                                   breaks,
                                   z_col = "z_obs",
                                   dev_col = "dev",
                                   id_col = "id",
                                   robust = TRUE) {
  dt <- data.table::copy(long_dt)
  dt[, tstart := as.numeric(tstart)]
  dt[, tstop  := as.numeric(tstop)]
  eps_dt <- 1e-6
  dt <- dt[is.finite(tstart) & is.finite(tstop) & (tstop - tstart) > eps_dt]

  # define interval from breaks using tstart
  dt[, interval := findInterval(tstart, vec = breaks, rightmost.closed = TRUE)]
  dt[interval < 1, interval := 1L]
  dt[, interval := as.integer(interval)]
  dt[, interval_f := factor(interval)]

  dt <- standardize_by_interval(dt, z_col = z_col, interval_col = "interval")

  # separate slope per interval, allow baseline hazard to vary by interval
  # - cluster(id) for robust SE by subject if requested
  fml <- stats::as.formula(paste0("survival::Surv(tstart, tstop, ", dev_col, ") ~ ",
                                 "0 + z_std:interval_f + strata(interval_f)",
                                 if (robust) paste0(" + cluster(", id_col, ")") else ""))

  fit <- tryCatch(survival::coxph(fml, data = dt), error = function(e) NULL)
  if (is.null(fit)) {
    # If the deviation Cox model fails (e.g., rare numerical zero-length intervals),
    # return an NA curve for the intervals present so the simulation can continue.
    intervals <- sort(unique(dt$interval))
    out <- data.table::data.table(interval = as.integer(intervals), gamma_hat = NA_real_, se = NA_real_)
    out[, HR_1SD := NA_real_]
    out[, `:=`(CI_low = NA_real_, CI_high = NA_real_)]
    out <- out[order(interval)]
    out[, tstart := breaks[pmax(interval, 1L)]]
    out[, tstop  := breaks[pmin(interval + 1L, length(breaks))]]
    out[, t_mid  := (tstart + tstop) / 2]
    return(list(fit = NULL, spd = out, dt_used = dt))
  }
  coefs <- stats::coef(fit)
  vc    <- stats::vcov(fit)
  nm    <- names(coefs)

  # parse interval index from coefficient names like "z_std:interval_f1"
  k <- as.integer(gsub(".*interval_f", "", nm))
  se <- sqrt(diag(vc))

  out <- data.table::data.table(
    interval = k,
    gamma_hat = as.numeric(coefs),
    se = as.numeric(se)
  )
  out[, HR_1SD := exp(gamma_hat)]
  z <- stats::qnorm(0.975)
  out[, `:=`(
    CI_low = exp(gamma_hat - z * se),
    CI_high= exp(gamma_hat + z * se)
  )]

  # add time midpoints
  out <- out[order(interval)]
  out[, tstart := breaks[pmax(interval, 1L)]]
  out[, tstop  := breaks[pmin(interval + 1L, length(breaks))]]
  out[, t_mid  := (tstart + tstop) / 2]

  list(fit = fit, spd = out, dt_used = dt)
}

compute_cum_pressure <- function(spd_dt) {
  dt <- data.table::copy(spd_dt)
  dt[, dt_len := (tstop - tstart)]
  dt[, Gamma_hat := cumsum(gamma_hat * dt_len)]
  dt[]
}