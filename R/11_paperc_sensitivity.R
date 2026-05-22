# R/11_paperc_sensitivity.R
# ------------------------------------------------------------
# Paper C: residual non-ignorability sensitivity (delta) and tipping summaries.
#
# Core idea:
# - We compute how the per-protocol risk contrast would MOVE as delta increases,
#   using a simple, transparent "virtual hazard shift" on the post-deviation period.
# - Then we ANCHOR that movement to any baseline estimator theta^0:
#     theta^(delta) = theta^0 + (theta_virtual(delta) - theta_virtual(0))
#   so that theta^(0) exactly matches the baseline estimator.
#
# Delta parameterization:
# - delta(t) is a (piecewise-constant) log hazard multiplier applied AFTER deviation,
#   on the remaining event hazard up to t_cut.
# - In terms of conditional survival beyond deviation time u:
#     S_delta(t_cut | u) = [S(t_cut)/S(u)]^{exp(delta_k)}  for the interval k containing u.
# ------------------------------------------------------------

# Build a piecewise-constant delta profile on the SPD breaks
make_delta_profile <- function(breaks,
                               shape = c("constant", "late_only", "sign_change"),
                               delta = 0.0,
                               t0 = NULL) {
  shape <- match.arg(shape, c("constant", "late_only", "sign_change"))
  tstart <- breaks[-length(breaks)]
  if (is.null(t0)) t0 <- 0.7 * max(breaks)

  if (shape == "constant") {
    return(rep(delta, length(tstart)))
  }
  if (shape == "late_only") {
    return(ifelse(tstart < t0, 0.0, delta))
  }
  # sign_change
  ifelse(tstart < t0, -delta, delta)
}

# Safe lookup of delta at a censoring time u
.delta_at_time <- function(u, breaks, delta_vec) {
  k <- findInterval(u, vec = breaks, rightmost.closed = TRUE)
  k <- pmax(pmin(k, length(delta_vec)), 1L)
  delta_vec[k]
}

# Virtual hazard-shift engine (model-based movement only)
# Returns marginal risks under A=1 and A=0 as well as RD and logRR.
virtual_risk_under_delta <- function(pp_dt,
                                     t_cut,
                                     breaks,
                                     delta_vec,
                                     outcome_formula = survival::Surv(time_pp, delta_pp) ~ A0 + age + sex + bmi + egfr + util + gall,
                                     w_floor = 1e-6) {
  d <- data.table::copy(pp_dt)

  # Fit Cox on observed per-protocol data (censored at deviation by design)
  fit <- tryCatch(
    survival::coxph(outcome_formula, data = d, x = TRUE, model = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(risk1 = NA_real_, risk0 = NA_real_, rd = NA_real_, logrr = NA_real_))
  }
  fit$bh0 <- survival::basehaz(fit, centered = FALSE)

  # Construct design matrices for A=0 and A=1
  # (We build from the model frame to match factor handling if any.)
  mf <- stats::model.frame(fit)
  # Extract the terms object. Use delete.response() for model.matrix();
  # otherwise some R versions include the Surv response and/or an intercept
  # column, which has no corresponding Cox coefficient and can propagate NAs.
  tt <- stats::delete.response(stats::terms(fit))

  # Build a helper to get lp for counterfactual A0.
  get_lp_for_A <- function(a) {
    mf2 <- mf
    if ("A0" %in% names(mf2)) mf2$A0 <- as.integer(a)
    X <- stats::model.matrix(tt, data = mf2)

    # Cox models have no intercept. Drop it if model.matrix() creates one.
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
    }

    b <- stats::coef(fit)
    common <- intersect(colnames(X), names(b))
    if (length(common) == 0L) {
      return(rep(0, nrow(X)))
    }

    # Missing coefficient columns are treated as zero contribution. This keeps
    # counterfactual prediction stable when a term is absent because of
    # rank-deficiency in a small quickcheck run.
    lp <- rep(0, nrow(X))
    lp <- lp + as.numeric(X[, common, drop = FALSE] %*% b[common])
    lp
  }

  lp1 <- get_lp_for_A(1L)
  lp0 <- get_lp_for_A(0L)

  # censoring (deviation) time used as missingness time
  u <- pmin(d$dev_time, t_cut)
  censored <- is.finite(d$dev_time) & (d$dev_time < t_cut) & (d$delta_pp == 0L)

  # helper: survival prediction at a vector of times given lp
  S_at <- function(time_vec, lp_vec) {
    bh <- fit$bh0
    bh <- bh[order(bh$time), ]
    idx <- findInterval(time_vec, bh$time)
    H0 <- ifelse(idx == 0, 0, bh$hazard[idx])
    S <- exp(- H0 * exp(lp_vec))
    pmax(S, w_floor)
  }

  # predicted survival at u and at t_cut under A=1 and A=0
  S1_u <- S_at(u, lp1); S1_t <- S_at(rep(t_cut, nrow(d)), lp1)
  S0_u <- S_at(u, lp0); S0_t <- S_at(rep(t_cut, nrow(d)), lp0)

  # delta at each subject's censoring time
  delta_i <- vapply(u, .delta_at_time, numeric(1), breaks = breaks, delta_vec = delta_vec)
  mult_i <- exp(delta_i)

  # apply virtual shift only to censored subjects; others keep baseline S(t_cut)
  adj_surv <- function(S_u, S_t) {
    S_adj <- S_t
    if (any(censored)) {
      S_cond <- pmin(pmax(S_t / pmax(S_u, w_floor), w_floor), 1)
      S_cond_adj <- S_cond ^ mult_i
      S_adj[censored] <- S_u[censored] * S_cond_adj[censored]
    }
    pmin(pmax(S_adj, w_floor), 1)
  }

  S1_adj <- adj_surv(S1_u, S1_t)
  S0_adj <- adj_surv(S0_u, S0_t)

  risk1 <- mean(1 - S1_adj)
  risk0 <- mean(1 - S0_adj)

  rd <- risk1 - risk0
  logrr <- log(risk1 / pmax(risk0, 1e-9))

  list(risk1 = risk1, risk0 = risk0, rd = rd, logrr = logrr)
}

# Compute anchored theta^(delta) curve for a given baseline estimator theta0
# - theta0_* can be IPCW or DR, etc.
anchored_theta_curve <- function(pp_dt,
                                 t_cut,
                                 breaks,
                                 delta_grid,
                                 shape = c("constant","late_only","sign_change"),
                                 t0 = NULL,
                                 theta0_logrr,
                                 theta0_rd) {
  shape <- match.arg(shape, c("constant","late_only","sign_change"))
  # virtual movement baseline at delta=0
  d0 <- make_delta_profile(breaks, shape = shape, delta = 0.0, t0 = t0)
  v0 <- virtual_risk_under_delta(pp_dt, t_cut = t_cut, breaks = breaks, delta_vec = d0)

  out <- data.table::data.table(delta = delta_grid, logrr = NA_real_, rd = NA_real_,
                                v_logrr = NA_real_, v_rd = NA_real_)
  for (j in seq_along(delta_grid)) {
    dvec <- make_delta_profile(breaks, shape = shape, delta = delta_grid[j], t0 = t0)
    v <- virtual_risk_under_delta(pp_dt, t_cut = t_cut, breaks = breaks, delta_vec = dvec)
    # anchored movement
    out$v_logrr[j] <- v$logrr
    out$v_rd[j]    <- v$rd
    out$logrr[j]   <- theta0_logrr + (v$logrr - v0$logrr)
    out$rd[j]      <- theta0_rd    + (v$rd    - v0$rd)
  }
  out
}

# Delta-star: minimal delta where decision flips relative to delta=0
# decision_mode:
#   - "sign": sign flip of point estimate
#   - "cross_zero": estimate crosses 0 (more stable near 0)
compute_delta_star <- function(delta_grid, theta_curve, measure = c("logrr","rd"),
                               decision_mode = c("cross_zero","sign")) {
  measure <- match.arg(measure, c("logrr","rd"))
  decision_mode <- match.arg(decision_mode, c("cross_zero","sign"))

  y <- theta_curve[[measure]]
  if (all(!is.finite(y))) return(NA_real_)

  y0 <- y[which.min(abs(delta_grid - 0))]
  if (!is.finite(y0)) y0 <- y[1]

  # define baseline "side"
  if (decision_mode == "sign") {
    side0 <- sign(y0)
    if (side0 == 0) return(NA_real_)
    side <- sign(y)
    flip_idx <- which(side != side0 & side != 0)
  } else {
    # cross_zero: change in sign or exact crossing
    flip_idx <- which(y0 * y <= 0 & is.finite(y))
    # remove delta=0 point itself
    flip_idx <- setdiff(flip_idx, which.min(abs(delta_grid - 0)))
  }

  if (length(flip_idx) == 0) return(Inf)
  min(delta_grid[flip_idx], na.rm = TRUE)
}

# Optional: time-resolved delta-star curve by moving the "late-only start" along breaks
delta_star_curve_by_start <- function(pp_dt,
                                      t_cut,
                                      breaks,
                                      delta_grid,
                                      theta0_logrr,
                                      theta0_rd,
                                      decision_mode = c("cross_zero","sign")) {
  decision_mode <- match.arg(decision_mode, c("cross_zero","sign"))
  tstart <- breaks[-length(breaks)]
  K <- length(tstart)

  # virtual baseline at delta=0
  d0 <- rep(0.0, K)
  v0 <- virtual_risk_under_delta(pp_dt, t_cut = t_cut, breaks = breaks, delta_vec = d0)

  res <- data.table::data.table(start_interval = seq_len(K),
                                t_start = tstart,
                                delta_star_logrr = NA_real_,
                                delta_star_rd = NA_real_)

  for (k0 in seq_len(K)) {
    # build curves where delta applies from k0 onwards
    curve <- data.table::data.table(delta = delta_grid, logrr = NA_real_, rd = NA_real_)
    for (j in seq_along(delta_grid)) {
      dvec <- rep(0.0, K)
      dvec[k0:K] <- delta_grid[j]
      v <- virtual_risk_under_delta(pp_dt, t_cut = t_cut, breaks = breaks, delta_vec = dvec)
      curve$logrr[j] <- theta0_logrr + (v$logrr - v0$logrr)
      curve$rd[j]    <- theta0_rd    + (v$rd    - v0$rd)
    }
    res$delta_star_logrr[k0] <- compute_delta_star(delta_grid, curve, measure = "logrr", decision_mode = decision_mode)
    res$delta_star_rd[k0]    <- compute_delta_star(delta_grid, curve, measure = "rd", decision_mode = decision_mode)
  }
  res
}
