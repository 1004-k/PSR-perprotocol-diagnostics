#!/usr/bin/env Rscript
# scripts/23_make_paperc_figures.R
# ------------------------------------------------------------
# Paper C figures (black/white friendly):
# - Figure 1: Operating map (SPD summary x min rESS) with symbols showing quartiles of median finite delta*
# - Figure 2: Bundle time-resolved example panels (SPD(t), rESS(t), tail-share(t)) for a representative scenario
# - Figure 3: Fragility curve theta^(delta) for a representative scenario (requires SAVE_SENS_CURVES=1)
# ------------------------------------------------------------

source("R/00_utils.R")
cfg <- init_project(n_cores = 1L, seed = 2026L, out_dir = Sys.getenv("OUT_DIR", "output_paperc"))
cfg$require_pkgs(c("data.table"))

source("R/01_scenarios.R")
source("R/07_plotting.R")

raw_dir <- file.path(cfg$out_dir, "raw")
rep_file <- file.path(raw_dir, "replicate_results_paperc.csv")
ds_file  <- file.path(raw_dir, "delta_star_paperc.csv")
spd_file <- file.path(raw_dir, "spd_curves_paperc.csv")
wd_file  <- file.path(raw_dir, "weight_diagnostics_paperc.csv")
sc_file  <- file.path(raw_dir, "sensitivity_curves_paperc.csv")

if (!file.exists(rep_file) || !file.exists(ds_file)) {
  stop("Missing required raw outputs. Run scripts/21_run_paperc_simulations.R first.")
}

rep <- data.table::fread(rep_file, showProgress = FALSE)
ds  <- data.table::fread(ds_file, showProgress = FALSE)
scn <- make_scenario_grid()[, .(scenario_id, scenario_num, panel_title)]
rep <- merge(rep, scn, by = "scenario_id", all.x = TRUE)

fig_dir <- file.path(cfg$out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_finite <- function(x) any(is.finite(x))

# Pick a valid sensitivity curve for Figure 3. Full manuscript runs should use
# the requested representative scenario. Small quickcheck runs can occasionally
# produce NA-only fragility curves in the requested replicate because a Cox model
# inside the virtual sensitivity engine does not converge. The fallback keeps the
# quickcheck useful by selecting the first available finite curve.
pick_sens_curve <- function(sc, sid, ms, dd, tr, rp, meth, shp) {
  sc1 <- sc[scenario_id == sid & mis_spec == ms & delta_dgm == dd & truth == tr &
            replicate == rp & method == meth & delta_shape == shp]
  if (nrow(sc1) > 0 && "logrr" %in% names(sc1) && has_finite(sc1$logrr)) return(sc1)

  cand <- sc[method == meth & delta_shape == shp & is.finite(logrr)]
  if (nrow(cand) == 0) cand <- sc[delta_shape == shp & is.finite(logrr)]
  if (nrow(cand) == 0) cand <- sc[is.finite(logrr)]
  if (nrow(cand) == 0) return(sc1[0])

  key <- cand[, .N, by = .(scenario_id, mis_spec, delta_dgm, truth, replicate, method, delta_shape)][1]
  sc[scenario_id == key$scenario_id & mis_spec == key$mis_spec &
     delta_dgm == key$delta_dgm & truth == key$truth & replicate == key$replicate &
     method == key$method & delta_shape == key$delta_shape]
}

# ----------------------------
# Figure 1: Operating map
# ----------------------------
# Use median delta* (logRR) by scenario and mis_spec for IPCW, constant delta-shape, truth=null.
map_truth <- Sys.getenv("MAP_TRUTH", "null")
ds1 <- ds[method == "IPCW" & delta_shape == "constant" & truth == map_truth]

# Summaries for operating map:
# - med_delta_star_finite: median delta* among finite values (if any)
# - p_flip: proportion of replicates with a finite delta* (i.e., conclusion flips within the explored delta grid)
ds_sum <- ds1[, .(
  n_total = .N,
  n_finite = sum(is.finite(delta_star_logrr)),
  p_flip = if (.N > 0) sum(is.finite(delta_star_logrr)) / .N else NA_real_,
  med_delta_star_finite = if (sum(is.finite(delta_star_logrr)) > 0)
    stats::median(delta_star_logrr[is.finite(delta_star_logrr)], na.rm = TRUE)
  else NA_real_
), by = .(scenario_id, mis_spec, delta_dgm)]

op <- rep[truth == map_truth, .(
  max_spd_med = stats::median(max_spd, na.rm = TRUE),
  min_rESS_med = stats::median(min_rESS, na.rm = TRUE)
), by = .(scenario_id, mis_spec, delta_dgm)]

op <- merge(op, ds_sum, by = c("scenario_id","mis_spec","delta_dgm"), all.x = TRUE)

# Quartiles of median finite delta* for visual coding. If no finite crossing is
# observed for a point, it is plotted as an open circle. In small quickcheck
# runs it is possible that all delta* values are infinite. In that case we still
# draw the operating map and label all points as "No finite crossing" rather
# than failing at the quantile step.
finite_vals <- op$med_delta_star_finite[is.finite(op$med_delta_star_finite)]
if (length(finite_vals) > 0) {
  finite_q <- stats::quantile(finite_vals, probs = c(0.25, 0.50, 0.75), na.rm = TRUE, type = 7)
  q1 <- finite_q[[1]]; q2 <- finite_q[[2]]; q3 <- finite_q[[3]]
  op[, qbin := data.table::fifelse(!is.finite(med_delta_star_finite), "No finite crossing",
                data.table::fifelse(med_delta_star_finite <= q1, "Q1",
                data.table::fifelse(med_delta_star_finite <= q2, "Q2",
                data.table::fifelse(med_delta_star_finite <= q3, "Q3", "Q4"))))]
  legend_labs <- c(
    Q1 = paste0("Q1: delta*<=", sprintf("%.2f", q1)),
    Q2 = paste0("Q2: delta*=", sprintf("%.2f", q1), "-", sprintf("%.2f", q2)),
    Q3 = paste0("Q3: delta*=", sprintf("%.2f", q2), "-", sprintf("%.2f", q3)),
    Q4 = paste0("Q4: delta*>", sprintf("%.2f", q3)),
    "No finite crossing" = "No finite crossing"
  )
} else {
  q1 <- q2 <- q3 <- NA_real_
  op[, qbin := "No finite crossing"]
  legend_labs <- c("No finite crossing" = "No finite crossing")
}
op[, qbin := factor(qbin, levels = c("Q1", "Q2", "Q3", "Q4", "No finite crossing"))]

pch_map <- c(Q1 = 21, Q2 = 22, Q3 = 24, Q4 = 25, "No finite crossing" = 1)
bg_map  <- c(Q1 = "grey80", Q2 = "grey60", Q3 = "grey45", Q4 = "grey25", "No finite crossing" = NA)

fig1 <- file.path(fig_dir, "Figure1.pdf")
mis_levels <- sort(unique(op$mis_spec))
n_panels <- length(mis_levels)
fig1_width <- if (n_panels >= 2) 10.5 else 6.2
open_pdf(fig1, width = fig1_width, height = 5.6, pointsize = 12)
par(mfrow = c(1, n_panels), mar = c(4.2,4.4,2.4,1.2), oma = c(0.6,0.6,1.2,0.4))

present <- levels(op$qbin)[levels(op$qbin) %in% as.character(unique(op$qbin))]
only_no_cross <- identical(present, "No finite crossing")

for (ms in mis_levels) {
  dms <- op[mis_spec == ms]

  plot(dms$max_spd_med, dms$min_rESS_med,
       type = "n",
       xlab = "Median max |SPD(t)|",
       ylab = "Median min rESS(t)",
       main = paste0("mis_spec=", ms))

  for (b in levels(op$qbin)) {
    dd <- dms[qbin == b]
    if (nrow(dd) == 0) next
    points(dd$max_spd_med, dd$min_rESS_med,
           pch = pch_map[[as.character(b)]],
           bg = bg_map[[as.character(b)]],
           col = "grey20", cex = 1.2)
  }

  abline(h = 0.25, lty = 3, col = "grey60")
  box()

  # In reduced quickcheck runs, all delta* values can be infinite because the
  # run uses few replicates and a narrow delta grid. In that case, avoid a large
  # legend over the plotting region and label the condition unobtrusively.
  if (only_no_cross) {
    mtext("No finite crossing within the explored grid", side = 3, line = -1.0,
          adj = 0.02, cex = 0.75, col = "grey35")
  } else if (ms == max(mis_levels)) {
    legend("topright",
           legend = unname(legend_labs[present]),
           pch = unname(pch_map[present]),
           pt.bg = unname(bg_map[present]),
           pt.cex = 1.0,
           bty = "n",
           cex = 0.75,
           title = "Median finite delta* quartiles")
  }
}

# Figure title is provided in the manuscript legend, not in the graphic file.
close_device()
message("Saved: ", fig1)

# ----------------------------
# Figure 2: Bundle time-resolved panels (one scenario)
# ----------------------------
if (file.exists(spd_file) && file.exists(wd_file)) {
  spd <- data.table::fread(spd_file, showProgress = FALSE)
  wd  <- data.table::fread(wd_file, showProgress = FALSE)

  # pick a representative scenario: S09, mis_spec=0, delta_dgm=D2_constant, truth=null, replicate=1
  sid <- Sys.getenv("EXAMPLE_SCENARIO", "S09")
  ms  <- as.integer(Sys.getenv("EXAMPLE_MIS_SPEC", "0"))
  dd  <- Sys.getenv("EXAMPLE_DELTA_DGM", "D2_constant")
  tr  <- Sys.getenv("EXAMPLE_TRUTH", "null")
  rp  <- as.integer(Sys.getenv("EXAMPLE_REP", "1"))

  spd1 <- spd[scenario_id == sid & mis_spec == ms & delta_dgm == dd & truth == tr & replicate == rp]
  wd1  <- wd[scenario_id == sid & mis_spec == ms & delta_dgm == dd & truth == tr & replicate == rp]

  if (nrow(spd1) > 0 && nrow(wd1) > 0) {
    fig2 <- file.path(fig_dir, "Figure2.pdf")
    open_pdf(fig2, width = 10.8, height = 5.0, pointsize = 12)
    par(mfrow = c(1,3), mar = c(4.4,4.6,3.6,1.0), oma = c(0.8,0.6,2.0,0.4))

    plot(spd1$t_mid, spd1$gamma_hat, type = "l", lwd = 3,
         xlab = "Time", ylab = "SPD(t) (log-HR per 1 SD)",
         main = "SPD(t)")
    abline(h = 0, lty = 3, col = "grey60")

    plot(wd1$t, wd1$rESS, type = "l", lwd = 3,
         xlab = "Time", ylab = "rESS(t)",
         main = "rESS(t)")
    abline(h = 0.25, lty = 3, col = "grey60")

    plot(wd1$t, wd1$tail_share, type = "l", lwd = 3,
         xlab = "Time", ylab = "Tail-share(t) (top 1%)",
         main = "Tail-share(t)")
    abline(h = 0.10, lty = 3, col = "grey60")

    # Figure title is provided in the manuscript legend, not in the graphic file.
    close_device()
    message("Saved: ", fig2)
  }
}

# ----------------------------
# Figure 3: Fragility curve (requires sensitivity_curves_paperc.csv)
# ----------------------------
if (file.exists(sc_file)) {
  sc <- data.table::fread(sc_file, showProgress = FALSE)

  sid <- Sys.getenv("FRAG_SCENARIO", "S09")
  ms  <- as.integer(Sys.getenv("FRAG_MIS_SPEC", "0"))
  dd  <- Sys.getenv("FRAG_DELTA_DGM", "D2_constant")
  tr  <- Sys.getenv("FRAG_TRUTH", "null")
  rp  <- as.integer(Sys.getenv("FRAG_REP", "1"))
  meth <- Sys.getenv("FRAG_METHOD", "IPCW")
  shp  <- Sys.getenv("FRAG_SHAPE", "constant")

  sc1 <- pick_sens_curve(sc, sid = sid, ms = ms, dd = dd, tr = tr, rp = rp,
                         meth = meth, shp = shp)

  if (nrow(sc1) > 0 && "delta" %in% names(sc1) && "logrr" %in% names(sc1) && has_finite(sc1$logrr)) {
    sc1 <- sc1[is.finite(delta) & is.finite(logrr)]
    data.table::setorder(sc1, delta)
    fig3 <- file.path(fig_dir, "Figure3.pdf")
    open_pdf(fig3, width = 7.6, height = 4.4, pointsize = 12)
    par(mfrow = c(1,1), mar = c(4.2,4.4,2.4,1.0), oma = c(0.4,0.4,0.6,0.2))

    plot(sc1$delta, sc1$logrr, type = "l", lwd = 3,
         xlab = expression(delta),
         ylab = "Anchored logRR(delta)",
         main = "")
    abline(h = 0, lty = 3, col = "grey60")

    close_device()
    message("Saved: ", fig3)
  } else {
    message("Skipped Figure3: no finite sensitivity curve was available in this quickcheck run.")
  }
}

cat("Done. Figures saved to: ", fig_dir, "\n", sep = "")
