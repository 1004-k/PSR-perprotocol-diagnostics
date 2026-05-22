#!/usr/bin/env Rscript
# scripts/24_make_paperc_method_compare_figures.R
# ------------------------------------------------------------
# Compare IPCW vs DR(glm) vs DR(ML) vs TMLE on:
# - baseline performance summaries (ci_excl0, rmse) if present
# - delta* summaries if present
#
# This generates Figure 4 in the main manuscript.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- Sys.getenv("OUT_DIR", "")
if (!nzchar(out_dir)) stop("OUT_DIR is empty. Set OUT_DIR and rerun.")

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# delta* summary (already aggregated)
ds_sum_file <- file.path(out_dir, "perf_summary_paperc_delta_star.csv")
if (!file.exists(ds_sum_file)) stop("Missing: ", ds_sum_file)

ds <- fread(ds_sum_file)
# Focus: truth=null, delta_shape=constant, measure logRR
ds <- ds[truth == "null" & delta_shape == "constant"]

# ensure scenario ordering
ds[, scenario_id := factor(scenario_id, levels = sort(unique(as.character(scenario_id))))]

p1 <- ggplot(ds, aes(x = scenario_id, y = med_delta_star_logrr, shape = method)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  facet_grid(delta_dgm ~ mis_spec, labeller = label_both) +
  labs(x = "Scenario", y = "Median delta* (logRR)", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Figure4.pdf"), p1, width = 10, height = 6)

message("Saved method comparison figure(s) to: ", fig_dir)
