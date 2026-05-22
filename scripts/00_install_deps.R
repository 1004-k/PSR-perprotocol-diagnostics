#!/usr/bin/env Rscript
# Install minimal dependencies for Paper C standalone repo.
# Usage:
#   Rscript scripts/00_install_deps.R
pkgs <- c(
  "data.table","survival","future","future.apply",
  "ggplot2","glmnet"
)
# Optional (only needed if you enable TMLE):
opt <- c("tmle","SuperLearner")

inst <- rownames(installed.packages())
need <- setdiff(pkgs, inst)
if (length(need) > 0) {
  install.packages(need, repos = "https://cloud.r-project.org")
}
# don't auto-install optional packages (TMLE can be heavy)
message("Installed core pkgs. Optional (TMLE): ", paste(setdiff(opt, inst), collapse = ", "))
