# scripts/26_run_mimic_demo_worked_example.R
# ------------------------------------------------------------
# Worked example: MIMIC-IV Clinical Database Demo (v2.2)
#
# v4 fixes two common Windows pitfalls:
#   (1) Working-directory drift (relative paths point to a different folder)
#   (2) exdir exists but is empty/partial, and previous scripts didn't re-unzip
#
# This script will:
#   - Re-unzip if hosp/admissions.csv.gz cannot be found under exdir
#   - Locate the true data root even if the zip contains an extra top folder
#   - Run the main SPD example, plus optional piecewise SPD(t)
# ------------------------------------------------------------

req_pkgs <- function(pkgs) {
  to_install <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(to_install)) install.packages(to_install)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# ---- Output directory ----
out_dir <- Sys.getenv("OUT_DIR", "output_paperc_mimic_demo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
  cfg <- init_project(n_cores = 1L, seed = 2026L, out_dir = out_dir)
  req_pkgs <- cfg$require_pkgs
}

req_pkgs(c("data.table", "survival"))


# ---- (1) Locate + ensure-extract MIMIC demo ----
zip_path <- Sys.getenv("MIMIC_ZIP", file.path("data", "mimic-iv-clinical-database-demo-2.2.zip"))
exdir    <- Sys.getenv("MIMIC_EXTRACT_DIR", file.path("data", "mimic_demo_2.2"))

if (!file.exists(zip_path)) {
  stop(
    "Cannot find MIMIC demo zip at: ", zip_path, "\n",
    "Fix: set an absolute path, e.g.\n",
    "  Sys.setenv(MIMIC_ZIP='C:/.../mimic-iv-clinical-database-demo-2.2.zip')\n",
    "Or set your working directory to the project root (where /data exists)."
  )
}

dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

find_admissions <- function(base_dir) {
  # robust cross-platform search (Windows uses \\ in returned paths)
  cand <- list.files(base_dir, recursive = TRUE, full.names = TRUE, pattern = "admissions\\.csv\\.gz$")
  if (!length(cand)) return(character(0))
  cand_norm <- gsub("\\\\", "/", cand)
  ok <- grepl("/hosp/admissions\\.csv\\.gz$", cand_norm)
  cand[ok]
}

adm_candidates <- find_admissions(exdir)
if (!length(adm_candidates)) {
  message("Could not find hosp/admissions.csv.gz under: ", normalizePath(exdir, winslash = "/", mustWork = FALSE))
  message("Re-unzipping the demo zip into exdir ...")
  utils::unzip(zipfile = zip_path, exdir = exdir)
  adm_candidates <- find_admissions(exdir)
}

if (!length(adm_candidates)) {
  stop(
    "Still could not locate hosp/admissions.csv.gz under: ", exdir, "\n",
    "Most common causes:\n",
    "  - Your working directory is not the project root (relative paths point elsewhere).\n",
    "  - The zip did not extract properly.\n\n",
    "Quick fixes:\n",
    "  1) setwd('C:/.../확장') to the project root, then re-run.\n",
    "  2) Delete the folder '", exdir, "' completely, then re-run.\n",
    "  3) Use absolute paths via Sys.setenv(MIMIC_ZIP=..., MIMIC_EXTRACT_DIR=...).\n"
  )
}

# Root is the folder that directly contains /hosp and /icu
root <- dirname(dirname(adm_candidates[1]))

paths <- list(
  icu_icustays = file.path(root, "icu",  "icustays.csv.gz"),
  hosp_adm     = file.path(root, "hosp", "admissions.csv.gz"),
  hosp_pat     = file.path(root, "hosp", "patients.csv.gz")
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) stop("Missing file: ", paths[[nm]])
}

# ---- (2) Read tables ----
icu <- data.table::fread(paths$icu_icustays)
adm <- data.table::fread(paths$hosp_adm)
pat <- data.table::fread(paths$hosp_pat)

# ---- (3) Build ICU cohort ----
tau_hours <- as.numeric(Sys.getenv("TAU_HOURS", "72"))

icu[, intime  := as.POSIXct(intime,  format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
icu[, outtime := as.POSIXct(outtime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
icu[, los_hours := as.numeric(difftime(outtime, intime, units = "hours"))]

setkey(icu, subject_id, hadm_id)
setkey(adm, subject_id, hadm_id)
setkey(pat, subject_id)

cohort <- icu[!is.na(los_hours) & los_hours > 0]
cohort <- adm[cohort]
cohort <- pat[cohort, on = .(subject_id)]

cohort[, female := as.integer(gender == "F")]
cohort[, age    := as.numeric(anchor_age)]
cohort[, death_in_hosp := as.integer(hospital_expire_flag == 1)]
cohort <- cohort[!is.na(age) & !is.na(female) & !is.na(death_in_hosp)]

# ---- (4) Prognostic model (tiny + stable) ----
prog_fit <- glm(death_in_hosp ~ scale(age) + female,
                data = cohort, family = binomial())

cohort[, s_hat := as.numeric(predict(prog_fit, type = "link"))]
cohort[, z     := as.numeric(scale(s_hat))]

# ---- (5) Deviation definition ----
cohort[, time := pmin(los_hours, tau_hours)]
cohort[, dev  := as.integer(los_hours <= tau_hours)]
cohort <- cohort[time > 0 & !is.na(dev) & !is.na(z)]

# ---- (6) Fit SPD (time-constant) ----
spd_fit <- coxph(Surv(time, dev) ~ z, data = cohort, ties = "efron")

gamma_std <- unname(coef(spd_fit)[["z"]])
hr_1sd    <- exp(gamma_std)
ci        <- exp(confint(spd_fit)["z", ])

res <- data.table(
  n_stays      = nrow(cohort),
  n_deviations = sum(cohort$dev == 1),
  tau_hours    = tau_hours,
  gamma_std    = gamma_std,
  HR_per_1SD   = hr_1sd,
  CI_low       = ci[[1]],
  CI_high      = ci[[2]]
)

print(res)

out_csv <- file.path(out_dir, "worked_example_mimic_demo_spd.csv")
fwrite(res, out_csv)

out_txt <- file.path(out_dir, "worked_example_mimic_demo_spd.txt")
cat(
  "MIMIC-IV Demo (v2.2) worked example (illustration only)\n",
  sprintf("Cohort: ICU stays (N=%d), deviation=ICU discharge within %.0f hours (events=%d)\n",
          res$n_stays, tau_hours, res$n_deviations),
  "Prognostic score: log-odds of in-hospital death from logistic model (age, sex), standardized (z)\n",
  sprintf("SPD (time-constant): gamma_std=%.3f; HR per 1 SD higher prognosis = %.3f (95%% CI %.3f to %.3f)\n",
          res$gamma_std, res$HR_per_1SD, res$CI_low, res$CI_high),
  file = out_txt
)

message("Saved: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
message("Saved: ", normalizePath(out_txt, winslash = "/", mustWork = FALSE))

# ---- (7) Optional: piecewise SPD(t) ----
do_piecewise <- as.integer(Sys.getenv("DO_PIECEWISE", "1"))
if (do_piecewise == 1L) {
  tryCatch({
    cut_points <- c(24, 48, tau_hours)
    cut_points <- cut_points[cut_points < tau_hours]

    long <- survival::survSplit(
      data = cohort,
      cut  = cut_points,
      end  = "time",
      start= "tstart",
      event= "dev",
      episode = "interval"
    )
    long <- data.table::as.data.table(long)
    long[, interval_f := factor(interval)]

    spdt_fit <- coxph(Surv(tstart, time, dev) ~ 0 + interval_f:z + strata(interval_f),
                      data = long, ties = "efron")

    coefs <- coef(spdt_fit)
    hr <- exp(coefs)
    ci_mat <- exp(confint(spdt_fit))

    out_piece <- data.table(
      interval  = gsub("interval_f", "", sub(":z$", "", names(coefs))),
      t_start_h = c(0, cut_points)[seq_along(coefs)],
      t_stop_h  = c(cut_points, tau_hours)[seq_along(coefs)],
      gamma_std = as.numeric(coefs),
      HR_per_1SD = as.numeric(hr),
      CI_low   = as.numeric(ci_mat[, 1]),
      CI_high  = as.numeric(ci_mat[, 2])
    )


    # ---- eTable 6 (manuscript-ready formatting) ----
    etab3 <- out_piece[, .(
      `Interval (hours)` = sprintf("%.0f-%.0f", t_start_h, t_stop_h),
      gamma_std = sprintf("%.3f", gamma_std),
      `HR per 1 SD` = sprintf("%.3f", HR_per_1SD),
      `95% CI` = sprintf("%.3f to %.3f", CI_low, CI_high)
    )]
    etab3_path <- file.path(out_dir, "tables", "eTable6_MIMIC_demo_piecewise_SPDt.csv")
    dir.create(dirname(etab3_path), showWarnings = FALSE, recursive = TRUE)
    fwrite(etab3, etab3_path)
    message("Saved: ", normalizePath(etab3_path, winslash = "/", mustWork = FALSE))

    # ---- Notes text for manuscript ----
    notes_path <- file.path(out_dir, "tables", "eTable6_MIMIC_demo_piecewise_SPDt_notes.txt")
    cat(
      "eTable 6. MIMIC-IV demo implementation vignette: piecewise SPD(t) for ICU discharge within 72 hours\n",
      "Notes: Cohort is ICU stays from the MIMIC-IV Clinical Database Demo (v2.2; N=", res$n_stays, "). ",
      "Protocol deviation was defined as ICU discharge within tau = ", tau_hours, " hours (events = ", res$n_deviations, "); ",
      "follow-up was administratively censored at ", tau_hours, " hours. ",
      "The prognostic score was the standardized (z) linear predictor (log-odds) from a logistic model for in-hospital death using age and sex. ",
      "SPD(t) was estimated within each time interval using a Cox model for the deviation hazard; HR per 1 SD corresponds to exp(gamma_std). ",
      "This worked example is provided to illustrate implementation and reporting only (no substantive clinical inference).\n",
      sep = "",
      file = notes_path
    )
    message("Saved: ", normalizePath(notes_path, winslash = "/", mustWork = FALSE))


    out_piece_csv <- file.path(out_dir, "worked_example_mimic_demo_spdt_piecewise.csv")
    fwrite(out_piece, out_piece_csv)
    message("Saved: ", normalizePath(out_piece_csv, winslash = "/", mustWork = FALSE))
  }, error = function(e) {
    message("Piecewise SPD(t) step skipped due to error: ", conditionMessage(e))
  })
}

message("Done.")