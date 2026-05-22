# RStudio wrapper: run MIMIC-IV demo implementation vignette (eTable 6 style)
# ------------------------------------------------------------
# Requirements:
# - You must have access to the MIMIC-IV Clinical Database Demo v2.2 ZIP
#   (or an extracted directory that contains hosp/ and icu/).
# ------------------------------------------------------------

OUT_DIR <- "output_paperc_mimic_demo"  # <-- change if needed

# Option A: point to the ZIP file
MIMIC_ZIP <- "~/mimic-iv-clinical-database-demo-2.2.zip"  # <-- change

# Option B: if you already extracted it, set MIMIC_DIR instead and leave MIMIC_ZIP as ""
MIMIC_DIR <- ""  # e.g., "~/mimic-iv-clinical-database-demo-2.2"

TAU_HOURS <- "72"
DO_PIECEWISE <- "1"

if (nzchar(MIMIC_DIR)) {
  stopifnot(dir.exists(path.expand(MIMIC_DIR)))
  Sys.setenv(MIMIC_DIR = path.expand(MIMIC_DIR))
} else {
  stopifnot(file.exists(path.expand(MIMIC_ZIP)))
  Sys.setenv(MIMIC_ZIP = path.expand(MIMIC_ZIP))
}

Sys.setenv(
  OUT_DIR = OUT_DIR,
  TAU_HOURS = TAU_HOURS,
  DO_PIECEWISE = DO_PIECEWISE
)

source("scripts/26_run_mimic_demo_worked_example.R")

message("Done. MIMIC demo worked example saved under: ", OUT_DIR)
