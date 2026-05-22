# RStudio wrapper: regenerate Paper C tables from an existing OUT_DIR
# ------------------------------------------------------------
# 1) Set OUT_DIR below to the folder that contains raw/ and logs/
#    e.g., "output_paperc_main" or an absolute path.
# 2) Run this script in RStudio (Source).
# ------------------------------------------------------------

OUT_DIR <- "output_paperc_main"   # <-- change if needed
stopifnot(dir.exists(OUT_DIR))

Sys.setenv(OUT_DIR = OUT_DIR)

# This script will create/overwrite OUT_DIR/tables/
source("scripts/22_make_paperc_tables.R")

message("Done. Tables saved under: ", file.path(OUT_DIR, "tables"))
