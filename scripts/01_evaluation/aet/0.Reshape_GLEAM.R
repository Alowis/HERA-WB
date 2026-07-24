###############################################################################
# AET DIEGO-STYLE - STEP 0: RESHAPE GLEAM DAILY EXTRACTION TO WIDE FORMAT
#
# The extract_gleam_aet_daily.R script produces a LONG CSV:
#   date, Id, n_pixels, n_valid, n_na, mean_w, sd, p2.5, p50, p97.5
# (one row per catchment per day, appended year by year).
#
# This step pivots it to WIDE format matching the LISFLOOD structure:
#   date, <catch_id_1>, <catch_id_2>, ...
# using only the `mean_w` column (area-weighted mean AET).
# The wide CSV is what the Homogenize step consumes.
###############################################################################

library(data.table)

base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
in_csv <- file.path(base_dir, "data", "gleam_aet_catchment_daily.csv")
out_csv <- file.path(base_dir, "data", "gleam_aet_daily_wide.csv")

message("Reading long extraction CSV...")
dt <- data.table::fread(in_csv, select = c("date", "Id", "mean_w"))
dt[, Id := as.character(as.integer(Id))] # normalise: strip leading zeros

message("Pivoting to wide (date x catchments)...")
wide <- data.table::dcast(dt, date ~ Id, value.var = "mean_w")
data.table::setorder(wide, date)

message("Writing wide CSV: ", out_csv)
data.table::fwrite(wide, out_csv)
message(
    "  rows (dates): ", nrow(wide),
    " | cols (catchments): ", ncol(wide) - 1
)
