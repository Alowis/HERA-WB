# =============================================================================
# Preprocess TSS data into daily and monthly aggregate CSVs
# Saves output in the same format as existing aggregates:
#   - {var}_daily_all_years.csv
#   - {var}_monthly_all_years.csv
# =============================================================================

library(data.table)
library(lubridate)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
agg_dir <- file.path(base_dir, "data", "aggregates")

# --- All available variables ---------------------------------------------------
# name: folder name for output
# file: TSS filename
# agg_method: "sum" for fluxes, "mean" for states
all_variables <- list(
    list(name = "rainfall", file = "rainUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "snowfall", file = "snowUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "snowmelt", file = "snowMeltUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "infiltration", file = "infUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "ActEvapo", file = "ActEvapo_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "runoff", file = "surfaceRunoffUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "prefflow", file = "prefFlowUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "dSubToUz", file = "dSubToUzUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "percolation", file = "percUZLZUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "quz", file = "qUzUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "qlz", file = "qLZUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "gwloss", file = "lossUpsX_nested_1951_2020.csv", agg_method = "sum"),
    list(name = "theta1", file = "theta1totalX_nested_1951_2020.csv", agg_method = "mean"),
    list(name = "theta2", file = "theta2totalX_nested_1951_2020.csv", agg_method = "mean"),
    list(name = "theta3", file = "theta3totalX_nested_1951_2020.csv", agg_method = "mean"),
    list(name = "snow_water_equivalent", file = "scovUps_nested_1951_2020.csv", agg_method = "mean"),
    list(name = "discharge", file = "disWin_nested_1951_2020.csv", agg_method = "mean")
)

# --- Choose which variables to process ----------------------------------------
# Set to NULL to process ALL, or provide a character vector of names to process
# Example: process_only <- c("rainfall", "snowfall", "theta3")
process_only <- c("prefflow") # NULL = process all

if (!is.null(process_only)) {
    variables <- all_variables[sapply(all_variables, function(v) v$name %in% process_only)]
    cat("Processing selected variables:", paste(process_only, collapse = ", "), "\n")
} else {
    variables <- all_variables
    cat("Processing ALL", length(variables), "variables\n")
}

# --- Date sequence ------------------------------------------------------------
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1
date_daily <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")

# --- Process each variable ----------------------------------------------------
for (v in variables) {
    var_name <- v$name
    var_file <- v$file
    agg_method <- v$agg_method

    out_folder <- file.path(agg_dir, var_name)
    dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

    daily_out <- file.path(out_folder, paste0(var_name, "_daily_all_years.csv"))
    monthly_out <- file.path(out_folder, paste0(var_name, "_monthly_all_years.csv"))

    # Skip if monthly already exists
    if (file.exists(monthly_out)) {
        cat("[", var_name, "] Monthly aggregate already exists, skipping.\n")
        next
    }

    cat("[", var_name, "] Reading", var_file, "...\n")
    dt <- fread(file.path(tss_dir, var_file), header = TRUE)

    n_rows <- nrow(dt)
    catch_cols <- names(dt)

    # Assign dates
    if (n_rows == n_days) {
        dt[, date := date_daily]
        cat("  -> Daily data (", n_rows, " rows)\n")
    } else if (n_rows == n_days * 4) {
        cat("  -> 6-hourly data, aggregating to daily...\n")
        dt[, day_idx := rep(seq_len(n_days), each = 4)]
        dt[, date := date_daily[day_idx]]
        catch_cols <- setdiff(names(dt), c("day_idx", "date"))
        if (agg_method == "sum") {
            dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = date, .SDcols = catch_cols]
        } else {
            dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = date, .SDcols = catch_cols]
        }
        catch_cols <- setdiff(names(dt), "date")
    } else {
        warning("  Unexpected row count (", n_rows, ") for ", var_name, ". Skipping.")
        next
    }

    # Save daily aggregate
    cat("  Saving daily aggregate...\n")
    daily_save <- data.table(
        period_start = dt$date,
        period_end = dt$date
    )
    daily_save <- cbind(daily_save, dt[, ..catch_cols])
    fwrite(daily_save, daily_out)
    cat("  -> Saved:", daily_out, "\n")

    # Compute monthly aggregates
    cat("  Computing monthly aggregates...\n")
    dt[, year := year(date)]
    dt[, month := month(date)]

    if (agg_method == "sum") {
        monthly <- dt[, lapply(.SD, sum, na.rm = TRUE),
            by = .(year, month), .SDcols = catch_cols
        ]
    } else {
        monthly <- dt[, lapply(.SD, mean, na.rm = TRUE),
            by = .(year, month), .SDcols = catch_cols
        ]
    }

    # Build period_start / period_end columns
    monthly[, period_start := as.Date(paste(year, month, "01", sep = "-"))]
    monthly[, period_end := period_start + days_in_month(period_start) - 1L]
    monthly[, month_idx := seq_len(.N)]

    # Reorder columns to match existing format
    meta <- c("month_idx", "period_start", "period_end")
    setcolorder(monthly, c(meta, catch_cols))
    monthly[, c("year", "month") := NULL]

    # Save
    fwrite(monthly, monthly_out)
    cat("  -> Saved:", monthly_out, "\n")

    # Free memory
    rm(dt, monthly, daily_save)
    gc()
}

cat("\nDone! Aggregates saved to:", agg_dir, "\n")
