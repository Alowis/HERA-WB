
####################################################################################
library(data.table)
library(purrr)
library(stringr)
library(dplyr)
library(tidyverse)
library(lubridate)

# 1. SET PATHS
path_out  <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\1.Diego_Merged\\"
file_lisf <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\0.Alois\\surface_soil_moisture_daily_all_years.csv"
esa_files <- list.files("Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\0.Diego_SM_Extacted\\daily\\", 
                        pattern = "\\.csv$", full.names = TRUE)

# 2. PROCESS AND RESTRUCTURE ESA CCI FILES (Fast multi-file bind & pivot)
message("Processing ESA CCI files...")

esa_raw <- rbindlist(lapply(esa_files, function(f) {
  # Extract numeric ID using regex and format it to match Lisflood columns (e.g., "X7402")
  id_num <- str_extract(basename(f), "(?<=SM_)[0-9]+(?=_daily)")
  id_clean <- paste0("X", as.integer(id_num))
  
  # Read only the necessary columns to save memory and time
  df <- fread(f, select = c("date", "mean_w"))
  df[, catchment_id := id_clean]
  return(df)
}))

# Pivot ESA data into the optimal Wide Format (Rows = Dates, Columns = Catchments)
esa_wide <- dcast(esa_raw, date ~ catchment_id, value.var = "mean_w")

# 3. LOAD AND ALIGN LISFLOOD
message("Processing Lisflood data...")
lisflood <- fread(file_lisf)

# Find common temporal coverage bounds based on available ESA data
common_dates <- intersect(lisflood$date, esa_wide$date)

# Find catchments present in both datasets
lisf_cols_with_X <- paste0("X", names(lisflood))
common_cols_with_X <- intersect(lisf_cols_with_X, names(esa_wide))
common_cols_lisflood <- gsub("^X", "", common_cols_with_X)

# 4. FILTER AND ALIGN BOTH DATASETS EXACTLY
# For Lisflood: use the original numeric names + 'date'
lisf_final <- lisflood[date %in% common_dates, c("date", common_cols_lisflood), with = FALSE]

# For ESA CCI: use the 'X' prefixed names + 'date'
esa_final  <- esa_wide[date %in% common_dates, c("date", common_cols_with_X), with = FALSE]

# Rename Lisflood columns to include 'X' so both matrices match perfectly
names(lisf_final) <- paste0("X", names(lisf_final))
names(lisf_final)[1] <- "date" # Keep date column named "date"

# Ensure strict chronological order
setorder(lisf_final, date)
setorder(esa_final, date)

# 5. SAVE HOMOGENIZED BIG DATAFRAMES
# Round Lisflood data to 3 decimal places
num_cols_lisf <- setdiff(names(lisf_final), "date")
for (col in num_cols_lisf) {
  set(lisf_final, j = col, value = round(lisf_final[[col]], 3))
}

# Round ESA CCI data to 3 decimal places
num_cols_esa <- setdiff(names(esa_final), "date")
for (col in num_cols_esa) {
  set(esa_final, j = col, value = round(esa_final[[col]], 3))
}

# Save the newly compressed big dataframes
fwrite(lisf_final, paste0(path_out, "lisflood_homogenized.csv"))
fwrite(esa_final,  paste0(path_out, "esacci_homogenized.csv"))


#

