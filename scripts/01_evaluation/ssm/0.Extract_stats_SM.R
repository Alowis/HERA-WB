library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(lubridate)
library(stringr)
library(future)
library(future.apply)
library(progressr)
library(purrr)
library(readr)


out <- "D:\\gomezdi\\My Data\\0. Projects\\Articles\\2026\\Alois_LisFlood_DataPaper\\out\\"

# ------------------------------------------------------------
# 1. LOAD DATA (same as before)
# ------------------------------------------------------------

my_files <- list.files(
  "D:\\gomezdi\\My Data\\0. Projects\\DataHub\\ESA_CCI_SM_v09\\",
  pattern = ".nc$", full.names = TRUE
)

v <- st_read("D:\\gomezdi\\My Data\\0. Projects\\Articles\\2026\\Alois_LisFlood_DataPaper\\catchments_analysis_final_v3.gpkg")
st_geometry(v) <- "geom"
v$Id <- as.character(v$catch_id)
v <- st_make_valid(v)
ref_crs <- terra::crs(terra::rast(my_files[1]))
v <- st_transform(v, ref_crs)

# Pre-convert to SpatVector ONCE (faster inside workers)
v_vect <- terra::vect(v)

# ------------------------------------------------------------
# 2. FUNCTIONS
# ------------------------------------------------------------
read_sm <- function(file, thr_unc = 0.04) {
  r   <- terra::rast(file)
  sm  <- r[["sm"]]
  unc <- r[["sm_uncertainty"]]
  sm[unc > thr_unc] <- NA
  sm
}

# Process ONE file across ALL catchments → returns a data.frame
process_one_file <- function(file, v_sf) {
  sm   <- read_sm(file)
  date <- lubridate::ymd(stringr::str_extract(basename(file), "\\d{8}"))
  
  stats <- exactextractr::exact_extract(
    x    = sm,
    y    = v_sf,
    fun = function(df) {
      ok <- !is.na(df$value)
      if (sum(ok) == 0) {
        return(data.frame(
          n_pixels = nrow(df), n_valid = 0L, n_na = nrow(df),
          mean_w = NA_real_, sd = NA_real_,
          p2.5 = NA_real_, p50 = NA_real_, p97.5 = NA_real_
        ))
      }
      v <- df$value[ok]
      w <- df$coverage_fraction[ok]
      w <- w / sum(w)
      mw <- weighted.mean(v, w)
      
      data.frame(
        n_pixels = nrow(df), n_valid  = sum(ok),  n_na = sum(!ok),
        mean_w   = round(mw, 4),
        sd       = round(sqrt(weighted.mean((v - mw)^2, w)), 4),
        p2.5     = round(quantile(v, 0.025, na.rm = TRUE), 4),
        p50      = round(quantile(v, 0.500, na.rm = TRUE), 4),
        p97.5    = round(quantile(v, 0.975, na.rm = TRUE), 4)
      )
    },
    summarize_df = TRUE,
    progress     = FALSE
  )
  
  stats$date <- date
  stats$Id   <- v_sf$Id
  stats
}

# ------------------------------------------------------------
# 3. PARALLEL OVER FILES  (each worker does all catchments)
# ------------------------------------------------------------
run_parallel_by_file <- function(v_sf, files, ncores = 4) {
  plan(multisession, workers = ncores)
  
  handlers(handler_progress(
    format = "[:bar] :percent | ETA: :eta | :message",
    clear  = FALSE
  ))
  
  with_progress({
    p <- progressr::progressor(steps = length(files))
    
    daily_list <- future.apply::future_lapply(
      files,
      function(f) {
        result <- tryCatch(
          process_one_file(f, v_sf),
          error = function(e) {
            message("ERROR in ", basename(f), ": ", e$message)
            NULL
          }
        )
        p(message = basename(f))
        result
      },
      future.seed = TRUE
    )
  })
  
  plan(sequential)
  
  daily_list <- purrr::compact(daily_list)
  
  if (length(daily_list) == 0) {
    stop("All file extractions failed. Check input files or polygon overlap.")
  }
  
  # Combine, sort, strip row names, and split into a list by polygon Id
  daily_df <- dplyr::bind_rows(daily_list) |>
    dplyr::arrange(Id, date)
  
  rownames(daily_df) <- NULL
  
  daily_by_id <- split(daily_df, daily_df$Id)
  
  # Remove row names from individual split data frames as well
  daily_by_id <- lapply(daily_by_id, function(df) {
    rownames(df) <- NULL
    df
  })
  
  daily_by_id
}

# ------------------------------------------------------------
# 4. RUN and save
# ------------------------------------------------------------

daily <- run_parallel_by_file(v, my_files, ncores = 4)

for (id in names(daily)) {
  file_name <- paste0(out,"daily\\SM_", id, "_daily.csv")
  write.csv(daily[[id]], file = file_name, row.names = FALSE)
}

# ------------------------------------------------------------
# 5. AGGREGATE AFTERWARDS (fast, in-memory)
# ------------------------------------------------------------
clean_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) return(NA_real_)
  round(m, 4)
}        # Helper function for mean aggregation and 4-digit rounding
agg_days <- function(df, n_days) {
  df %>%
    arrange(date) %>%
    mutate(group = as.numeric(date - min(date)) %/% n_days) %>%
    group_by(group) %>%
    summarise(
      Ini = min(date),
      End = max(date),
      Id = first(Id),
      across(c(n_pixels, n_valid, n_na, mean_w, sd, p2.5, p50, p97.5), clean_mean),
      .groups = "drop"
    ) %>%
    select(Ini, End, Id, everything(), -group)
} # 1. Function for N-day window aggregation
agg_monthly <- function(df) {
  df %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(Id, date) %>%
    summarise(
      across(c(n_pixels, n_valid, n_na, mean_w, sd, p2.5, p50, p97.5), clean_mean),
      .groups = "drop"
    ) %>%
    select(date, Id, everything())
}      # 2. Function for monthly aggregation (natural calendar months)
iwalk(daily, function(df, id) {
  write_csv(agg_days(df, 7),  paste0(out,"Agg_7d\\SM_", id, "_7day.csv"))
  write_csv(agg_days(df, 15), paste0(out,"Agg_15d\\SM_", id, "_15day.csv"))
  write_csv(agg_monthly(df),  paste0(out,"Agg_monthly\\SM_", id, "_monthly.csv"))
})   # 3. Process and export all data frames in the list







