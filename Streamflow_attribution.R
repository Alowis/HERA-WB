#############################################
# streamflow_attribution.R
# -------------------------------------------------
# 0️⃣ Packages & helpers
#############################################
library(stats)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)

`%!in%` <- Negate(`%in%`)

log_diff <- function(x) {
  if (any(x <= 0, na.rm = TRUE))
    stop("All values must be > 0 to compute log‑differences.")
  diff(log(x))
}

#############################################
# 1️⃣ Data preparation
#############################################
prepare_data <- function(file_path,
                         var_names = c(Q = "Q", P = "P", PET = "PET", S = "S"),
                         year_col = "year") {
  df_raw <- read_csv(file_path, col_types = cols())
  needed <- c(year_col, var_names)
  df_raw <- df_raw %>% select(all_of(needed)) %>% arrange(.data[[year_col]])
  df_diff <- df_raw %>%
    mutate(across(all_of(var_names), ~log_diff(.x), .names = "d_{col}")) %>%
    slice(-1) %>%
    mutate(!!year_col := .data[[year_col]][-1])
  return(df_diff)
}

#############################################
# 2️⃣ Elasticity estimation (stepwise)
#############################################
estimate_elasticities <- function(df_diff,
                                  response = "d_Q",
                                  predictors = c("d_P", "d_PET", "d_S"),
                                  r2_inc_thresh = 0.05,
                                  verbose = TRUE) {
  formula_full <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  formula_null <- as.formula(paste(response, "~ 1"))
  fit_null <- lm(formula_null, data = df_diff)
  
  step_fit <- step(fit_null,
                   scope = list(lower = formula_null, upper = formula_full),
                   direction = "both",
                   trace = ifelse(verbose, 1, 0),
                   k = log(nrow(df_diff)))   # BIC penalty
  
  final_terms <- attr(terms(step_fit), "term.labels")
  if (verbose) {
    cat("\nFinal model terms:", paste(final_terms, collapse = ", "), "\n")
    cat("R²:", summary(step_fit)$r.squared, "\n")
  }
  
  eps <- setNames(rep(0, length(predictors)), predictors)
  eps[final_terms] <- coef(step_fit)[final_terms]
  
  list(
    model = step_fit,
    elasticities = eps,
    final_predictors = final_terms,
    r2 = summary(step_fit)$r.squared
  )
}

#############################################
# 3️⃣ Predict streamflow change
#############################################
predict_Q_change <- function(elasticities,
                             Q_hist_mean,
                             X_hist_means,
                             X_fut_means) {
  delta_X <- X_fut_means - X_hist_means
  contrib <- elasticities * (delta_X / X_hist_means) * Q_hist_mean
  dQ_hat  <- sum(contrib)
  list(
    contribution = contrib,
    total_change = dQ_hat,
    fractional_change = dQ_hat / Q_hist_mean
  )
}

#############################################
# 4️⃣ Paired (or unpaired) t‑test
#############################################
compare_predictions <- function(pred_A, pred_B,
                                paired = TRUE, alpha = 0.05) {
  if (length(pred_A) != length(pred_B))
    stop("Vectors must be of equal length.")
  test_res <- t.test(pred_A, pred_B,
                     paired = paired,
                     alternative = "two.sided")
  list(
    t_statistic = test_res$statistic,
    df = test_res$parameter,
    p_value = test_res$p.value,
    conf_int = test_res$conf.int,
    reject_null = test_res$p.value < alpha,
    method = test_res$method
  )
}

#############################################
# 📊 Demo with synthetic data (no external files)
#############################################
set.seed(42)

# ---- 1) Simulate a 30‑year historical record (1985‑2014) ----
years <- 198
