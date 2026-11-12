# ============================================================
# fit_nonstationary_gpd_sensitivity_matrix.R
#
# Fast sensitivity matrix for the rainfall GPD benchmark.
#
# Specifications:
#   P95   / 3-day
#   P97.5 / 2-day
#   P97.5 / 3-day  (primary)
#   P97.5 / 5-day
#   P99   / 3-day
#
# For each specification fit:
#   M0: stationary GPD scale, pooled constant shape
#   M1: linear scale trend per decade, pooled constant shape
#   M2: three-period scale, pooled constant shape
#
# Assumes helper functions have been defined in:
#   R/models/fit_nonstationary_gpd_benchmark.R
#
# This script sources that file in FUNCTIONS-ONLY mode. Add the guard shown
# below to the benchmark script before its "Load data" section if not already
# present:
#
# if (identical(Sys.getenv("GPD_FUNCTIONS_ONLY"), "true")) {
#   return(invisible(TRUE))
# }
# ============================================================

library(data.table)
library(arrow)

analysis_start_year <- 1950L
analysis_end_year <- 2025L
centre_year <- mean(c(analysis_start_year, analysis_end_year))
accepted_rainfall_quality <- c(0L, 1L, 2L, 3L, 4L)

all_events_file <- file.path(
  "data-derived", "rainfall", "pot",
  "rainfall_pot_events_all_thresholds_runs_1950_2025.parquet"
)

helper_file <- file.path(
  "R",
  "models",
  "gpd_helpers.R"
)

output_dir <- file.path(
  "tables", "rainfall_models", "gpd"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(all_events_file)) {
  stop("Missing POT event file: ", all_events_file, call. = FALSE)
}

if (!file.exists(helper_file)) {
  stop("Missing non-stationary GPD benchmark script: ", helper_file, call. = FALSE)
}

# ------------------------------------------------------------
# Load functions without rerunning the primary benchmark
# ------------------------------------------------------------

source(
  helper_file
)

required_functions <- c(
  "fit_gpd_m0", "fit_gpd_m1", "fit_gpd_m2",
  "extract_m1_parameters", "shape_from_raw",
  "assign_climate_period"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Required helper functions were not loaded: ",
    paste(missing_functions, collapse = ", "),
    "\nAdd the GPD_FUNCTIONS_ONLY guard to the benchmark script before its data-loading section.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Load all POT specifications
# ------------------------------------------------------------

pot <- as.data.table(arrow::read_parquet(all_events_file))

pot <- pot[
  rainfall_quality_peak %in% accepted_rainfall_quality &
    is.finite(excess_mm) & excess_mm > 0 &
    is.finite(threshold_mm) & threshold_mm > 0
]

pot[, year := as.integer(format(event_peak_date, "%Y"))]
pot <- pot[year >= analysis_start_year & year <= analysis_end_year]
pot[, year_decade := (year - centre_year) / 10]
pot[, climate_period := assign_climate_period(year)]
pot[, normalized_excess := excess_mm / threshold_mm]

# ------------------------------------------------------------
# Compact sensitivity design
# ------------------------------------------------------------

sensitivity_design <- data.table(
  specification = c(
    "P95 / run3",
    "P97.5 / run2",
    "P97.5 / run3",
    "P97.5 / run5",
    "P99 / run3"
  ),
  threshold_id = c(
    "p95", "p975", "p975", "p975", "p99"
  ),
  run_length_days = c(
    3L, 2L, 3L, 5L, 3L
  ),
  primary = c(
    FALSE, FALSE, TRUE, FALSE, FALSE
  )
)

# ------------------------------------------------------------
# Helpers to extract sensitivity metrics
# ------------------------------------------------------------

extract_sensitivity_result <- function(
    specification,
    threshold_id,
    run_length_days,
    primary,
    x,
    m0,
    m1,
    m2) {

  m1_par <- extract_m1_parameters(m1)
  m1_time <- m1_par[term == "scale_time_per_decade"]

  beta_time <- as.numeric(m1_time$estimate)
  beta_low <- as.numeric(m1_time$conf_low)
  beta_high <- as.numeric(m1_time$conf_high)

  middle_log_diff <- as.numeric(
    m2$parameters_raw["period_1976_2000"]
  )

  recent_log_diff <- as.numeric(
    m2$parameters_raw["period_2001_2025"]
  )

  shape_m0 <- shape_from_raw(
    m0$parameters_raw["raw_shape"]
  )

  shape_m1 <- shape_from_raw(
    m1$parameters_raw["raw_shape"]
  )

  shape_m2 <- shape_from_raw(
    m2$parameters_raw["raw_shape"]
  )

  lr_m1 <- 2 * (m1$log_likelihood - m0$log_likelihood)
  lr_m2 <- 2 * (m2$log_likelihood - m0$log_likelihood)

  data.table(
    specification = specification,
    threshold_id = threshold_id,
    run_length_days = as.integer(run_length_days),
    primary = primary,
    events = as.integer(nrow(x)),
    stations = as.integer(uniqueN(x$station_number)),

    m0_valid = m0$fit_valid,
    m1_valid = m1$fit_valid,
    m2_valid = m2$fit_valid,

    shape_m0 = as.numeric(shape_m0),
    shape_m1 = as.numeric(shape_m1),
    shape_m2 = as.numeric(shape_m2),

    m1_log_scale_change_per_decade = beta_time,
    m1_percent_scale_change_per_decade =
      100 * (exp(beta_time) - 1),
    m1_percent_change_conf_low =
      100 * (exp(beta_low) - 1),
    m1_percent_change_conf_high =
      100 * (exp(beta_high) - 1),
    m1_interval_excludes_zero =
      beta_low > 0 || beta_high < 0,

    m1_delta_AIC_vs_m0 = as.numeric(m1$AIC - m0$AIC),
    m1_LR_vs_m0 = as.numeric(lr_m1),
    m1_LR_df = as.integer(m1$parameter_count - m0$parameter_count),
    m1_LR_p = as.numeric(
      pchisq(
        lr_m1,
        df = m1$parameter_count - m0$parameter_count,
        lower.tail = FALSE
      )
    ),

    m2_middle_percent_vs_early =
      100 * (exp(middle_log_diff) - 1),
    m2_recent_percent_vs_early =
      100 * (exp(recent_log_diff) - 1),

    m2_delta_AIC_vs_m0 = as.numeric(m2$AIC - m0$AIC),
    m2_LR_vs_m0 = as.numeric(lr_m2),
    m2_LR_df = as.integer(m2$parameter_count - m0$parameter_count),
    m2_LR_p = as.numeric(
      pchisq(
        lr_m2,
        df = m2$parameter_count - m0$parameter_count,
        lower.tail = FALSE
      )
    ),

    m0_AIC = as.numeric(m0$AIC),
    m1_AIC = as.numeric(m1$AIC),
    m2_AIC = as.numeric(m2$AIC),

    m0_max_gradient = as.numeric(m0$maximum_gradient),
    m1_max_gradient = as.numeric(m1$maximum_gradient),
    m2_max_gradient = as.numeric(m2$maximum_gradient)
  )
}

# ------------------------------------------------------------
# Fit the five specifications
# ------------------------------------------------------------

result_list <- vector("list", nrow(sensitivity_design))
model_list <- vector("list", nrow(sensitivity_design))

for (i in seq_len(nrow(sensitivity_design))) {

  spec <- sensitivity_design[i]

  message(
    "\n========================================",
    "\nSensitivity specification ", i, " of ", nrow(sensitivity_design),
    ": ", spec$specification,
    "\n========================================"
  )

  x <- pot[
    threshold_id == spec$threshold_id &
      run_length_days == spec$run_length_days
  ]

  if (nrow(x) < 1000L) {
    stop(
      "Unexpectedly small event sample for ",
      spec$specification,
      ": ",
      nrow(x),
      call. = FALSE
    )
  }

  message("  Fitting M0")
  m0 <- fit_gpd_m0(
    excess = x$normalized_excess
  )

  message("  Fitting M1")
  m1 <- fit_gpd_m1(
    excess = x$normalized_excess,
    year_decade = x$year_decade
  )

  message("  Fitting M2")
  m2 <- fit_gpd_m2(
    excess = x$normalized_excess,
    climate_period = x$climate_period
  )

  result_list[[i]] <- extract_sensitivity_result(
    specification = spec$specification,
    threshold_id = spec$threshold_id,
    run_length_days = spec$run_length_days,
    primary = spec$primary,
    x = x,
    m0 = m0,
    m1 = m1,
    m2 = m2
  )

  model_list[[i]] <- list(
    M0 = m0,
    M1 = m1,
    M2 = m2
  )

  names(model_list)[i] <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    spec$specification
  )
}

sensitivity_results <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

# Preserve requested order.
sensitivity_results[, specification_order := match(
  specification,
  sensitivity_design$specification
)]
setorder(sensitivity_results, specification_order)
sensitivity_results[, specification_order := NULL]

# ------------------------------------------------------------
# Robustness flags
# ------------------------------------------------------------

sensitivity_results[, no_linear_trend :=
  !m1_interval_excludes_zero &
    m1_LR_p >= 0.05
]

sensitivity_results[, m2_preferred_over_m0 :=
  m2_delta_AIC_vs_m0 < -2 &
    m2_LR_p < 0.05
]

sensitivity_results[, mid_period_lower :=
  m2_middle_percent_vs_early < 0
]

sensitivity_results[, recent_nearer_early_than_middle :=
  abs(m2_recent_percent_vs_early) <
    abs(m2_middle_percent_vs_early)
]

sensitivity_results[, all_models_valid :=
  m0_valid & m1_valid & m2_valid
]

robustness_summary <- data.table(
  specifications = nrow(sensitivity_results),
  specifications_all_models_valid = sum(sensitivity_results$all_models_valid),
  specifications_with_no_linear_trend = sum(sensitivity_results$no_linear_trend),
  specifications_preferring_m2 = sum(sensitivity_results$m2_preferred_over_m0),
  specifications_with_lower_middle_period = sum(sensitivity_results$mid_period_lower),
  specifications_with_recent_recovery = sum(
    sensitivity_results$recent_nearer_early_than_middle
  ),
  minimum_linear_change_pct_decade = min(
    sensitivity_results$m1_percent_scale_change_per_decade
  ),
  maximum_linear_change_pct_decade = max(
    sensitivity_results$m1_percent_scale_change_per_decade
  ),
  minimum_shape = min(sensitivity_results$shape_m0),
  maximum_shape = max(sensitivity_results$shape_m0)
)

# ------------------------------------------------------------
# Compact publication-style table
# ------------------------------------------------------------

sensitivity_compact <- sensitivity_results[, .(
  specification,
  events,
  stations,
  shape = shape_m0,

  linear_change_pct_decade =
    m1_percent_scale_change_per_decade,
  linear_ci_low =
    m1_percent_change_conf_low,
  linear_ci_high =
    m1_percent_change_conf_high,
  linear_LR_p =
    m1_LR_p,
  linear_delta_AIC =
    m1_delta_AIC_vs_m0,

  middle_vs_early_pct =
    m2_middle_percent_vs_early,
  recent_vs_early_pct =
    m2_recent_percent_vs_early,
  period_LR_p =
    m2_LR_p,
  period_delta_AIC =
    m2_delta_AIC_vs_m0,

  all_models_valid,
  no_linear_trend,
  m2_preferred_over_m0
)]

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  list(
    design = sensitivity_design,
    results = sensitivity_results,
    models = model_list,
    robustness_summary = robustness_summary
  ),
  file.path(
    "data-derived", "rainfall", "models", "gpd",
    "nonstationary_gpd_sensitivity_matrix.rds"
  ),
  compress = "xz"
)

fwrite(
  sensitivity_results,
  file.path(
    output_dir,
    "nonstationary_gpd_sensitivity_matrix.csv"
  )
)

fwrite(
  sensitivity_compact,
  file.path(
    output_dir,
    "nonstationary_gpd_sensitivity_matrix_compact.csv"
  )
)

fwrite(
  robustness_summary,
  file.path(
    output_dir,
    "nonstationary_gpd_sensitivity_robustness_summary.csv"
  )
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

message(
  "\n========================================",
  "\nNON-STATIONARY GPD SENSITIVITY COMPLETE",
  "\n========================================"
)

print(sensitivity_compact)
print(robustness_summary)
