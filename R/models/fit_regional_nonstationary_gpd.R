# ============================================================
# fit_regional_nonstationary_gpd.R
#
# Lightweight regional non-stationary GPD analysis.
#
# Primary event definition:
#   - station-specific P97.5 wet-day threshold
#   - 3-day runs declustering
#   - PRCP_QUAL codes 0-4
#   - 1950-2025
#
# Fits separately by state/territory:
#   M0: stationary GPD scale, constant regional shape
#   M1: linear scale trend per decade, constant regional shape
#   M2: three-period scale, constant regional shape
#
# No spatial random fields, MCMC, bootstrap, or station-by-station fits.
# Reuses tested functions from fit_nonstationary_gpd_benchmark.R.
# ============================================================

library(data.table)
library(arrow)

# ============================================================
# Configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L
centre_year <- mean(c(analysis_start_year, analysis_end_year))

accepted_rainfall_quality <- c(0L, 1L, 2L, 3L, 4L)

minimum_events_region <- 500L
minimum_events_period <- 100L
minimum_stations_region <- 5L

pot_file <- file.path(
  "data-derived", "rainfall", "pot",
  "rainfall_pot_events_p975_run3_1950_2025.parquet"
)

helper_file <- file.path(
  "R",
  "models",
  "gpd_helpers.R"
)

model_output_dir <- file.path(
  "data-derived", "rainfall", "models", "gpd"
)

table_output_dir <- file.path(
  "tables", "rainfall_models", "gpd"
)

dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pot_file)) {
  stop("Missing primary POT event file: ", pot_file, call. = FALSE)
}

if (!file.exists(helper_file)) {
  stop("Missing GPD helper script: ", helper_file, call. = FALSE)
}

# ============================================================
# Load helper functions only
# ============================================================

source(
  helper_file
)

required_functions <- c(
  "fit_gpd_m0",
  "fit_gpd_m1",
  "fit_gpd_m2",
  "extract_m1_parameters",
  "shape_from_raw",
  "assign_climate_period"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Missing helper functions: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}

# ============================================================
# Load and prepare primary events
# ============================================================

events <- as.data.table(arrow::read_parquet(pot_file))

required_columns <- c(
  "station_number", "station_name", "state",
  "event_peak_date", "threshold_mm", "excess_mm",
  "rainfall_quality_peak"
)

missing_columns <- setdiff(required_columns, names(events))

if (length(missing_columns) > 0L) {
  stop(
    "Primary POT events are missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

events <- events[
  rainfall_quality_peak %in% accepted_rainfall_quality &
    is.finite(excess_mm) & excess_mm > 0 &
    is.finite(threshold_mm) & threshold_mm > 0 &
    !is.na(state) & nzchar(state)
]

events[, year := as.integer(format(event_peak_date, "%Y"))]
events <- events[year >= analysis_start_year & year <= analysis_end_year]
events[, year_decade := (year - centre_year) / 10]
events[, climate_period := assign_climate_period(year)]
events[, normalized_excess := excess_mm / threshold_mm]

# Combine ACT with NSW if ACT exists as a separate state code. This avoids a
# very small regional fit while preserving all events.
events[state == "ACT", state := "NSW"]

# ============================================================
# Regional sample audit and eligibility
# ============================================================

regional_sample_audit <- events[, .(
  events = as.integer(.N),
  stations = as.integer(uniqueN(station_number)),
  first_year = as.integer(min(year)),
  last_year = as.integer(max(year)),
  events_1950_1975 = as.integer(sum(climate_period == "1950-1975")),
  events_1976_2000 = as.integer(sum(climate_period == "1976-2000")),
  events_2001_2025 = as.integer(sum(climate_period == "2001-2025"))
), by = state]

regional_sample_audit[, eligible :=
  events >= minimum_events_region &
    stations >= minimum_stations_region &
    events_1950_1975 >= minimum_events_period &
    events_1976_2000 >= minimum_events_period &
    events_2001_2025 >= minimum_events_period
]

setorder(regional_sample_audit, state)

regions_to_fit <- regional_sample_audit[eligible == TRUE, state]

if (length(regions_to_fit) == 0L) {
  stop("No regions meet the minimum sample criteria.", call. = FALSE)
}

# ============================================================
# Extract compact regional result
# ============================================================

extract_regional_result <- function(region, x, m0, m1, m2) {

  m1_parameters <- extract_m1_parameters(m1)
  m1_time <- m1_parameters[term == "scale_time_per_decade"]

  beta_time <- as.numeric(m1_time$estimate)
  beta_low <- as.numeric(m1_time$conf_low)
  beta_high <- as.numeric(m1_time$conf_high)

  middle_log_difference <- as.numeric(
    m2$parameters_raw["period_1976_2000"]
  )

  recent_log_difference <- as.numeric(
    m2$parameters_raw["period_2001_2025"]
  )

  shape_m0 <- as.numeric(
    shape_from_raw(m0$parameters_raw["raw_shape"])
  )

  shape_m1 <- as.numeric(
    shape_from_raw(m1$parameters_raw["raw_shape"])
  )

  shape_m2 <- as.numeric(
    shape_from_raw(m2$parameters_raw["raw_shape"])
  )

  lr_m1 <- 2 * (m1$log_likelihood - m0$log_likelihood)
  lr_m2 <- 2 * (m2$log_likelihood - m0$log_likelihood)

  data.table(
    region = region,
    events = as.integer(nrow(x)),
    stations = as.integer(uniqueN(x$station_number)),

    m0_valid = m0$fit_valid,
    m1_valid = m1$fit_valid,
    m2_valid = m2$fit_valid,

    shape_m0 = shape_m0,
    shape_m1 = shape_m1,
    shape_m2 = shape_m2,

    linear_change_pct_decade = 100 * (exp(beta_time) - 1),
    linear_ci_low = 100 * (exp(beta_low) - 1),
    linear_ci_high = 100 * (exp(beta_high) - 1),
    linear_interval_excludes_zero = beta_low > 0 || beta_high < 0,

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

    middle_vs_early_pct = 100 * (exp(middle_log_difference) - 1),
    recent_vs_early_pct = 100 * (exp(recent_log_difference) - 1),

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

# ============================================================
# Fit M0, M1, M2 by region
# ============================================================

regional_result_list <- vector("list", length(regions_to_fit))
regional_model_list <- vector("list", length(regions_to_fit))

for (i in seq_along(regions_to_fit)) {

  region_value <- regions_to_fit[i]

  message(
    "\n========================================",
    "\nRegional GPD ", i, " of ", length(regions_to_fit),
    ": ", region_value,
    "\n========================================"
  )

  x <- events[state == region_value]

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

  regional_result_list[[i]] <- extract_regional_result(
    region = region_value,
    x = x,
    m0 = m0,
    m1 = m1,
    m2 = m2
  )

  regional_model_list[[i]] <- list(
    M0 = m0,
    M1 = m1,
    M2 = m2
  )

  names(regional_model_list)[i] <- region_value
}

regional_results <- rbindlist(
  regional_result_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(regional_results, region)

# ============================================================
# Interpretation flags
# ============================================================

regional_results[, all_models_valid :=
  m0_valid & m1_valid & m2_valid
]

regional_results[, no_linear_trend :=
  !linear_interval_excludes_zero &
    m1_LR_p >= 0.05
]

regional_results[, linear_trend_preferred :=
  linear_interval_excludes_zero &
    m1_delta_AIC_vs_m0 < -2 &
    m1_LR_p < 0.05
]

regional_results[, m2_preferred_over_m0 :=
  m2_delta_AIC_vs_m0 < -2 &
    m2_LR_p < 0.05
]

regional_results[, mid_period_lower :=
  middle_vs_early_pct < 0
]

regional_results[, recent_recovery :=
  abs(recent_vs_early_pct) < abs(middle_vs_early_pct)
]

regional_results[, recent_above_early :=
  recent_vs_early_pct > 0
]

# ============================================================
# Compact results and cross-region summaries
# ============================================================

regional_compact <- regional_results[, .(
  region,
  events,
  stations,
  shape = shape_m0,

  linear_change_pct_decade,
  linear_ci_low,
  linear_ci_high,
  m1_LR_p,
  m1_delta_AIC_vs_m0,

  middle_vs_early_pct,
  recent_vs_early_pct,
  m2_LR_p,
  m2_delta_AIC_vs_m0,

  all_models_valid,
  no_linear_trend,
  linear_trend_preferred,
  m2_preferred_over_m0,
  mid_period_lower,
  recent_recovery,
  recent_above_early
)]

regional_summary <- data.table(
  regions_fitted = nrow(regional_results),
  regions_all_models_valid = sum(regional_results$all_models_valid),
  regions_no_linear_trend = sum(regional_results$no_linear_trend),
  regions_with_preferred_linear_trend = sum(
    regional_results$linear_trend_preferred
  ),
  regions_preferring_m2 = sum(regional_results$m2_preferred_over_m0),
  regions_with_lower_middle_period = sum(regional_results$mid_period_lower),
  regions_with_recent_recovery = sum(regional_results$recent_recovery),
  regions_recent_above_early = sum(regional_results$recent_above_early),
  minimum_shape = min(regional_results$shape_m0),
  maximum_shape = max(regional_results$shape_m0),
  minimum_linear_change_pct_decade = min(
    regional_results$linear_change_pct_decade
  ),
  maximum_linear_change_pct_decade = max(
    regional_results$linear_change_pct_decade
  )
)

# ============================================================
# Save outputs
# ============================================================

saveRDS(
  list(
    sample_audit = regional_sample_audit,
    results = regional_results,
    models = regional_model_list,
    configuration = list(
      analysis_start_year = analysis_start_year,
      analysis_end_year = analysis_end_year,
      centre_year = centre_year,
      threshold = "P97.5",
      run_length_days = 3L,
      accepted_rainfall_quality = accepted_rainfall_quality,
      minimum_events_region = minimum_events_region,
      minimum_events_period = minimum_events_period,
      minimum_stations_region = minimum_stations_region
    )
  ),
  file.path(
    model_output_dir,
    "regional_nonstationary_gpd_models.rds"
  ),
  compress = "xz"
)

fwrite(
  regional_sample_audit,
  file.path(
    table_output_dir,
    "regional_nonstationary_gpd_sample_audit.csv"
  )
)

fwrite(
  regional_results,
  file.path(
    table_output_dir,
    "regional_nonstationary_gpd_results.csv"
  )
)

fwrite(
  regional_compact,
  file.path(
    table_output_dir,
    "regional_nonstationary_gpd_results_compact.csv"
  )
)

fwrite(
  regional_summary,
  file.path(
    table_output_dir,
    "regional_nonstationary_gpd_summary.csv"
  )
)

# ============================================================
# Console output
# ============================================================

message(
  "\n========================================",
  "\nREGIONAL NON-STATIONARY GPD COMPLETE",
  "\n========================================"
)

print(regional_sample_audit)
print(regional_compact)
print(regional_summary)
