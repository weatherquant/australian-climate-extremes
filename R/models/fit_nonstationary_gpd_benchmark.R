# ============================================================
# fit_nonstationary_gpd_benchmark.R
#
# Lightweight non-stationary GPD benchmark for Australian
# extreme daily rainfall.
#
# Primary event definition:
#   Threshold:       station-specific P97.5
#   Declustering:    3-day runs
#   Quality:         PRCP_QUAL codes 0-4
#   Period:          1950-2025
#
# Models:
#   M0 stationary:
#      log(scale_ratio) = beta_0
#      shape = constant
#
#   M1 linear time trend in scale:
#      log(scale_ratio_t) = beta_0 + beta_1 * year_decade
#      shape = constant
#
#   M2 categorical climate-period effects in scale:
#      log(scale_ratio_t) = beta_0 + period effects
#      shape = constant
#
# year_decade = (year - 1987.5) / 10, so beta_1 is the
# log multiplicative change in GPD scale per decade.
#
# Computational design:
#   - no mesh
#   - no random field
#   - no MCMC
#   - no bootstrap in this first benchmark
#   - only 2, 3, or 4 fitted parameters
#
# Return levels combine:
#   - station-specific P97.5 threshold
#   - station-specific observed event rate over 1950-2025
#   - modelled time-varying GPD scale
#   - pooled constant GPD shape
# ============================================================


# ============================================================
# Packages and reusable GPD functions
# ============================================================

library(data.table)
library(arrow)

source(
  "R/models/gpd_helpers.R"
)


# ============================================================
# Analysis configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L

analysis_years <-
  analysis_end_year -
  analysis_start_year +
  1L

centre_year <- mean(
  c(
    analysis_start_year,
    analysis_end_year
  )
)

primary_threshold_id <- "p975"
primary_run_length <- 3L

accepted_rainfall_quality <- c(
  0L,
  1L,
  2L,
  3L,
  4L
)

representative_years <- c(
  1960L,
  1990L,
  2020L
)

return_periods_years <- c(
  10,
  20,
  50,
  100
)

pot_file <- file.path(
  "data-derived",
  "rainfall",
  "pot",
  "rainfall_pot_events_p975_run3_1950_2025.parquet"
)

model_output_dir <- file.path(
  "data-derived",
  "rainfall",
  "models",
  "gpd"
)

table_output_dir <- file.path(
  "tables",
  "rainfall_models",
  "gpd"
)

dir.create(
  model_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load data
# ============================================================

if (!file.exists(pot_file)) {
  stop(
    "Primary P97.5 POT event file does not exist: ",
    pot_file,
    call. = FALSE
  )
}

events <- as.data.table(
  arrow::read_parquet(pot_file)
)

required_columns <- c(
  "station_number",
  "station_name",
  "state",
  "latitude",
  "longitude",
  "event_peak_date",
  "event_peak_mm",
  "threshold_mm",
  "excess_mm",
  "rainfall_quality_peak"
)

missing_columns <- setdiff(required_columns, names(events))

if (length(missing_columns) > 0L) {
  stop(
    "POT event data are missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

events <- events[
  rainfall_quality_peak %in% accepted_rainfall_quality &
    is.finite(excess_mm) &
    excess_mm > 0 &
    is.finite(threshold_mm) &
    threshold_mm > 0
]

events[, year := as.integer(format(event_peak_date, "%Y"))]

events <- events[
  year >= analysis_start_year &
    year <= analysis_end_year
]

events[, year_decade := (year - centre_year) / 10]

events[, climate_period := assign_climate_period(year)]

events[, normalized_excess := excess_mm / threshold_mm]

if (nrow(events) != 84351L) {
  warning(
    "Expected 84,351 primary events but found ",
    nrow(events),
    ". Continue only after confirming the QC/event definition.",
    call. = FALSE
  )
}


# ============================================================
# Fit M0, M1, and M2
# ============================================================

message("Fitting M0: stationary GPD")

fit_m0 <- fit_gpd_m0(
  excess = events$normalized_excess
)

message("Fitting M1: linear time trend in GPD scale")

fit_m1 <- fit_gpd_m1(
  excess = events$normalized_excess,
  year_decade = events$year_decade
)

message("Fitting M2: three-period GPD scale")

fit_m2 <- fit_gpd_m2(
  excess = events$normalized_excess,
  climate_period = events$climate_period
)


# ============================================================
# Model comparison
# ============================================================

model_comparison <- rbindlist(
  lapply(
    list(fit_m0, fit_m1, fit_m2),
    function(fit) {
      data.table(
        model = fit$model_name,
        observations = fit$n,
        parameters = fit$parameter_count,
        log_likelihood = fit$log_likelihood,
        negative_log_likelihood = fit$negative_log_likelihood,
        AIC = fit$AIC,
        convergence = fit$convergence,
        maximum_gradient = fit$maximum_gradient,
        hessian_positive_definite = fit$hessian_positive_definite,
        fit_valid = fit$fit_valid
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

model_comparison[, delta_AIC := AIC - min(AIC)]
setorder(model_comparison, AIC)

lrt_m0_m1 <- data.table(
  comparison = "M0 stationary vs M1 linear scale trend",
  likelihood_ratio = 2 * (fit_m1$log_likelihood - fit_m0$log_likelihood),
  degrees_of_freedom = fit_m1$parameter_count - fit_m0$parameter_count
)

lrt_m0_m1[, p_value := pchisq(
  likelihood_ratio,
  df = degrees_of_freedom,
  lower.tail = FALSE
)]

lrt_m0_m2 <- data.table(
  comparison = "M0 stationary vs M2 three-period scale",
  likelihood_ratio = 2 * (fit_m2$log_likelihood - fit_m0$log_likelihood),
  degrees_of_freedom = fit_m2$parameter_count - fit_m0$parameter_count
)

lrt_m0_m2[, p_value := pchisq(
  likelihood_ratio,
  df = degrees_of_freedom,
  lower.tail = FALSE
)]

likelihood_ratio_tests <- rbindlist(
  list(lrt_m0_m1, lrt_m0_m2),
  use.names = TRUE
)


# ============================================================
# Interpretable parameter tables
# ============================================================

parameters_m0 <- extract_m0_parameters(fit_m0)
parameters_m1 <- extract_m1_parameters(fit_m1)
parameters_m2 <- extract_m2_parameters(fit_m2)

parameter_estimates <- rbindlist(
  list(
    parameters_m0,
    parameters_m1,
    parameters_m2
  ),
  use.names = TRUE,
  fill = TRUE
)

m1_time <- parameters_m1[
  term == "scale_time_per_decade"
]

m1_time_effect <- data.table(
  estimate_log_scale_per_decade = m1_time$estimate,
  standard_error = m1_time$std_error,
  conf_low_log_scale = m1_time$conf_low,
  conf_high_log_scale = m1_time$conf_high,
  scale_ratio_per_decade = exp(m1_time$estimate),
  percent_scale_change_per_decade = 100 * (exp(m1_time$estimate) - 1),
  percent_change_conf_low = 100 * (exp(m1_time$conf_low) - 1),
  percent_change_conf_high = 100 * (exp(m1_time$conf_high) - 1),
  interval_excludes_zero =
    m1_time$conf_low > 0 ||
    m1_time$conf_high < 0
)

m2_period_comparison <- data.table(
  comparison = c(
    "1976-2000 vs 1950-1975",
    "2001-2025 vs 1950-1975"
  ),
  log_scale_difference = c(
    fit_m2$parameters_raw["period_1976_2000"],
    fit_m2$parameters_raw["period_2001_2025"]
  )
)

m2_period_comparison[, percent_scale_difference :=
  100 * (exp(log_scale_difference) - 1)
]


# ============================================================
# Time-varying scale under M1
# ============================================================

m1_scale_by_year <- data.table(
  year = seq.int(
    analysis_start_year,
    analysis_end_year
  )
)

m1_scale_by_year[, year_decade := (year - centre_year) / 10]

m1_scale_by_year[, scale_ratio := exp(
  fit_m1$parameters_raw["log_scale_intercept"] +
    fit_m1$parameters_raw["scale_time_per_decade"] * year_decade
)]

m1_shape <- shape_from_raw(
  fit_m1$parameters_raw["raw_shape"]
)

representative_scale <- m1_scale_by_year[
  year %in% representative_years
]

representative_scale[, shape := m1_shape]


# ============================================================
# Station event rates
# ============================================================

station_event_rates <- events[, .(
  events = as.integer(.N),
  event_rate_per_year = as.numeric(.N / analysis_years),
  threshold_mm = as.numeric(first(threshold_mm))
), by = .(
  station_number,
  station_name,
  state,
  latitude,
  longitude
)]


# ============================================================
# M1 return levels by representative year
# ============================================================

station_return_levels_m1 <- rbindlist(
  lapply(
    representative_years,
    function(target_year) {

      target_scale_ratio <- m1_scale_by_year[
        year == target_year,
        scale_ratio
      ]

      rbindlist(
        lapply(
          return_periods_years,
          function(return_period) {

            station_event_rates[, .(
              station_number,
              station_name,
              state,
              latitude,
              longitude,
              events,
              event_rate_per_year,
              threshold_mm,
              model_year = as.integer(target_year),
              scale_ratio = as.numeric(target_scale_ratio),
              shape = as.numeric(m1_shape),
              return_period_years = as.numeric(return_period),
              return_level_mm = pot_return_level(
                threshold_mm = threshold_mm,
                scale_ratio = target_scale_ratio,
                shape = m1_shape,
                event_rate_per_year = event_rate_per_year,
                return_period_years = return_period
              )
            )]
          }
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

return_level_summary_m1 <- station_return_levels_m1[
  is.finite(return_level_mm),
  .(
    stations = as.integer(.N),
    median_return_level_mm = as.numeric(median(return_level_mm)),
    q10_return_level_mm = as.numeric(quantile(return_level_mm, .10, names = FALSE)),
    q90_return_level_mm = as.numeric(quantile(return_level_mm, .90, names = FALSE))
  ),
  by = .(
    model_year,
    return_period_years
  )
]

setorder(
  return_level_summary_m1,
  return_period_years,
  model_year
)


# ============================================================
# Return-level changes: 1960 to 2020
# ============================================================

return_1960 <- return_level_summary_m1[
  model_year == 1960L,
  .(
    return_period_years,
    median_1960_mm = median_return_level_mm
  )
]

return_2020 <- return_level_summary_m1[
  model_year == 2020L,
  .(
    return_period_years,
    median_2020_mm = median_return_level_mm
  )
]

return_level_change_1960_2020 <- merge(
  return_1960,
  return_2020,
  by = "return_period_years",
  all = TRUE
)

return_level_change_1960_2020[, `:=`(
  absolute_change_mm = median_2020_mm - median_1960_mm,
  percent_change = 100 * (median_2020_mm / median_1960_mm - 1)
)]


# ============================================================
# Validation audit
# ============================================================

benchmark_audit <- data.table(
  analysis_start_year = analysis_start_year,
  analysis_end_year = analysis_end_year,
  centre_year = centre_year,
  events = nrow(events),
  stations = uniqueN(events$station_number),
  accepted_quality_codes = paste(
    accepted_rainfall_quality,
    collapse = ","
  ),
  m0_valid = fit_m0$fit_valid,
  m1_valid = fit_m1$fit_valid,
  m2_valid = fit_m2$fit_valid
)

if (!all(benchmark_audit[, c(m0_valid, m1_valid, m2_valid)])) {
  warning(
    "At least one non-stationary GPD benchmark did not pass all fit checks.",
    call. = FALSE
  )
}


# ============================================================
# Save outputs
# ============================================================

saveRDS(
  list(
    M0 = fit_m0,
    M1 = fit_m1,
    M2 = fit_m2,
    configuration = list(
      analysis_start_year = analysis_start_year,
      analysis_end_year = analysis_end_year,
      centre_year = centre_year,
      threshold_id = primary_threshold_id,
      run_length = primary_run_length,
      accepted_rainfall_quality = accepted_rainfall_quality
    )
  ),
  file.path(
    model_output_dir,
    "nonstationary_gpd_benchmark_models.rds"
  ),
  compress = "xz"
)

fwrite(
  benchmark_audit,
  file.path(
    table_output_dir,
    "nonstationary_gpd_benchmark_audit.csv"
  )
)

fwrite(
  model_comparison,
  file.path(
    table_output_dir,
    "nonstationary_gpd_model_comparison.csv"
  )
)

fwrite(
  likelihood_ratio_tests,
  file.path(
    table_output_dir,
    "nonstationary_gpd_likelihood_ratio_tests.csv"
  )
)

fwrite(
  parameter_estimates,
  file.path(
    table_output_dir,
    "nonstationary_gpd_parameter_estimates.csv"
  )
)

fwrite(
  m1_time_effect,
  file.path(
    table_output_dir,
    "nonstationary_gpd_linear_time_effect.csv"
  )
)

fwrite(
  m2_period_comparison,
  file.path(
    table_output_dir,
    "nonstationary_gpd_period_comparison.csv"
  )
)

fwrite(
  m1_scale_by_year,
  file.path(
    table_output_dir,
    "nonstationary_gpd_m1_scale_by_year.csv"
  )
)

fwrite(
  representative_scale,
  file.path(
    table_output_dir,
    "nonstationary_gpd_representative_scale.csv"
  )
)

fwrite(
  station_return_levels_m1,
  file.path(
    table_output_dir,
    "nonstationary_gpd_station_return_levels.csv"
  )
)

fwrite(
  return_level_summary_m1,
  file.path(
    table_output_dir,
    "nonstationary_gpd_return_level_summary.csv"
  )
)

fwrite(
  return_level_change_1960_2020,
  file.path(
    table_output_dir,
    "nonstationary_gpd_return_level_change_1960_2020.csv"
  )
)


# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nNON-STATIONARY GPD BENCHMARK COMPLETE",
  "\n========================================"
)

print(benchmark_audit)
print(model_comparison)
print(likelihood_ratio_tests)
print(parameter_estimates)
print(m1_time_effect)
print(m2_period_comparison)
print(representative_scale)
print(return_level_summary_m1)
print(return_level_change_1960_2020)
