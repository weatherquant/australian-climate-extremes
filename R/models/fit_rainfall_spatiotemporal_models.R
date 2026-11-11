# ============================================================
# fit_rainfall_spatiotemporal_models.R
#
# Spatial and spatiotemporal inferential models for the frozen
# 1950-2025 continental Australian rainfall network.
#
# Responses:
#   1. Wet-day occurrence: binomial proportion
#   2. P95 exceedance frequency: negative binomial with exposure
#   3. P99 exceedance frequency: negative binomial with exposure
#   4. Mean wet-day intensity: Gamma with log link
#   5. Rx1day: Gaussian model for log(Rx1day)
#
# Each response is fitted with:
#   M3: spatial field only
#   M4: spatial field + AR(1) spatiotemporal field
#
# Coordinates:
#   EPSG:3577 Australian Albers, converted to kilometres.
#
# Time coefficient:
#   year_decade = (year - 1987.5) / 10
#   Therefore the fixed-effect coefficient is change per decade.
#
# This script is restartable. Existing model RDS files are loaded
# unless overwrite_models is set to TRUE.
# ============================================================


# ============================================================
# Configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L
annual_completeness_threshold <- 0.95
mesh_cutoff_km <- 30
overwrite_models <- FALSE

metrics_file <- file.path(
  "data-derived",
  "rainfall",
  "extremes",
  "rainfall_annual_extreme_metrics_1950_2025.parquet"
)

model_output_dir <- file.path("data-derived", "rainfall", "models", "spatiotemporal")

summary_output_dir <- file.path("tables", "rainfall_models", "spatiotemporal")


# ============================================================
# Packages
# ============================================================

required_packages <- c("data.table", "arrow", "sf", "sdmTMB")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

library(data.table)
library(arrow)
library(sf)
library(sdmTMB)


dir.create(model_output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(summary_output_dir,
           recursive = TRUE,
           showWarnings = FALSE)


# ============================================================
# Helpers
# ============================================================

safe_aic <- function(model) {
  value <- tryCatch(
    stats::AIC(model),
    error = function(e)
      NA_real_
  )
  as.numeric(value)
}


safe_loglik <- function(model) {
  value <- tryCatch(
    stats::logLik(model),
    error = function(e)
      NA
  )
  
  if (length(value) == 0L || is.na(value[[1]])) {
    return(NA_real_)
  }
  
  as.numeric(value)
}


extract_fixed_effects <- function(model, model_name, response_name) {
  coefficient_table <- tryCatch(
    summary(model)$coefficients,
    error = function(e)
      NULL
  )
  
  if (is.null(coefficient_table)) {
    return(
      data.table(
        model = model_name,
        response = response_name,
        term = character(),
        estimate = numeric(),
        standard_error = numeric(),
        statistic = numeric(),
        p_value = numeric(),
        conf_low = numeric(),
        conf_high = numeric()
      )
    )
  }
  
  coefficient_table <- as.data.frame(coefficient_table)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  
  column_names <- tolower(names(coefficient_table))
  
  estimate_column <- which(column_names %in% c("estimate", "est"))[1]
  se_column <- which(grepl("std|se", column_names))[1]
  statistic_column <- which(grepl("z value|t value|statistic", column_names))[1]
  p_column <- which(grepl("pr\\(|p.value|p_value", column_names))[1]
  
  estimate <- if (!is.na(estimate_column)) {
    as.numeric(coefficient_table[[estimate_column]])
  } else {
    rep(NA_real_, nrow(coefficient_table))
  }
  
  standard_error <- if (!is.na(se_column)) {
    as.numeric(coefficient_table[[se_column]])
  } else {
    rep(NA_real_, nrow(coefficient_table))
  }
  
  statistic <- if (!is.na(statistic_column)) {
    as.numeric(coefficient_table[[statistic_column]])
  } else {
    estimate / standard_error
  }
  
  p_value <- if (!is.na(p_column)) {
    as.numeric(coefficient_table[[p_column]])
  } else {
    2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }
  
  data.table(
    model = model_name,
    response = response_name,
    term = as.character(coefficient_table$term),
    estimate = estimate,
    standard_error = standard_error,
    statistic = statistic,
    p_value = p_value,
    conf_low = estimate - 1.96 * standard_error,
    conf_high = estimate + 1.96 * standard_error
  )
}


fit_or_load_sdmTMB <- function(model_name,
                               formula,
                               data,
                               mesh,
                               family,
                               spatial = "on",
                               spatiotemporal = "off",
                               time = NULL,
                               offset = NULL,
                               weights = NULL,
                               overwrite = FALSE) {
  output_file <- file.path(model_output_dir, paste0(model_name, ".rds"))
  
  if (file.exists(output_file) && !overwrite) {
    message("Loading existing model: ", model_name)
    return(readRDS(output_file))
  }
  
  message(
    "\n========================================",
    "\nFitting model: ",
    model_name,
    "\n========================================"
  )
  
  fit_arguments <- list(
    formula = formula,
    data = data,
    mesh = mesh,
    family = family,
    spatial = spatial,
    spatiotemporal = spatiotemporal,
    silent = FALSE
  )
  
  if (!is.null(time)) {
    fit_arguments$time <- time
  }
  
  if (!is.null(offset)) {
    fit_arguments$offset <- offset
  }
  
  if (!is.null(weights)) {
    fit_arguments$weights <- weights
  }
  
  fitted_model <- do.call(sdmTMB::sdmTMB, fit_arguments)
  
  saveRDS(fitted_model, output_file, compress = "xz")
  
  fitted_model
}


build_model_comparison <- function(models, response_name) {
  rbindlist(lapply(names(models), function(model_name) {
    model <- models[[model_name]]
    
    data.table(
      response = response_name,
      model = model_name,
      n = nrow(model$data),
      log_likelihood = safe_loglik(model),
      AIC = safe_aic(model)
    )
  }),
  use.names = TRUE,
  fill = TRUE)[, delta_AIC := AIC - min(AIC, na.rm = TRUE), by = response]
}


# ============================================================
# Input data
# ============================================================

if (!file.exists(metrics_file)) {
  stop("Annual extreme-rainfall metrics do not exist: ",
       metrics_file,
       call. = FALSE)
}

rain <- as.data.table(arrow::read_parquet(metrics_file))

required_columns <- c(
  "station_number",
  "station_name",
  "year",
  "latitude",
  "longitude",
  "n_valid_1day",
  "wet_days_valid",
  "mean_wet_day_mm",
  "rx1day_mm",
  "days_above_reference_p95",
  "days_above_reference_p99",
  "annual_complete_95",
  "reference_threshold_valid"
)

missing_columns <- setdiff(required_columns, names(rain))

if (length(missing_columns) > 0L) {
  stop("Annual metrics are missing columns: ",
       paste(missing_columns, collapse = ", "),
       call. = FALSE)
}


# ============================================================
# Analysis population
# ============================================================

model_data <- rain[year >= analysis_start_year &
                     year <= analysis_end_year &
                     annual_complete_95 == TRUE &
                     reference_threshold_valid == TRUE &
                     !is.na(latitude) &
                     !is.na(longitude) &
                     n_valid_1day > 0]

model_data[, year_decade :=
             (year - mean(c(analysis_start_year, analysis_end_year))) / 10]

# Consecutive integer needed for AR(1) time indexing.
model_data[, year_id := year - analysis_start_year + 1L]

model_data[, wet_day_proportion := wet_days_valid / n_valid_1day]

model_data[, log_rx1day_mm := log(rx1day_mm)]

if (any(!is.finite(model_data$log_rx1day_mm))) {
  warning("Non-finite log Rx1day values detected; affected rows will be excluded.",
          call. = FALSE)
}


# ============================================================
# Coordinate transformation
# ============================================================

station_metadata <- unique(model_data[, .(station_number, station_name, latitude, longitude)], by = "station_number")

station_sf <- sf::st_as_sf(
  station_metadata,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

station_albers <- sf::st_transform(station_sf, 3577)

coordinates_m <- sf::st_coordinates(station_albers)

station_metadata[, `:=`(x_km = coordinates_m[, "X"] / 1000, y_km = coordinates_m[, "Y"] / 1000)]

model_data <- merge(
  model_data,
  station_metadata[, .(station_number, x_km, y_km)],
  by = "station_number",
  all.x = TRUE,
  sort = FALSE
)

setorder(model_data, year, station_number)


# ============================================================
# Base spatial mesh geometry
#
# Build the triangulation once from unique station locations.
# Response-specific sdmTMB mesh objects are subsequently
# constructed from this same triangulation but using exactly
# the rows supplied to each fitted model.
# ============================================================

mesh_data <- unique(model_data[, .(station_number, x_km, y_km)], by = "station_number")

base_rainfall_mesh <- sdmTMB::make_mesh(
  mesh_data,
  xy_cols = c("x_km", "y_km"),
  cutoff = mesh_cutoff_km,
  type = "cutoff"
)

saveRDS(base_rainfall_mesh,
        file.path(
          model_output_dir,
          paste0("rainfall_base_mesh_", mesh_cutoff_km, "km.rds")
        ),
        compress = "xz")

mesh_plot <- file.path(summary_output_dir,
                       paste0("rainfall_mesh_", mesh_cutoff_km, "km.png"))

png(mesh_plot,
    width = 1800,
    height = 1500,
    res = 180)

plot(base_rainfall_mesh)

points(
  mesh_data$x_km,
  mesh_data$y_km,
  pch = 16,
  cex = 0.35,
  col = grDevices::adjustcolor("black", alpha.f = 0.45)
)

dev.off()


# ============================================================
# Response-specific datasets
# ============================================================

wet_occurrence_data <- model_data[!is.na(wet_day_proportion) &
                                    n_valid_1day > 0]

p95_count_data <- model_data[!is.na(days_above_reference_p95) &
                               n_valid_1day > 0]

p95_count_data[, log_exposure := log(n_valid_1day)]

p99_count_data <- model_data[!is.na(days_above_reference_p99) &
                               n_valid_1day > 0]

p99_count_data[, log_exposure := log(n_valid_1day)]

intensity_data <- model_data[!is.na(mean_wet_day_mm) &
                               mean_wet_day_mm > 0 &
                               wet_days_valid > 0]

rx1day_data <- model_data[!is.na(log_rx1day_mm) &
                            is.finite(log_rx1day_mm)]

# ============================================================
# Response-specific sdmTMB mesh objects
#
# The underlying triangulation is identical across responses,
# but the observation projection matrices must correspond
# exactly to the rows supplied to each model.
# ============================================================

wet_occurrence_mesh <- sdmTMB::make_mesh(wet_occurrence_data,
                                         xy_cols = c("x_km", "y_km"),
                                         mesh = base_rainfall_mesh$mesh)

p95_mesh <- sdmTMB::make_mesh(p95_count_data,
                              xy_cols = c("x_km", "y_km"),
                              mesh = base_rainfall_mesh$mesh)

p99_mesh <- sdmTMB::make_mesh(p99_count_data,
                              xy_cols = c("x_km", "y_km"),
                              mesh = base_rainfall_mesh$mesh)

intensity_mesh <- sdmTMB::make_mesh(intensity_data,
                                    xy_cols = c("x_km", "y_km"),
                                    mesh = base_rainfall_mesh$mesh)

rx1day_mesh <- sdmTMB::make_mesh(rx1day_data,
                                 xy_cols = c("x_km", "y_km"),
                                 mesh = base_rainfall_mesh$mesh)

# ============================================================
# 1. Wet-day occurrence
# ============================================================

wet_occurrence_models <- list()

wet_occurrence_models$wet_occurrence_spatial <-
  fit_or_load_sdmTMB(
    model_name = "wet_occurrence_spatial",
    formula = wet_day_proportion ~ year_decade,
    data = wet_occurrence_data,
    mesh = wet_occurrence_mesh,
    family = stats::binomial(link = "logit"),
    spatial = "on",
    spatiotemporal = "off",
    weights = wet_occurrence_data$n_valid_1day,
    overwrite = overwrite_models
  )

wet_occurrence_models$wet_occurrence_ar1 <-
  fit_or_load_sdmTMB(
    model_name = "wet_occurrence_ar1",
    formula = wet_day_proportion ~ year_decade,
    data = wet_occurrence_data,
    mesh = wet_occurrence_mesh,
    family = stats::binomial(link = "logit"),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year_id",
    weights = wet_occurrence_data$n_valid_1day,
    overwrite = overwrite_models
  )


# ============================================================
# 2. P95 exceedance frequency
# ============================================================

p95_models <- list()

p95_models$p95_frequency_spatial <-
  fit_or_load_sdmTMB(
    model_name = "p95_frequency_spatial",
    formula = days_above_reference_p95 ~ year_decade,
    data = p95_count_data,
    mesh = p95_mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = "on",
    spatiotemporal = "off",
    offset = p95_count_data$log_exposure,
    overwrite = overwrite_models
  )

p95_models$p95_frequency_ar1 <-
  fit_or_load_sdmTMB(
    model_name = "p95_frequency_ar1",
    formula = days_above_reference_p95 ~ year_decade,
    data = p95_count_data,
    mesh = p95_mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year_id",
    offset = p95_count_data$log_exposure,
    overwrite = overwrite_models
  )


# ============================================================
# 3. P99 exceedance frequency
# ============================================================

p99_models <- list()

p99_models$p99_frequency_spatial <-
  fit_or_load_sdmTMB(
    model_name = "p99_frequency_spatial",
    formula = days_above_reference_p99 ~ year_decade,
    data = p99_count_data,
    mesh = p99_mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = "on",
    spatiotemporal = "off",
    offset = p99_count_data$log_exposure,
    overwrite = overwrite_models
  )

p99_models$p99_frequency_ar1 <-
  fit_or_load_sdmTMB(
    model_name = "p99_frequency_ar1",
    formula = days_above_reference_p99 ~ year_decade,
    data = p99_count_data,
    mesh = p99_mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year_id",
    offset = p99_count_data$log_exposure,
    overwrite = overwrite_models
  )


# ============================================================
# 4. Wet-day intensity
# ============================================================

intensity_models <- list()

intensity_models$wet_intensity_spatial <-
  fit_or_load_sdmTMB(
    model_name = "wet_intensity_spatial",
    formula = mean_wet_day_mm ~ year_decade,
    data = intensity_data,
    mesh = intensity_mesh,
    family = stats::Gamma(link = "log"),
    spatial = "on",
    spatiotemporal = "off",
    weights = intensity_data$wet_days_valid,
    overwrite = overwrite_models
  )

intensity_models$wet_intensity_ar1 <-
  fit_or_load_sdmTMB(
    model_name = "wet_intensity_ar1",
    formula = mean_wet_day_mm ~ year_decade,
    data = intensity_data,
    mesh = intensity_mesh,
    family = stats::Gamma(link = "log"),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year_id",
    weights = intensity_data$wet_days_valid,
    overwrite = overwrite_models
  )


# ============================================================
# 5. Log Rx1day
# ============================================================

rx1day_models <- list()

rx1day_models$log_rx1day_spatial <-
  fit_or_load_sdmTMB(
    model_name = "log_rx1day_spatial",
    formula = log_rx1day_mm ~ year_decade,
    data = rx1day_data,
    mesh = rx1day_mesh,
    family = stats::gaussian(link = "identity"),
    spatial = "on",
    spatiotemporal = "off",
    overwrite = overwrite_models
  )

rx1day_models$log_rx1day_ar1 <-
  fit_or_load_sdmTMB(
    model_name = "log_rx1day_ar1",
    formula = log_rx1day_mm ~ year_decade,
    data = rx1day_data,
    mesh = rx1day_mesh,
    family = stats::gaussian(link = "identity"),
    spatial = "on",
    spatiotemporal = "ar1",
    time = "year_id",
    overwrite = overwrite_models
  )


# ============================================================
# Model comparisons
# ============================================================

model_comparison <- rbindlist(
  list(
    build_model_comparison(wet_occurrence_models, "Wet-day occurrence"),
    build_model_comparison(p95_models, "P95 exceedance frequency"),
    build_model_comparison(p99_models, "P99 exceedance frequency"),
    build_model_comparison(intensity_models, "Mean wet-day intensity"),
    build_model_comparison(rx1day_models, "Log Rx1day")
  ),
  use.names = TRUE,
  fill = TRUE
)

setorder(model_comparison, response, AIC)


# ============================================================
# Fixed-effect summaries
# ============================================================

all_models <- c(wet_occurrence_models,
                p95_models,
                p99_models,
                intensity_models,
                rx1day_models)

response_lookup <- c(
  wet_occurrence_spatial = "Wet-day occurrence",
  wet_occurrence_ar1 = "Wet-day occurrence",
  p95_frequency_spatial = "P95 exceedance frequency",
  p95_frequency_ar1 = "P95 exceedance frequency",
  p99_frequency_spatial = "P99 exceedance frequency",
  p99_frequency_ar1 = "P99 exceedance frequency",
  wet_intensity_spatial = "Mean wet-day intensity",
  wet_intensity_ar1 = "Mean wet-day intensity",
  log_rx1day_spatial = "Log Rx1day",
  log_rx1day_ar1 = "Log Rx1day"
)

fixed_effects <- rbindlist(lapply(names(all_models), function(model_name) {
  extract_fixed_effects(model = all_models[[model_name]],
                        model_name = model_name,
                        response_name = response_lookup[[model_name]])
}),
use.names = TRUE,
fill = TRUE)


# ============================================================
# Convert time effects to interpretable multiplicative changes
# ============================================================

trend_effects <- fixed_effects[term == "year_decade"]

trend_effects[, `:=`(
  multiplicative_change_per_decade = exp(estimate),
  multiplicative_conf_low = exp(conf_low),
  multiplicative_conf_high = exp(conf_high),
  percent_change_per_decade = 100 * (exp(estimate) - 1),
  percent_change_conf_low = 100 * (exp(conf_low) - 1),
  percent_change_conf_high = 100 * (exp(conf_high) - 1)
)]

# For Gaussian log-Rx1day, exponentiating the coefficient also gives
# the multiplicative change in Rx1day per decade.
# For wet-day occurrence, this is an odds ratio rather than a direct
# percentage-point change in probability.


# ============================================================
# Model-data audit
# ============================================================

model_data_audit <- data.table(
  model = c(
    "wet_occurrence",
    "p95_frequency",
    "p99_frequency",
    "wet_intensity",
    "log_rx1day"
  ),
  rows = c(
    nrow(wet_occurrence_data),
    nrow(p95_count_data),
    nrow(p99_count_data),
    nrow(intensity_data),
    nrow(rx1day_data)
  ),
  stations = c(
    uniqueN(wet_occurrence_data$station_number),
    uniqueN(p95_count_data$station_number),
    uniqueN(p99_count_data$station_number),
    uniqueN(intensity_data$station_number),
    uniqueN(rx1day_data$station_number)
  ),
  first_year = c(
    min(wet_occurrence_data$year),
    min(p95_count_data$year),
    min(p99_count_data$year),
    min(intensity_data$year),
    min(rx1day_data$year)
  ),
  last_year = c(
    max(wet_occurrence_data$year),
    max(p95_count_data$year),
    max(p99_count_data$year),
    max(intensity_data$year),
    max(rx1day_data$year)
  )
)


# ============================================================
# Save outputs
# ============================================================

fwrite(model_data_audit,
       file.path(summary_output_dir, "model_data_audit.csv"))

fwrite(model_comparison,
       file.path(summary_output_dir, "model_comparison.csv"))

fwrite(fixed_effects,
       file.path(summary_output_dir, "fixed_effects.csv"))

fwrite(trend_effects,
       file.path(summary_output_dir, "trend_effects_per_decade.csv"))

saveRDS(
  list(
    wet_occurrence = wet_occurrence_models,
    p95_frequency = p95_models,
    p99_frequency = p99_models,
    wet_intensity = intensity_models,
    log_rx1day = rx1day_models
  ),
  file.path(
    model_output_dir,
    "rainfall_spatiotemporal_model_bundle.rds"
  ),
  compress = "xz"
)


# ============================================================
# Final output
# ============================================================

message(
  "\n========================================",
  "\nSPATIOTEMPORAL MODELS COMPLETE",
  "\n========================================"
)

print(model_data_audit)
print(model_comparison)
print(trend_effects)
