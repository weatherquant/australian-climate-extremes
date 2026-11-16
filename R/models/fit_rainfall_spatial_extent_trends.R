# ============================================================
# fit_rainfall_spatial_extent_trends.R
#
# Model changes in the size and duration of coherent Australian
# P97.5 rainfall episodes, using the validated 200-km clusters.
#
# Primary inferential period: 1950-2022
# Sensitivity period:          1950-2023
# 2024-2025:                  descriptive only
#
# Primary responses:
#   1. maximum_simultaneous_stations
#   2. unique_stations_affected
#   3. duration_days
#
# Secondary response:
#   4. maximum_cluster_area_km2
#
# Candidate temporal structures:
#   M0 stationary
#   M1 linear trend per decade
#   M2 three periods: 1950-1975, 1976-2000, 2001-2022
# ============================================================

library(data.table)
library(arrow)
library(MASS)

# ============================================================
# Configuration
# ============================================================

episode_file <- file.path(
  "data-derived",
  "rainfall",
  "spatial_extremes",
  "rainfall_spatial_extreme_episodes_p975_200km_1950_2025.parquet"
)

output_dir <- file.path("data-derived", "rainfall", "spatial_extremes")

table_output_dir <- file.path("tables", "rainfall_spatial_extremes")

model_output_file <- file.path(output_dir, "rainfall_spatial_extent_models.rds")

primary_start_year <- 1950L
primary_end_year <- 2022L
sensitivity_end_year <- 2023L
centre_year <- mean(c(primary_start_year, primary_end_year))

minimum_episode_stations <- 3L
minimum_duration_days <- 1L

# ============================================================
# Input checks and preparation
# ============================================================

if (!file.exists(episode_file)) {
  stop("Missing episode file: ", episode_file, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

episodes <- as.data.table(read_parquet(episode_file))

required_columns <- c(
  "episode_id",
  "start_date",
  "start_year",
  "all_days_primary_coverage",
  "all_days_strict_coverage",
  "maximum_simultaneous_stations",
  "unique_stations_affected",
  "duration_days",
  "maximum_cluster_area_km2"
)

missing_columns <- setdiff(required_columns, names(episodes))
if (length(missing_columns) > 0L) {
  stop("Episode file is missing: ",
       paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

# ============================================================
# Duration response
#
# Every episode lasts at least one day. Model additional days
# beyond the mandatory first day:
#
#   1-day episode -> 0
#   2-day episode -> 1
#   3-day episode -> 2
# ============================================================

episodes[, duration_extra_days :=
           as.integer(duration_days - 1L)]

if (any(is.na(episodes$duration_extra_days)) ||
    any(episodes$duration_extra_days < 0L)) {
  stop("Invalid duration_extra_days values.", call. = FALSE)
}

episodes[, year_decade := (start_year - centre_year) / 10]

episodes[, spatial_climate_period := fcase(
  start_year >= 1950L & start_year <= 1975L,
  "1950-1975",
  start_year >= 1976L & start_year <= 2000L,
  "1976-2000",
  start_year >= 2001L & start_year <= 2022L,
  "2001-2022",
  default = NA_character_
)]

episodes[, spatial_climate_period := factor(spatial_climate_period,
                                            levels = c("1950-1975", "1976-2000", "2001-2022"))]

episodes[, primary_analysis_period :=
           start_year >= primary_start_year &
           start_year <= primary_end_year]

episodes[, sensitivity_analysis_period :=
           start_year >= primary_start_year &
           start_year <= sensitivity_end_year]

analysis_data <- episodes[primary_analysis_period == TRUE &
                            all_days_primary_coverage == TRUE &
                            maximum_simultaneous_stations >= minimum_episode_stations &
                            is.finite(unique_stations_affected) &
                            is.finite(maximum_cluster_area_km2) &
                            is.finite(duration_extra_days) &
                            maximum_cluster_area_km2 > 0]

if (nrow(analysis_data) == 0L) {
  stop("No eligible primary-analysis episodes.", call. = FALSE)
}

# ============================================================
# Model helpers
# ============================================================

fit_nb_triplet <- function(data, response) {
  formula_m0 <- as.formula(paste(response, "~ 1"))
  formula_m1 <- as.formula(paste(response, "~ year_decade"))
  formula_m2 <- as.formula(paste(response, "~ spatial_climate_period"))
  
  list(
    M0 = MASS::glm.nb(formula_m0, data = data, control = glm.control(maxit = 100)),
    M1 = MASS::glm.nb(formula_m1, data = data, control = glm.control(maxit = 100)),
    M2 = MASS::glm.nb(formula_m2, data = data, control = glm.control(maxit = 100))
  )
}

fit_poisson_triplet <- function(data, response) {
  formula_m0 <-
    as.formula(paste(response, "~ 1"))
  
  formula_m1 <-
    as.formula(paste(response, "~ year_decade"))
  
  formula_m2 <-
    as.formula(paste(response, "~ spatial_climate_period"))
  
  list(
    M0 =
      glm(formula_m0, family = poisson(link = "log"), data = data),
    
    M1 =
      glm(formula_m1, family = poisson(link = "log"), data = data),
    
    M2 =
      glm(formula_m2, family = poisson(link = "log"), data = data)
  )
}

fit_lognormal_triplet <- function(data, response) {
  transformed <- paste0("log(", response, ")")
  
  list(
    M0 = lm(as.formula(paste(transformed, "~ 1")), data = data),
    M1 = lm(as.formula(paste(
      transformed, "~ year_decade"
    )), data = data),
    M2 = lm(as.formula(
      paste(transformed, "~ spatial_climate_period")
    ), data = data)
  )
}

model_comparison <- function(fits, response, family_name) {
  out <- rbindlist(lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    data.table(
      response = response,
      family = family_name,
      model = model_name,
      parameters = attr(logLik(fit), "df"),
      log_likelihood = as.numeric(logLik(fit)),
      AIC = AIC(fit)
    )
  }))
  
  out[, delta_AIC := AIC - min(AIC)]
  setorder(out, AIC)
  out
}

extract_linear_effect <- function(fit,
                                  response,
                                  family_name,
                                  analysis_period) {
  beta <- coef(fit)["year_decade"]
  se <- sqrt(vcov(fit)["year_decade", "year_decade"])
  z <- qnorm(0.975)
  p <- summary(fit)$coefficients["year_decade", ncol(summary(fit)$coefficients)]
  
  data.table(
    response = response,
    family = family_name,
    analysis_period = analysis_period,
    log_change_per_decade = as.numeric(beta),
    standard_error = as.numeric(se),
    ratio_per_decade = as.numeric(exp(beta)),
    percent_change_per_decade = as.numeric(100 * (exp(beta) - 1)),
    conf_low_percent = as.numeric(100 * (exp(beta - z * se) - 1)),
    conf_high_percent = as.numeric(100 * (exp(beta + z * se) - 1)),
    p_value = as.numeric(p)
  )
}

extract_period_effects <- function(fit,
                                   response,
                                   family_name,
                                   analysis_period) {
  terms <- c("spatial_climate_period1976-2000",
             "spatial_climate_period2001-2022")
  
  beta <- coef(fit)[terms]
  se <- sqrt(diag(vcov(fit)))[terms]
  z <- qnorm(0.975)
  p <- summary(fit)$coefficients[terms, ncol(summary(fit)$coefficients)]
  
  data.table(
    response = response,
    family = family_name,
    analysis_period = analysis_period,
    comparison = c("1976-2000 vs 1950-1975", "2001-2022 vs 1950-1975"),
    log_difference = as.numeric(beta),
    standard_error = as.numeric(se),
    ratio = as.numeric(exp(beta)),
    percent_difference = as.numeric(100 * (exp(beta) - 1)),
    conf_low_percent = as.numeric(100 * (exp(beta - z * se) - 1)),
    conf_high_percent = as.numeric(100 * (exp(beta + z * se) - 1)),
    p_value = as.numeric(p)
  )
}

# ============================================================
# Primary models, 1950-2022
# ============================================================

response_configuration <- list(
  maximum_simultaneous_stations = "negative_binomial",
  unique_stations_affected = "negative_binomial",
  duration_extra_days = "poisson",
  maximum_cluster_area_km2 = "lognormal"
)

fits_primary <- list()
comparison_list <- list()
linear_list <- list()
period_list <- list()

for (response in names(response_configuration)) {
  family_name <- response_configuration[[response]]
  
  message("Fitting spatial-extent models for: ", response)
  
  fits <-
    if (family_name ==
        "negative_binomial") {
      fit_nb_triplet(analysis_data, response)
      
    } else if (family_name ==
               "poisson") {
      fit_poisson_triplet(analysis_data, response)
      
    } else {
      fit_lognormal_triplet(analysis_data, response)
    }
  
  fits_primary[[response]] <- fits
  
  comparison_list[[response]] <- model_comparison(fits, response, family_name)
  
  linear_list[[response]] <- extract_linear_effect(fits$M1, response, family_name, "1950-2022")
  
  period_list[[response]] <- extract_period_effects(fits$M2, response, family_name, "1950-2022")
}

model_comparisons <- rbindlist(comparison_list)
linear_effects <- rbindlist(linear_list)
period_effects <- rbindlist(period_list)

# ============================================================
# 1950-2023 linear sensitivity
# ============================================================

sensitivity_data <- episodes[sensitivity_analysis_period == TRUE &
                               all_days_primary_coverage == TRUE &
                               maximum_simultaneous_stations >= minimum_episode_stations &
                               is.finite(unique_stations_affected) &
                               is.finite(maximum_cluster_area_km2) &
                               is.finite(duration_extra_days) &
                               maximum_cluster_area_km2 > 0]

sensitivity_list <- list()

for (response in names(response_configuration)) {
  family_name <- response_configuration[[response]]
  
  fit <-
    if (family_name ==
        "negative_binomial") {
      MASS::glm.nb(as.formula(paste(response, "~ year_decade")),
                   data =
                     sensitivity_data,
                   control =
                     glm.control(maxit = 100))
      
    } else if (family_name ==
               "poisson") {
      glm(as.formula(paste(response, "~ year_decade")),
          family =
            poisson(link = "log"),
          data =
            sensitivity_data)
      
    } else {
      lm(as.formula(paste0("log(", response, ") ~ year_decade")), data =
           sensitivity_data)
    }
  
  sensitivity_list[[response]] <- extract_linear_effect(fit, response, family_name, "1950-2023_sensitivity")
}

linear_sensitivity_2023 <- rbindlist(sensitivity_list)

# ============================================================
# Descriptive period summaries
# ============================================================

period_descriptives <- analysis_data[, .(
  episodes = .N,
  median_maximum_simultaneous_stations = median(maximum_simultaneous_stations),
  q10_maximum_simultaneous_stations = quantile(maximum_simultaneous_stations, 0.10, names = FALSE),
  q90_maximum_simultaneous_stations = quantile(maximum_simultaneous_stations, 0.90, names = FALSE),
  median_unique_stations_affected = median(unique_stations_affected),
  q10_unique_stations_affected = quantile(unique_stations_affected, 0.10, names = FALSE),
  q90_unique_stations_affected = quantile(unique_stations_affected, 0.90, names = FALSE),
  median_duration_days = median(duration_days),
  q90_duration_days = quantile(duration_days, 0.90, names = FALSE),
  median_maximum_cluster_area_km2 = median(maximum_cluster_area_km2),
  q90_maximum_cluster_area_km2 = quantile(maximum_cluster_area_km2, 0.90, names = FALSE)
), by = spatial_climate_period]

# ============================================================
# Save outputs
# ============================================================

fwrite(
  model_comparisons,
  file.path(
    table_output_dir,
    "rainfall_spatial_extent_model_comparison.csv"
  )
)

fwrite(
  linear_effects,
  file.path(
    table_output_dir,
    "rainfall_spatial_extent_linear_effects.csv"
  )
)

fwrite(
  period_effects,
  file.path(
    table_output_dir,
    "rainfall_spatial_extent_period_effects.csv"
  )
)

fwrite(
  linear_sensitivity_2023,
  file.path(
    table_output_dir,
    "rainfall_spatial_extent_linear_sensitivity_2023.csv"
  )
)

fwrite(
  period_descriptives,
  file.path(
    table_output_dir,
    "rainfall_spatial_extent_period_descriptives.csv"
  )
)

saveRDS(
  list(
    configuration = list(
      primary_analysis_period = "1950-2022",
      sensitivity_analysis_period = "1950-2023",
      threshold = "p975",
      connectivity_distance_km = 200,
      minimum_coherent_stations = minimum_episode_stations,
      primary_coverage = "80/80"
    ),
    primary_fits = fits_primary,
    model_comparisons = model_comparisons,
    linear_effects = linear_effects,
    period_effects = period_effects,
    linear_sensitivity_2023 = linear_sensitivity_2023,
    period_descriptives = period_descriptives
  ),
  model_output_file,
  compress = "xz"
)

# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nRAINFALL SPATIAL EXTENT MODELS COMPLETE",
  "\n========================================"
)

message("\nLinear effects, 1950-2022:")
print(linear_effects)

message("\nPeriod effects, 1950-2022:")
print(period_effects)

message("\nModel comparisons:")
print(model_comparisons)

message("\n1950-2023 sensitivity:")
print(linear_sensitivity_2023)

message("\nPeriod descriptives:")
print(period_descriptives)
