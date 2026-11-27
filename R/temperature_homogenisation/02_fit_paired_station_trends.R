# ============================================================
# 02_fit_paired_station_trends.R
#
# Fit paired raw and homogenised annual temperature trends for
# ACORN-SAT v2.3 locations.
#
# Key design:
#   - raw and homogenised trends use exactly the same annual values
#   - trend_difference = homogenised trend - raw trend
#   - the trend difference is fitted directly from the annual
#     adjustment series, giving the correct paired uncertainty
#   - OLS slopes are primary
#   - Newey-West HAC standard errors are used for inference
#   - Theil-Sen slopes are retained as a robust sensitivity
#
# Periods:
#   1910-2021
#   1950-2021
#   1970-2021
#   1990-2021
#
# Eligibility for each station x element x period:
#   - first paired year <= period start + 5 years
#   - last paired year  >= period end - 5 years
#   - >= 80% of requested years have paired annual values
#   - >= 20 paired annual values
# ============================================================

# ============================================================
# Configuration
# ============================================================

annual_input_file <- file.path(
  "data-derived",
  "temperature_homogenisation",
  "temperature_raw_homogenised_paired_annual_v2.3.parquet"
)

output_dir <- file.path(
  "data-derived",
  "temperature_homogenisation"
)

table_output_dir <- file.path(
  "tables",
  "temperature_homogenisation"
)

period_definitions <- data.table::data.table(
  period_id = c(
    "1910_2021",
    "1950_2021",
    "1970_2021",
    "1990_2021"
  ),
  period_start = c(1910L, 1950L, 1970L, 1990L),
  period_end = c(2021L, 2021L, 2021L, 2021L)
)

minimum_coverage_fraction <- 0.80
minimum_valid_years <- 20L
endpoint_tolerance_years <- 5L
negligible_trend_effect <- 0.01
confidence_level <- 0.95

# Newey-West automatic lag when NULL. This uses the standard
# floor(4 * (n / 100)^(2/9)) rule and caps the lag at n - 2.
hac_lag_override <- NULL

# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "data.table",
  "arrow",
  "sandwich"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(annual_input_file)) {
  stop(
    "Validated paired annual input file not found: ",
    annual_input_file,
    call. = FALSE
  )
}

# ============================================================
# Helper functions
# ============================================================

automatic_hac_lag <- function(n) {
  if (!is.null(hac_lag_override)) {
    return(
      as.integer(
        min(
          max(hac_lag_override, 0L),
          max(n - 2L, 0L)
        )
      )
    )
  }
  
  as.integer(
    min(
      floor(4 * (n / 100)^(2 / 9)),
      max(n - 2L, 0L)
    )
  )
}

safe_theil_sen <- function(year, value) {
  valid <- is.finite(year) & is.finite(value)
  year <- as.numeric(year[valid])
  value <- as.numeric(value[valid])
  
  n <- length(year)
  
  if (n < 2L) {
    return(NA_real_)
  }
  
  # Annual series contain at most 112 observations, so the complete
  # set of pairwise slopes is small and can be evaluated directly.
  slopes <- vector(
    mode = "numeric",
    length = n * (n - 1L) / 2L
  )
  
  k <- 0L
  
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    count <- length(j)
    
    slopes[(k + 1L):(k + count)] <-
      (value[j] - value[i]) /
      (year[j] - year[i])
    
    k <- k + count
  }
  
  as.numeric(median(slopes, na.rm = TRUE) * 10)
}

fit_one_trend <- function(year, value, confidence_level = 0.95) {
  valid <- is.finite(year) & is.finite(value)
  
  dat <- data.frame(
    year = as.numeric(year[valid]),
    value = as.numeric(value[valid])
  )
  
  n <- nrow(dat)
  
  if (n < 3L || length(unique(dat$year)) < 2L) {
    return(list(
      n = as.integer(n),
      slope_c_per_decade = NA_real_,
      hac_se_c_per_decade = NA_real_,
      conf_low_c_per_decade = NA_real_,
      conf_high_c_per_decade = NA_real_,
      p_value = NA_real_,
      theil_sen_c_per_decade = NA_real_,
      hac_lag = NA_integer_,
      intercept = NA_real_,
      residual_sd = NA_real_,
      r_squared = NA_real_,
      model = NULL
    ))
  }
  
  fit <- stats::lm(
    value ~ year,
    data = dat
  )
  
  hac_lag <- automatic_hac_lag(n)
  
  covariance_hac <- sandwich::NeweyWest(
    fit,
    lag = hac_lag,
    prewhite = FALSE,
    adjust = TRUE
  )
  
  beta <- stats::coef(fit)[["year"]]
  se <- sqrt(covariance_hac["year", "year"])
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  
  statistic <- beta / se
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  
  list(
    n = as.integer(n),
    slope_c_per_decade = as.numeric(beta * 10),
    hac_se_c_per_decade = as.numeric(se * 10),
    conf_low_c_per_decade = as.numeric((beta - z * se) * 10),
    conf_high_c_per_decade = as.numeric((beta + z * se) * 10),
    p_value = as.numeric(p_value),
    theil_sen_c_per_decade = safe_theil_sen(dat$year, dat$value),
    hac_lag = as.integer(hac_lag),
    intercept = as.numeric(stats::coef(fit)[["(Intercept)"]]),
    residual_sd = as.numeric(stats::sigma(fit)),
    r_squared = as.numeric(summary(fit)$r.squared),
    model = fit
  )
}

classify_trend_effect <- function(trend_difference) {
  data.table::fcase(
    !is.finite(trend_difference),
    NA_character_,
    
    trend_difference > negligible_trend_effect,
    "raw_understates_warming",
    
    trend_difference < -negligible_trend_effect,
    "raw_overstates_warming",
    
    default = "negligible_effect"
  )
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  as.numeric(
    stats::quantile(
      x,
      probability,
      names = FALSE,
      type = 7
    )
  )
}

# ============================================================
# Load and validate annual paired data
# ============================================================

annual <- data.table::as.data.table(
  arrow::read_parquet(annual_input_file)
)

required_columns <- c(
  "acorn_station_number",
  "station_name",
  "year",
  "tmax_annual_raw",
  "tmax_annual_hom",
  "tmin_annual_raw",
  "tmin_annual_hom",
  "tmean_annual_raw",
  "tmean_annual_hom"
)

missing_columns <- setdiff(required_columns, names(annual))

if (length(missing_columns) > 0L) {
  stop(
    "Paired annual product is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (annual[, anyDuplicated(paste(acorn_station_number, year))] > 0L) {
  stop(
    "Duplicated ACORN location-year rows in the annual paired product.",
    call. = FALSE
  )
}

setorder(annual, acorn_station_number, year)

# Convert wide paired annual values to a compact long table.
annual_long <- rbindlist(
  list(
    annual[, .(
      acorn_station_number,
      station_name,
      year,
      element = "tmax",
      raw_value = tmax_annual_raw,
      hom_value = tmax_annual_hom
    )],
    
    annual[, .(
      acorn_station_number,
      station_name,
      year,
      element = "tmin",
      raw_value = tmin_annual_raw,
      hom_value = tmin_annual_hom
    )],
    
    annual[, .(
      acorn_station_number,
      station_name,
      year,
      element = "tmean",
      raw_value = tmean_annual_raw,
      hom_value = tmean_annual_hom
    )]
  ),
  use.names = TRUE
)

annual_long[, paired_value_available :=
              is.finite(raw_value) & is.finite(hom_value)
]

annual_long[, adjustment_value :=
              data.table::fifelse(
                paired_value_available,
                hom_value - raw_value,
                NA_real_
              )
]

setorder(annual_long, acorn_station_number, element, year)

# ============================================================
# Fit paired station trends
# ============================================================

trend_results <- vector(
  mode = "list",
  length = nrow(period_definitions) *
    uniqueN(annual_long$acorn_station_number) *
    uniqueN(annual_long$element)
)

model_store <- list()
result_index <- 0L

for (period_row in seq_len(nrow(period_definitions))) {
  period_id <- period_definitions$period_id[period_row]
  period_start <- period_definitions$period_start[period_row]
  period_end <- period_definitions$period_end[period_row]
  requested_years <- period_end - period_start + 1L
  
  period_data <- annual_long[
    year >= period_start &
      year <= period_end
  ]
  
  group_keys <- unique(
    period_data[, .(
      acorn_station_number,
      station_name,
      element
    )]
  )
  
  for (group_row in seq_len(nrow(group_keys))) {
    station_number <- group_keys$acorn_station_number[group_row]
    station_name <- group_keys$station_name[group_row]
    element <- group_keys$element[group_row]
    
    dat <- period_data[
      acorn_station_number == station_number &
        element == group_keys$element[group_row]
    ]
    
    paired <- dat[paired_value_available == TRUE]
    
    valid_years <- nrow(paired)
    coverage_fraction <- valid_years / requested_years
    
    first_valid_year <- if (valid_years > 0L) {
      min(paired$year)
    } else {
      NA_integer_
    }
    
    last_valid_year <- if (valid_years > 0L) {
      max(paired$year)
    } else {
      NA_integer_
    }
    
    start_covered <-
      is.finite(first_valid_year) &&
      first_valid_year <= period_start + endpoint_tolerance_years
    
    end_covered <-
      is.finite(last_valid_year) &&
      last_valid_year >= period_end - endpoint_tolerance_years
    
    eligible <-
      valid_years >= minimum_valid_years &&
      coverage_fraction >= minimum_coverage_fraction &&
      start_covered &&
      end_covered
    
    result_index <- result_index + 1L
    
    base_result <- data.table(
      acorn_station_number = station_number,
      station_name = station_name,
      element = element,
      period_id = period_id,
      period_start = period_start,
      period_end = period_end,
      requested_years = requested_years,
      valid_years = as.integer(valid_years),
      coverage_fraction = as.numeric(coverage_fraction),
      first_valid_year = as.integer(first_valid_year),
      last_valid_year = as.integer(last_valid_year),
      start_covered = start_covered,
      end_covered = end_covered,
      eligible = eligible
    )
    
    if (!eligible) {
      trend_results[[result_index]] <- base_result[, `:=`(
        raw_trend_c_per_decade = NA_real_,
        raw_hac_se_c_per_decade = NA_real_,
        raw_conf_low_c_per_decade = NA_real_,
        raw_conf_high_c_per_decade = NA_real_,
        raw_p_value = NA_real_,
        raw_theil_sen_c_per_decade = NA_real_,
        raw_hac_lag = NA_integer_,
        raw_r_squared = NA_real_,
        
        hom_trend_c_per_decade = NA_real_,
        hom_hac_se_c_per_decade = NA_real_,
        hom_conf_low_c_per_decade = NA_real_,
        hom_conf_high_c_per_decade = NA_real_,
        hom_p_value = NA_real_,
        hom_theil_sen_c_per_decade = NA_real_,
        hom_hac_lag = NA_integer_,
        hom_r_squared = NA_real_,
        
        trend_difference_c_per_decade = NA_real_,
        difference_hac_se_c_per_decade = NA_real_,
        difference_conf_low_c_per_decade = NA_real_,
        difference_conf_high_c_per_decade = NA_real_,
        difference_p_value = NA_real_,
        difference_theil_sen_c_per_decade = NA_real_,
        difference_hac_lag = NA_integer_,
        difference_r_squared = NA_real_,
        trend_effect_direction = NA_character_,
        ols_difference_identity_error = NA_real_
      )]
      
      next
    }
    
    raw_fit <- fit_one_trend(
      paired$year,
      paired$raw_value,
      confidence_level
    )
    
    hom_fit <- fit_one_trend(
      paired$year,
      paired$hom_value,
      confidence_level
    )
    
    difference_fit <- fit_one_trend(
      paired$year,
      paired$adjustment_value,
      confidence_level
    )
    
    trend_difference <-
      hom_fit$slope_c_per_decade -
      raw_fit$slope_c_per_decade
    
    identity_error <-
      trend_difference -
      difference_fit$slope_c_per_decade
    
    model_key <- paste(
      station_number,
      element,
      period_id,
      sep = "__"
    )
    
    model_store[[model_key]] <- list(
      raw = raw_fit$model,
      homogenised = hom_fit$model,
      difference = difference_fit$model,
      valid_years = paired$year,
      hac_lag = difference_fit$hac_lag
    )
    
    trend_results[[result_index]] <- base_result[, `:=`(
      raw_trend_c_per_decade = raw_fit$slope_c_per_decade,
      raw_hac_se_c_per_decade = raw_fit$hac_se_c_per_decade,
      raw_conf_low_c_per_decade = raw_fit$conf_low_c_per_decade,
      raw_conf_high_c_per_decade = raw_fit$conf_high_c_per_decade,
      raw_p_value = raw_fit$p_value,
      raw_theil_sen_c_per_decade = raw_fit$theil_sen_c_per_decade,
      raw_hac_lag = raw_fit$hac_lag,
      raw_r_squared = raw_fit$r_squared,
      
      hom_trend_c_per_decade = hom_fit$slope_c_per_decade,
      hom_hac_se_c_per_decade = hom_fit$hac_se_c_per_decade,
      hom_conf_low_c_per_decade = hom_fit$conf_low_c_per_decade,
      hom_conf_high_c_per_decade = hom_fit$conf_high_c_per_decade,
      hom_p_value = hom_fit$p_value,
      hom_theil_sen_c_per_decade = hom_fit$theil_sen_c_per_decade,
      hom_hac_lag = hom_fit$hac_lag,
      hom_r_squared = hom_fit$r_squared,
      
      trend_difference_c_per_decade =
        difference_fit$slope_c_per_decade,
      difference_hac_se_c_per_decade =
        difference_fit$hac_se_c_per_decade,
      difference_conf_low_c_per_decade =
        difference_fit$conf_low_c_per_decade,
      difference_conf_high_c_per_decade =
        difference_fit$conf_high_c_per_decade,
      difference_p_value = difference_fit$p_value,
      difference_theil_sen_c_per_decade =
        difference_fit$theil_sen_c_per_decade,
      difference_hac_lag = difference_fit$hac_lag,
      difference_r_squared = difference_fit$r_squared,
      trend_effect_direction =
        classify_trend_effect(
          difference_fit$slope_c_per_decade
        ),
      ols_difference_identity_error = as.numeric(identity_error)
    )]
  }
}

paired_station_trends <- rbindlist(
  trend_results[seq_len(result_index)],
  use.names = TRUE,
  fill = TRUE
)

setorder(
  paired_station_trends,
  period_start,
  element,
  acorn_station_number
)

maximum_identity_error <- paired_station_trends[
  eligible == TRUE,
  max(abs(ols_difference_identity_error), na.rm = TRUE)
]

if (
  !is.finite(maximum_identity_error) ||
  maximum_identity_error > 1e-10
) {
  stop(
    "Paired trend identity failed. Maximum error = ",
    maximum_identity_error,
    call. = FALSE
  )
}

# ============================================================
# Summaries by period and element
# ============================================================

period_summary <- paired_station_trends[
  eligible == TRUE,
  .(
    eligible_stations = as.integer(.N),
    
    median_raw_trend_c_per_decade =
      as.numeric(median(raw_trend_c_per_decade, na.rm = TRUE)),
    
    median_hom_trend_c_per_decade =
      as.numeric(median(hom_trend_c_per_decade, na.rm = TRUE)),
    
    median_trend_difference_c_per_decade =
      as.numeric(median(trend_difference_c_per_decade, na.rm = TRUE)),
    
    mean_trend_difference_c_per_decade =
      as.numeric(mean(trend_difference_c_per_decade, na.rm = TRUE)),
    
    q10_trend_difference_c_per_decade =
      safe_quantile(trend_difference_c_per_decade, 0.10),
    
    q25_trend_difference_c_per_decade =
      safe_quantile(trend_difference_c_per_decade, 0.25),
    
    q75_trend_difference_c_per_decade =
      safe_quantile(trend_difference_c_per_decade, 0.75),
    
    q90_trend_difference_c_per_decade =
      safe_quantile(trend_difference_c_per_decade, 0.90),
    
    stations_raw_understates_warming =
      as.integer(sum(
        trend_effect_direction == "raw_understates_warming"
      )),
    
    stations_raw_overstates_warming =
      as.integer(sum(
        trend_effect_direction == "raw_overstates_warming"
      )),
    
    stations_negligible_effect =
      as.integer(sum(
        trend_effect_direction == "negligible_effect"
      )),
    
    proportion_raw_understates_warming =
      as.numeric(mean(
        trend_effect_direction == "raw_understates_warming"
      )),
    
    proportion_raw_overstates_warming =
      as.numeric(mean(
        trend_effect_direction == "raw_overstates_warming"
      )),
    
    proportion_negligible_effect =
      as.numeric(mean(
        trend_effect_direction == "negligible_effect"
      )),
    
    stations_difference_ci_excludes_zero =
      as.integer(sum(
        difference_conf_low_c_per_decade > 0 |
          difference_conf_high_c_per_decade < 0,
        na.rm = TRUE
      )),
    
    correlation_raw_hom_trends =
      as.numeric(cor(
        raw_trend_c_per_decade,
        hom_trend_c_per_decade,
        use = "complete.obs"
      )),
    
    median_absolute_trend_effect =
      as.numeric(median(
        abs(trend_difference_c_per_decade),
        na.rm = TRUE
      )),
    
    maximum_absolute_trend_effect =
      as.numeric(max(
        abs(trend_difference_c_per_decade),
        na.rm = TRUE
      ))
  ),
  by = .(
    period_id,
    period_start,
    period_end,
    element
  )
]

setorder(period_summary, period_start, element)

direction_summary <- paired_station_trends[
  eligible == TRUE,
  .(
    stations = as.integer(.N),
    median_trend_difference_c_per_decade =
      as.numeric(median(
        trend_difference_c_per_decade,
        na.rm = TRUE
      )),
    median_absolute_trend_effect =
      as.numeric(median(
        abs(trend_difference_c_per_decade),
        na.rm = TRUE
      ))
  ),
  by = .(
    period_id,
    period_start,
    period_end,
    element,
    trend_effect_direction
  )
]

setorder(
  direction_summary,
  period_start,
  element,
  trend_effect_direction
)

# Largest positive and negative trend effects for audit/case studies.
largest_trend_effects <- rbindlist(
  lapply(
    unique(paired_station_trends$period_id),
    function(selected_period) {
      rbindlist(
        lapply(
          c("tmax", "tmin", "tmean"),
          function(selected_element) {
            x <- paired_station_trends[
              period_id == selected_period &
                element == selected_element &
                eligible == TRUE
            ]
            
            positive <- x[
              order(-trend_difference_c_per_decade)
            ][1:min(10L, .N)]
            
            negative <- x[
              order(trend_difference_c_per_decade)
            ][1:min(10L, .N)]
            
            positive[, effect_rank := "largest_positive"]
            negative[, effect_rank := "largest_negative"]
            
            rbindlist(
              list(positive, negative),
              use.names = TRUE,
              fill = TRUE
            )
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

# ============================================================
# Save outputs
# ============================================================

fwrite(
  paired_station_trends,
  file.path(
    table_output_dir,
    "temperature_paired_station_trends.csv"
  )
)

fwrite(
  period_summary,
  file.path(
    table_output_dir,
    "temperature_paired_trend_period_summary.csv"
  )
)

fwrite(
  direction_summary,
  file.path(
    table_output_dir,
    "temperature_paired_trend_direction_summary.csv"
  )
)

fwrite(
  largest_trend_effects,
  file.path(
    table_output_dir,
    "temperature_paired_largest_trend_effects.csv"
  )
)

arrow::write_parquet(
  paired_station_trends,
  file.path(
    output_dir,
    "temperature_paired_station_trends.parquet"
  ),
  compression = "zstd"
)

saveRDS(
  list(
    models = model_store,
    configuration = list(
      periods = period_definitions,
      minimum_coverage_fraction = minimum_coverage_fraction,
      minimum_valid_years = minimum_valid_years,
      endpoint_tolerance_years = endpoint_tolerance_years,
      negligible_trend_effect = negligible_trend_effect,
      confidence_level = confidence_level,
      hac_lag_override = hac_lag_override
    ),
    maximum_ols_difference_identity_error = maximum_identity_error
  ),
  file.path(
    output_dir,
    "temperature_paired_station_trend_models.rds"
  ),
  compress = "xz"
)

# ============================================================
# Console summary
# ============================================================

message("\n========================================")
message("PAIRED TEMPERATURE STATION TRENDS COMPLETE")
message("========================================\n")

message("Period and element summary:")
print(period_summary)

message("\nTrend-effect direction summary:")
print(direction_summary)

message("\nEligibility summary:")
print(
  paired_station_trends[
    ,
    .(
      station_element_combinations = as.integer(.N),
      eligible_combinations = as.integer(sum(eligible)),
      ineligible_combinations = as.integer(sum(!eligible))
    ),
    by = .(
      period_id,
      element
    )
  ][
    order(period_id, element)
  ]
)

message("\nMaximum paired OLS identity error:")
print(maximum_identity_error)
