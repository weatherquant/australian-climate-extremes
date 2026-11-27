# ============================================================
# 03_build_paired_network_series.R
#
# Build paired raw and homogenised ACORN-SAT network temperature
# anomaly series and estimate the effect of homogenisation on
# Australian long-term temperature trends.
#
# Design principles:
#   1. Raw and homogenised values always use identical station-years.
#   2. Station anomalies use a common 1961-1990 climatology computed
#      separately for raw and homogenised values from identical years.
#   3. The primary network is frozen separately for each requested
#      analysis period using explicit endpoint and coverage rules.
#   4. A latitude-band area-weighted national mean is primary.
#   5. An equal-station-weight national mean is retained as sensitivity.
#   6. Paired station bootstrap intervals quantify uncertainty in the
#      homogenisation-induced difference in national trends.
#
# Important scope:
#   This script estimates the effect of homogenisation within the
#   matched ACORN-SAT station network. It does not reconstruct the
#   Bureau's operational gridded national temperature product.
# ============================================================

# ============================================================
# Configuration
# ============================================================

annual_input_file <- file.path(
  "data-derived",
  "temperature_homogenisation",
  "temperature_raw_homogenised_paired_annual_v2.3.parquet"
)

station_metadata_file <- file.path(
  "data-derived",
  "station_metadata",
  "adam_stations.parquet"
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

climatology_start_year <- 1961L
climatology_end_year <- 1990L
minimum_climatology_years <- 24L

minimum_period_coverage_fraction <- 0.80
minimum_period_valid_years <- 20L
endpoint_tolerance_years <- 5L

latitude_band_width_degrees <- 5
minimum_annual_station_fraction <- 0.80

bootstrap_replicates <- 1000L
random_seed <- 20261127L
confidence_level <- 0.95

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
    call. = FALSE
  )
}

library(data.table)
library(arrow)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  annual_input_file,
  station_metadata_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

# ============================================================
# Helpers
# ============================================================

automatic_hac_lag <- function(n) {
  as.integer(
    min(
      floor(4 * (n / 100)^(2 / 9)),
      max(n - 2L, 0L)
    )
  )
}

fit_trend <- function(year, value, confidence_level = 0.95) {
  valid <- is.finite(year) & is.finite(value)
  
  dat <- data.frame(
    year = as.numeric(year[valid]),
    value = as.numeric(value[valid])
  )
  
  n <- nrow(dat)
  
  if (n < 3L || length(unique(dat$year)) < 2L) {
    return(list(
      n_years = as.integer(n),
      trend_c_per_decade = NA_real_,
      hac_se_c_per_decade = NA_real_,
      conf_low_c_per_decade = NA_real_,
      conf_high_c_per_decade = NA_real_,
      p_value = NA_real_,
      hac_lag = NA_integer_,
      model = NULL
    ))
  }
  
  fit <- stats::lm(value ~ year, data = dat)
  lag <- automatic_hac_lag(n)
  
  covariance <- sandwich::NeweyWest(
    fit,
    lag = lag,
    prewhite = FALSE,
    adjust = TRUE
  )
  
  beta <- stats::coef(fit)[["year"]]
  se <- sqrt(covariance["year", "year"])
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  
  list(
    n_years = as.integer(n),
    trend_c_per_decade = as.numeric(beta * 10),
    hac_se_c_per_decade = as.numeric(se * 10),
    conf_low_c_per_decade = as.numeric((beta - z * se) * 10),
    conf_high_c_per_decade = as.numeric((beta + z * se) * 10),
    p_value = as.numeric(
      2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
    ),
    hac_lag = as.integer(lag),
    model = fit
  )
}

weighted_mean_safe <- function(x, w) {
  valid <- is.finite(x) & is.finite(w) & w > 0
  
  if (!any(valid)) {
    return(NA_real_)
  }
  
  as.numeric(stats::weighted.mean(x[valid], w[valid]))
}

normalise_weights <- function(x) {
  if (length(x) == 0L || !all(is.finite(x)) || sum(x) <= 0) {
    return(rep(NA_real_, length(x)))
  }
  
  as.numeric(x / sum(x))
}

resolve_metadata_column <- function(candidates, available, label) {
  match <- candidates[candidates %in% available]
  
  if (length(match) == 0L) {
    stop(
      "Could not identify ", label,
      " in station metadata. Available columns: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  
  match[[1L]]
}

build_network_series <- function(
    station_year_data,
    station_network,
    weighting_method,
    period_start,
    period_end) {
  
  if (!weighting_method %in% c("latitude_band", "equal_station")) {
    stop("Unknown weighting method.", call. = FALSE)
  }
  
  network <- copy(station_network)
  
  if (weighting_method == "equal_station") {
    network[, station_weight := 1 / .N]
  }
  
  if (weighting_method == "latitude_band") {
    band_counts <- network[, .N, by = latitude_band_id]
    
    network[band_counts, on = "latitude_band_id", stations_in_band := i.N]
    
    network[, station_weight :=
              latitude_band_area_weight /
              stations_in_band
    ]
    
    network[, station_weight := normalise_weights(station_weight)]
  }
  
  x <- station_year_data[
    network[, .(
      acorn_station_number,
      station_weight
    )],
    on = "acorn_station_number",
    nomatch = 0L
  ]
  
  expected_stations <- nrow(network)
  expected_weight <- sum(network$station_weight)
  
  result <- x[
    year >= period_start & year <= period_end,
    .(
      expected_stations = as.integer(expected_stations),
      available_stations = as.integer(sum(
        is.finite(raw_anomaly) & is.finite(hom_anomaly)
      )),
      station_fraction = as.numeric(
        sum(is.finite(raw_anomaly) & is.finite(hom_anomaly)) /
          expected_stations
      ),
      available_weight = as.numeric(sum(
        station_weight[
          is.finite(raw_anomaly) & is.finite(hom_anomaly)
        ],
        na.rm = TRUE
      )),
      expected_weight = as.numeric(expected_weight),
      raw_anomaly = weighted_mean_safe(raw_anomaly, station_weight),
      hom_anomaly = weighted_mean_safe(hom_anomaly, station_weight),
      adjustment_anomaly = weighted_mean_safe(
        hom_anomaly - raw_anomaly,
        station_weight
      )
    ),
    by = year
  ]
  
  result[, `:=`(
    weighting_method = weighting_method,
    annual_coverage_adequate =
      station_fraction >= minimum_annual_station_fraction &
      available_weight >= minimum_annual_station_fraction
  )]
  
  setorder(result, year)
  result
}

build_replicate_series <- function(
    station_year_data,
    station_network,
    sampled_station_numbers,
    weighting_method,
    period_start,
    period_end) {
  
  multiplicities <- data.table(
    acorn_station_number = sampled_station_numbers
  )[, .N, by = acorn_station_number]
  
  network <- station_network[
    multiplicities,
    on = "acorn_station_number",
    nomatch = 0L
  ]
  
  if (weighting_method == "equal_station") {
    network[, replicate_weight := N]
  }
  
  if (weighting_method == "latitude_band") {
    # Preserve the fixed latitude-band area contribution, while
    # resampling stations within the full paired network.
    band_multiplicity <- network[
      , sum(N),
      by = latitude_band_id
    ]
    
    network[
      band_multiplicity,
      on = "latitude_band_id",
      band_sample_size := i.V1
    ]
    
    network[, replicate_weight :=
              latitude_band_area_weight *
              N /
              band_sample_size
    ]
  }
  
  network[, replicate_weight := normalise_weights(replicate_weight)]
  
  x <- station_year_data[
    network[, .(
      acorn_station_number,
      replicate_weight
    )],
    on = "acorn_station_number",
    nomatch = 0L
  ]
  
  replicate_series <- x[
    year >= period_start & year <= period_end,
    .(
      raw_anomaly = weighted_mean_safe(raw_anomaly, replicate_weight),
      hom_anomaly = weighted_mean_safe(hom_anomaly, replicate_weight),
      adjustment_anomaly = weighted_mean_safe(
        hom_anomaly - raw_anomaly,
        replicate_weight
      )
    ),
    by = year
  ]
  
  setorder(replicate_series, year)
  replicate_series
}

# ============================================================
# Read paired annual temperatures and station metadata
# ============================================================

annual <- as.data.table(read_parquet(annual_input_file))
metadata <- as.data.table(read_parquet(station_metadata_file))

required_annual_columns <- c(
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

missing_annual_columns <- setdiff(required_annual_columns, names(annual))

if (length(missing_annual_columns) > 0L) {
  stop(
    "Annual paired file is missing required columns: ",
    paste(missing_annual_columns, collapse = ", "),
    call. = FALSE
  )
}

station_id_column <- resolve_metadata_column(
  c("station_number", "station_id", "site_number"),
  names(metadata),
  "station identifier"
)

latitude_column <- resolve_metadata_column(
  c("latitude", "lat", "station_latitude"),
  names(metadata),
  "latitude"
)

longitude_column <- resolve_metadata_column(
  c("longitude", "lon", "long", "station_longitude"),
  names(metadata),
  "longitude"
)

state_column <- resolve_metadata_column(
  c("state", "state_code", "state_abbreviation"),
  names(metadata),
  "state"
)

metadata_core <- metadata[
  , .(
    acorn_station_number = as.character(get(station_id_column)),
    latitude = as.numeric(get(latitude_column)),
    longitude = as.numeric(get(longitude_column)),
    state = as.character(get(state_column))
  )
]

metadata_core <- metadata_core[
  is.finite(latitude) & is.finite(longitude)
]

metadata_duplicates <- metadata_core[
  , .N,
  by = acorn_station_number
][N > 1L]

if (nrow(metadata_duplicates) > 0L) {
  stop(
    "Station metadata contains duplicated station identifiers for ACORN locations.",
    call. = FALSE
  )
}

acorn_locations <- unique(
  annual[, .(
    acorn_station_number,
    station_name
  )]
)

station_metadata <- metadata_core[
  acorn_locations,
  on = "acorn_station_number"
]

missing_coordinates <- station_metadata[
  !is.finite(latitude) | !is.finite(longitude)
]

if (nrow(missing_coordinates) > 0L) {
  fwrite(
    missing_coordinates,
    file.path(
      table_output_dir,
      "temperature_network_missing_coordinates.csv"
    )
  )
  
  stop(
    nrow(missing_coordinates),
    " ACORN locations lack coordinates in station metadata.",
    call. = FALSE
  )
}

# Five-degree latitude bands spanning all ACORN locations.
minimum_latitude_edge <-
  floor(min(station_metadata$latitude) / latitude_band_width_degrees) *
  latitude_band_width_degrees

maximum_latitude_edge <-
  ceiling(max(station_metadata$latitude) / latitude_band_width_degrees) *
  latitude_band_width_degrees

latitude_breaks <- seq(
  minimum_latitude_edge,
  maximum_latitude_edge,
  by = latitude_band_width_degrees
)

station_metadata[, latitude_band_id :=
                   cut(
                     latitude,
                     breaks = latitude_breaks,
                     include.lowest = TRUE,
                     right = FALSE,
                     labels = FALSE
                   )
]

latitude_band_table <- data.table(
  latitude_band_id = seq_len(length(latitude_breaks) - 1L),
  latitude_lower = latitude_breaks[-length(latitude_breaks)],
  latitude_upper = latitude_breaks[-1L]
)

latitude_band_table[, latitude_band_area_weight :=
                      abs(
                        sin(latitude_upper * pi / 180) -
                          sin(latitude_lower * pi / 180)
                      )
]

# Only retain bands represented by ACORN locations, then normalise.
represented_bands <- unique(station_metadata$latitude_band_id)
latitude_band_table <- latitude_band_table[
  latitude_band_id %in% represented_bands
]

latitude_band_table[, latitude_band_area_weight :=
                      normalise_weights(latitude_band_area_weight)
]

station_metadata <- latitude_band_table[
  station_metadata,
  on = "latitude_band_id"
]

# ============================================================
# Build annual station anomalies
# ============================================================

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

climatology <- annual_long[
  paired_value_available == TRUE &
    year >= climatology_start_year &
    year <= climatology_end_year,
  .(
    climatology_years = as.integer(.N),
    raw_climatology = as.numeric(mean(raw_value)),
    hom_climatology = as.numeric(mean(hom_value))
  ),
  by = .(
    acorn_station_number,
    station_name,
    element
  )
]

climatology[, climatology_valid :=
              climatology_years >= minimum_climatology_years
]

station_year_anomalies <- climatology[
  annual_long,
  on = .(
    acorn_station_number,
    station_name,
    element
  )
]

station_year_anomalies[, `:=`(
  raw_anomaly = data.table::fifelse(
    paired_value_available & climatology_valid,
    raw_value - raw_climatology,
    NA_real_
  ),
  hom_anomaly = data.table::fifelse(
    paired_value_available & climatology_valid,
    hom_value - hom_climatology,
    NA_real_
  )
)]

station_year_anomalies[, adjustment_anomaly :=
                         hom_anomaly - raw_anomaly
]

station_year_anomalies <- station_metadata[
  station_year_anomalies,
  on = "acorn_station_number"
]

setorder(
  station_year_anomalies,
  element,
  acorn_station_number,
  year
)

# ============================================================
# Freeze eligible station networks by period and element
# ============================================================

network_list <- list()
network_index <- 0L

for (period_row in seq_len(nrow(period_definitions))) {
  period_id <- period_definitions$period_id[period_row]
  period_start <- period_definitions$period_start[period_row]
  period_end <- period_definitions$period_end[period_row]
  requested_years <- period_end - period_start + 1L
  
  for (selected_element in c("tmax", "tmin", "tmean")) {
    eligibility <- station_year_anomalies[
      element == selected_element &
        year >= period_start &
        year <= period_end,
      {
        paired_years <- year[
          is.finite(raw_anomaly) & is.finite(hom_anomaly)
        ]
        
        valid_years <- length(paired_years)
        
        list(
          valid_years = as.integer(valid_years),
          coverage_fraction = as.numeric(valid_years / requested_years),
          first_valid_year = if (valid_years > 0L) {
            min(paired_years)
          } else {
            NA_integer_
          },
          last_valid_year = if (valid_years > 0L) {
            max(paired_years)
          } else {
            NA_integer_
          }
        )
      },
      by = .(
        acorn_station_number,
        station_name,
        latitude,
        longitude,
        state,
        latitude_band_id,
        latitude_lower,
        latitude_upper,
        latitude_band_area_weight
      )
    ]
    
    eligibility[, `:=`(
      period_id = period_id,
      period_start = period_start,
      period_end = period_end,
      element = selected_element,
      requested_years = requested_years,
      start_covered =
        is.finite(first_valid_year) &
        first_valid_year <= period_start + endpoint_tolerance_years,
      end_covered =
        is.finite(last_valid_year) &
        last_valid_year >= period_end - endpoint_tolerance_years
    )]
    
    eligibility[, eligible_network :=
                  valid_years >= minimum_period_valid_years &
                  coverage_fraction >= minimum_period_coverage_fraction &
                  start_covered &
                  end_covered
    ]
    
    network_index <- network_index + 1L
    network_list[[network_index]] <- eligibility
  }
}

network_audit <- rbindlist(network_list, use.names = TRUE, fill = TRUE)
setorder(network_audit, period_start, element, acorn_station_number)

frozen_network <- network_audit[eligible_network == TRUE]

# ============================================================
# Build paired national network series
# ============================================================

series_list <- list()
series_index <- 0L

for (period_row in seq_len(nrow(period_definitions))) {
  period_id <- period_definitions$period_id[period_row]
  period_start <- period_definitions$period_start[period_row]
  period_end <- period_definitions$period_end[period_row]
  
  for (selected_element in c("tmax", "tmin", "tmean")) {
    network <- frozen_network[
      period_id == period_definitions$period_id[period_row] &
        element == selected_element
    ]
    
    if (nrow(network) == 0L) {
      next
    }
    
    station_data <- station_year_anomalies[
      element == selected_element
    ]
    
    for (weighting_method in c("latitude_band", "equal_station")) {
      series_index <- series_index + 1L
      
      series <- build_network_series(
        station_year_data = station_data,
        station_network = network,
        weighting_method = weighting_method,
        period_start = period_start,
        period_end = period_end
      )
      
      series[, `:=`(
        period_id = period_id,
        period_start = period_start,
        period_end = period_end,
        element = selected_element
      )]
      
      series_list[[series_index]] <- series
    }
  }
}

network_series <- rbindlist(series_list, use.names = TRUE, fill = TRUE)
setorder(network_series, period_start, element, weighting_method, year)

# ============================================================
# Fit raw, homogenised and paired-difference trends
# ============================================================

trend_results <- list()
trend_models <- list()
trend_index <- 0L

trend_groups <- unique(
  network_series[, .(
    period_id,
    period_start,
    period_end,
    element,
    weighting_method
  )]
)

for (group_row in seq_len(nrow(trend_groups))) {
  group <- trend_groups[group_row]
  
  dat <- network_series[
    period_id == group$period_id &
      element == group$element &
      weighting_method == group$weighting_method &
      annual_coverage_adequate == TRUE
  ]
  
  raw_fit <- fit_trend(dat$year, dat$raw_anomaly, confidence_level)
  hom_fit <- fit_trend(dat$year, dat$hom_anomaly, confidence_level)
  difference_fit <- fit_trend(
    dat$year,
    dat$adjustment_anomaly,
    confidence_level
  )
  
  trend_difference_identity_error <-
    hom_fit$trend_c_per_decade -
    raw_fit$trend_c_per_decade -
    difference_fit$trend_c_per_decade
  
  trend_index <- trend_index + 1L
  
  trend_results[[trend_index]] <- data.table(
    period_id = group$period_id,
    period_start = group$period_start,
    period_end = group$period_end,
    element = group$element,
    weighting_method = group$weighting_method,
    expected_stations = unique(dat$expected_stations),
    valid_network_years = raw_fit$n_years,
    minimum_annual_station_fraction = min(dat$station_fraction),
    median_annual_station_fraction = median(dat$station_fraction),
    minimum_annual_available_weight = min(dat$available_weight),
    median_annual_available_weight = median(dat$available_weight),
    
    raw_trend_c_per_decade = raw_fit$trend_c_per_decade,
    raw_hac_se_c_per_decade = raw_fit$hac_se_c_per_decade,
    raw_conf_low_c_per_decade = raw_fit$conf_low_c_per_decade,
    raw_conf_high_c_per_decade = raw_fit$conf_high_c_per_decade,
    raw_p_value = raw_fit$p_value,
    
    hom_trend_c_per_decade = hom_fit$trend_c_per_decade,
    hom_hac_se_c_per_decade = hom_fit$hac_se_c_per_decade,
    hom_conf_low_c_per_decade = hom_fit$conf_low_c_per_decade,
    hom_conf_high_c_per_decade = hom_fit$conf_high_c_per_decade,
    hom_p_value = hom_fit$p_value,
    
    trend_difference_c_per_decade =
      difference_fit$trend_c_per_decade,
    difference_hac_se_c_per_decade =
      difference_fit$hac_se_c_per_decade,
    difference_conf_low_c_per_decade =
      difference_fit$conf_low_c_per_decade,
    difference_conf_high_c_per_decade =
      difference_fit$conf_high_c_per_decade,
    difference_p_value = difference_fit$p_value,
    difference_hac_lag = difference_fit$hac_lag,
    trend_difference_identity_error =
      trend_difference_identity_error
  )
  
  model_key <- paste(
    group$period_id,
    group$element,
    group$weighting_method,
    sep = "__"
  )
  
  trend_models[[model_key]] <- list(
    raw = raw_fit$model,
    homogenised = hom_fit$model,
    difference = difference_fit$model
  )
}

network_trends <- rbindlist(trend_results, use.names = TRUE, fill = TRUE)
setorder(network_trends, period_start, element, weighting_method)

maximum_identity_error <- network_trends[
  , max(abs(trend_difference_identity_error), na.rm = TRUE)
]

if (!is.finite(maximum_identity_error) || maximum_identity_error > 1e-10) {
  stop(
    "Network paired trend identity failed. Maximum error = ",
    maximum_identity_error,
    call. = FALSE
  )
}

# ============================================================
# Paired station bootstrap of national trend differences
# ============================================================

set.seed(random_seed)

bootstrap_results <- list()
bootstrap_index <- 0L

for (group_row in seq_len(nrow(trend_groups))) {
  group <- trend_groups[group_row]
  
  network <- frozen_network[
    period_id == group$period_id &
      element == group$element
  ]
  
  station_data <- station_year_anomalies[
    element == group$element
  ]
  
  station_numbers <- network$acorn_station_number
  station_count <- length(station_numbers)
  
  replicate_effects <- rep(NA_real_, bootstrap_replicates)
  
  for (replicate in seq_len(bootstrap_replicates)) {
    if (
      group$weighting_method ==
      "latitude_band"
    ) {
      
      sampled_station_numbers <- network[
        ,
        sample(
          acorn_station_number,
          size = .N,
          replace = TRUE
        ),
        by = latitude_band_id
      ]$V1
      
    } else {
      
      sampled_station_numbers <- sample(
        station_numbers,
        size = station_count,
        replace = TRUE
      )
    }
    
    replicate_series <- build_replicate_series(
      station_year_data = station_data,
      station_network = network,
      sampled_station_numbers = sampled_station_numbers,
      weighting_method = group$weighting_method,
      period_start = group$period_start,
      period_end = group$period_end
    )
    
    replicate_fit <- fit_trend(
      replicate_series$year,
      replicate_series$adjustment_anomaly,
      confidence_level
    )
    
    replicate_effects[replicate] <-
      replicate_fit$trend_c_per_decade
  }
  
  point_estimate <- network_trends[
    period_id == group$period_id &
      element == group$element &
      weighting_method == group$weighting_method,
    trend_difference_c_per_decade
  ]
  
  successful <- replicate_effects[is.finite(replicate_effects)]
  
  bootstrap_index <- bootstrap_index + 1L
  
  bootstrap_results[[bootstrap_index]] <- data.table(
    period_id = group$period_id,
    period_start = group$period_start,
    period_end = group$period_end,
    element = group$element,
    weighting_method = group$weighting_method,
    stations = as.integer(station_count),
    requested_replicates = bootstrap_replicates,
    successful_replicates = as.integer(length(successful)),
    point_estimate_c_per_decade = as.numeric(point_estimate),
    bootstrap_mean_c_per_decade = as.numeric(mean(successful)),
    bootstrap_median_c_per_decade = as.numeric(median(successful)),
    bootstrap_sd_c_per_decade = as.numeric(stats::sd(successful)),
    bootstrap_conf_low_c_per_decade = as.numeric(
      stats::quantile(successful, 0.025, names = FALSE)
    ),
    bootstrap_conf_high_c_per_decade = as.numeric(
      stats::quantile(successful, 0.975, names = FALSE)
    ),
    bootstrap_bias_c_per_decade =
      as.numeric(
        mean(successful) -
          point_estimate
      ),
    bootstrap_median_bias_c_per_decade =
      as.numeric(
        median(successful) -
          point_estimate
      ),
    proportion_positive = as.numeric(mean(successful > 0)),
    proportion_negative = as.numeric(mean(successful < 0))
  )
}

network_bootstrap <- rbindlist(
  bootstrap_results,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  network_bootstrap,
  period_start,
  element,
  weighting_method
)

# ============================================================
# Coverage and network audit summaries
# ============================================================

network_summary <- network_audit[
  , .(
    stations_considered = as.integer(.N),
    eligible_stations = as.integer(sum(eligible_network)),
    median_valid_years = as.numeric(median(valid_years)),
    minimum_valid_years = as.integer(min(valid_years)),
    median_coverage_fraction = as.numeric(median(coverage_fraction)),
    latitude_bands = as.integer(uniqueN(latitude_band_id))
  ),
  by = .(
    period_id,
    period_start,
    period_end,
    element
  )
]

setorder(network_summary, period_start, element)

climatology_audit <- climatology[
  , .(
    stations = as.integer(.N),
    valid_climatologies = as.integer(sum(climatology_valid)),
    minimum_climatology_years = as.integer(min(climatology_years)),
    median_climatology_years = as.numeric(median(climatology_years))
  ),
  by = element
]

# ============================================================
# Save outputs
# ============================================================

arrow::write_parquet(
  station_year_anomalies,
  file.path(
    output_dir,
    "temperature_paired_station_year_anomalies_v2.3.parquet"
  ),
  compression = "zstd"
)

arrow::write_parquet(
  network_series,
  file.path(
    output_dir,
    "temperature_paired_network_series_v2.3.parquet"
  ),
  compression = "zstd"
)

arrow::write_parquet(
  network_trends,
  file.path(
    output_dir,
    "temperature_paired_network_trends_v2.3.parquet"
  ),
  compression = "zstd"
)

fwrite(
  network_series,
  file.path(
    table_output_dir,
    "temperature_paired_network_series_v2.3.csv"
  )
)

fwrite(
  network_trends,
  file.path(
    table_output_dir,
    "temperature_paired_network_trends_v2.3.csv"
  )
)

fwrite(
  network_bootstrap,
  file.path(
    table_output_dir,
    "temperature_paired_network_trend_bootstrap_v2.3.csv"
  )
)

fwrite(
  network_summary,
  file.path(
    table_output_dir,
    "temperature_paired_network_definition_v2.3.csv"
  )
)

fwrite(
  network_audit,
  file.path(
    table_output_dir,
    "temperature_paired_network_station_audit_v2.3.csv"
  )
)

fwrite(
  climatology_audit,
  file.path(
    table_output_dir,
    "temperature_paired_climatology_audit_v2.3.csv"
  )
)

fwrite(
  latitude_band_table,
  file.path(
    table_output_dir,
    "temperature_latitude_band_weights_v2.3.csv"
  )
)

saveRDS(
  list(
    models = trend_models,
    configuration = list(
      periods = period_definitions,
      climatology_start_year = climatology_start_year,
      climatology_end_year = climatology_end_year,
      minimum_climatology_years = minimum_climatology_years,
      minimum_period_coverage_fraction =
        minimum_period_coverage_fraction,
      minimum_period_valid_years = minimum_period_valid_years,
      endpoint_tolerance_years = endpoint_tolerance_years,
      latitude_band_width_degrees = latitude_band_width_degrees,
      minimum_annual_station_fraction =
        minimum_annual_station_fraction,
      bootstrap_replicates = bootstrap_replicates,
      random_seed = random_seed,
      confidence_level = confidence_level
    ),
    maximum_paired_trend_identity_error = maximum_identity_error
  ),
  file.path(
    output_dir,
    "temperature_paired_network_models_v2.3.rds"
  ),
  compress = "xz"
)

# ============================================================
# Console summary
# ============================================================

message("\n========================================")
message("PAIRED TEMPERATURE NETWORK SERIES COMPLETE")
message("========================================\n")

message("Climatology audit:")
print(climatology_audit)

message("\nFrozen network summary:")
print(network_summary)

message("\nNational paired trends:")
print(network_trends)

message("\nPaired station-bootstrap trend differences:")
print(network_bootstrap)

message("\nMaximum paired trend identity error:")
print(maximum_identity_error)
