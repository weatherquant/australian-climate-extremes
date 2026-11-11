# ============================================================
# build_extreme_rainfall_metrics.R
#
# Build annual extreme-rainfall metrics for the frozen
# continental Australian rainfall network.
#
# Default analysis:
#   Period:                  1950-2025
#   Reference period:        1961-1990
#   Wet-day definition:      rainfall >= 1 mm
#   Annual completeness:     reported, not filtered here
#
# One-day extreme eligibility:
#
#   rainfall is observed
#   AND rainfall_period is missing or equals 1
#   AND rainfall_day_count is missing or <= 1
#
# Multi-day and non-standard-period records are retained in
# audit counts but excluded from Rx1day, percentile thresholds
# and fixed daily-threshold exceedances.
#
# Important:
#   This uses fixed station-specific reference thresholds.
#   It is not yet the bootstrap-adjusted ETCCDI implementation.
# ============================================================


# ============================================================
# Package check
# ============================================================

check_extreme_rainfall_packages <- function() {
  
  required_packages <- c(
    "data.table",
    "arrow"
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
      paste(
        missing_packages,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Calendar days in year
# ============================================================

rainfall_days_in_year <- function(year) {
  
  year <- as.integer(year)
  
  as.integer(
    as.Date(
      paste0(year + 1L, "-01-01")
    ) -
      as.Date(
        paste0(year, "-01-01")
      )
  )
}


# ============================================================
# Validate and prepare core network
# ============================================================

prepare_extreme_rainfall_network <- function(
    core_network) {
  
  required <- c(
    "station_number",
    "station_name",
    "latitude",
    "longitude",
    "elevation_m",
    "state",
    "rainfall_district",
    "selection_start_year",
    "selection_end_year"
  )
  
  missing_columns <- setdiff(
    required,
    names(core_network)
  )
  
  if (length(missing_columns) > 0L) {
    
    stop(
      "Core rainfall network is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  network <- data.table::copy(
    core_network
  )
  
  if (
    "continental_core" %in% names(network)
  ) {
    
    network <- network[
      continental_core == TRUE
    ]
  }
  
  duplicate_stations <- network[
    ,
    .N,
    by = station_number
  ][
    N > 1L
  ]
  
  if (nrow(duplicate_stations) > 0L) {
    
    stop(
      "Core network contains duplicated station numbers.",
      call. = FALSE
    )
  }
  
  if (
    any(
      is.na(network$station_number)
    )
  ) {
    
    stop(
      "Core network contains missing station numbers.",
      call. = FALSE
    )
  }
  
  data.table::setorder(
    network,
    station_number
  )
  
  network
}


# ============================================================
# Daily one-day eligibility
# ============================================================

add_rainfall_eligibility_flags <- function(
    rainfall,
    wet_day_threshold_mm = 1) {
  
  required <- c(
    "station_number",
    "date",
    "rainfall_mm",
    "observation_quality_flag",
    "rainfall_period",
    "rainfall_day_count",
    "rainfall_quality",
    "aws_flag"
  )
  
  missing_columns <- setdiff(
    required,
    names(rainfall)
  )
  
  if (length(missing_columns) > 0L) {
    
    stop(
      "Daily rainfall data are missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  x <- data.table::copy(
    rainfall
  )
  
  x[
    ,
    rainfall_observed :=
      !is.na(rainfall_mm)
  ]
  
  x[
    ,
    multiday_record :=
      rainfall_observed &
      !is.na(rainfall_day_count) &
      rainfall_day_count > 1
  ]
  
  x[
    ,
    nonstandard_period :=
      rainfall_observed &
      !is.na(rainfall_period) &
      rainfall_period != 1
  ]
  
  x[
    ,
    eligible_1day :=
      rainfall_observed &
      (
        is.na(rainfall_period) |
          rainfall_period == 1
      ) &
      (
        is.na(rainfall_day_count) |
          rainfall_day_count <= 1
      )
  ]
  
  x[
    ,
    valid_wet_day :=
      eligible_1day &
      rainfall_mm >= wet_day_threshold_mm
  ]
  
  x[
    ,
    quality_flag_nonzero :=
      rainfall_observed &
      (
        (
          !is.na(observation_quality_flag) &
            observation_quality_flag != 0
        ) |
          (
            !is.na(rainfall_quality) &
              rainfall_quality != 0
          )
      )
  ]
  
  x
}


# ============================================================
# Discover Parquet files overlapping a year interval
# ============================================================

rainfall_files_for_period <- function(
    parquet_files,
    start_year,
    end_year) {
  
  if (length(parquet_files) == 0L) {
    
    stop(
      "No rainfall Parquet files supplied.",
      call. = FALSE
    )
  }
  
  filenames <- basename(
    parquet_files
  )
  
  match <- regexec(
    "^adam_prcp_days_([0-9]{4})_([0-9]{4})[.]parquet$",
    filenames
  )
  
  parts <- regmatches(
    filenames,
    match
  )
  
  valid <- lengths(parts) == 3L
  
  file_periods <- data.table::data.table(
    path = parquet_files[valid],
    filename = filenames[valid],
    file_start_year = as.integer(
      vapply(
        parts[valid],
        `[[`,
        character(1),
        2L
      )
    ),
    file_end_year = as.integer(
      vapply(
        parts[valid],
        `[[`,
        character(1),
        3L
      )
    )
  )
  
  selected <- file_periods[
    file_end_year >= start_year &
      file_start_year <= end_year
  ]
  
  data.table::setorder(
    selected,
    file_start_year
  )
  
  selected$path
}


# ============================================================
# Read one rainfall Parquet chunk for selected stations
# ============================================================

read_extreme_rainfall_chunk <- function(
    path,
    station_numbers,
    start_year,
    end_year) {
  
  if (!file.exists(path)) {
    
    stop(
      "Rainfall Parquet file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  rainfall <- data.table::as.data.table(
    arrow::read_parquet(
      path,
      col_select = c(
        "station_number",
        "date",
        "rainfall_mm",
        "observation_quality_flag",
        "rainfall_period",
        "rainfall_type",
        "rainfall_quality",
        "rainfall_day_count",
        "aws_flag",
        "source_file"
      )
    )
  )
  
  rainfall <- rainfall[
    station_number %in% station_numbers
  ]
  
  rainfall[
    ,
    year :=
      as.integer(
        format(
          date,
          "%Y"
        )
      )
  ]
  
  rainfall <- rainfall[
    year >= start_year &
      year <= end_year
  ]
  
  rainfall
}


# ============================================================
# Build station-specific reference thresholds
# ============================================================

build_rainfall_reference_thresholds <- function(
    parquet_files,
    core_network,
    reference_start_year = 1961L,
    reference_end_year = 1990L,
    wet_day_threshold_mm = 1,
    minimum_reference_wet_days = 300L,
    quantile_type = 8L) {
  
  check_extreme_rainfall_packages()
  
  network <- prepare_extreme_rainfall_network(
    core_network
  )
  
  station_numbers <- network$station_number
  
  reference_files <- rainfall_files_for_period(
    parquet_files = parquet_files,
    start_year = reference_start_year,
    end_year = reference_end_year
  )
  
  if (length(reference_files) == 0L) {
    
    stop(
      "No rainfall Parquet files overlap the reference period.",
      call. = FALSE
    )
  }
  
  wet_day_tables <- vector(
    mode = "list",
    length = length(reference_files)
  )
  
  for (i in seq_along(reference_files)) {
    
    message(
      "\nReference file ",
      i,
      " of ",
      length(reference_files),
      ": ",
      basename(reference_files[[i]])
    )
    
    rainfall <- read_extreme_rainfall_chunk(
      path = reference_files[[i]],
      station_numbers = station_numbers,
      start_year = reference_start_year,
      end_year = reference_end_year
    )
    
    rainfall <- add_rainfall_eligibility_flags(
      rainfall = rainfall,
      wet_day_threshold_mm = wet_day_threshold_mm
    )
    
    wet_day_tables[[i]] <- rainfall[
      valid_wet_day == TRUE,
      .(
        station_number,
        rainfall_mm
      )
    ]
    
    rm(rainfall)
    
    invisible(
      gc()
    )
  }
  
  reference_wet_days <- data.table::rbindlist(
    wet_day_tables,
    use.names = TRUE,
    fill = TRUE
  )
  
  thresholds <- reference_wet_days[
    ,
    .(
      n_reference_wet_days =
        .N,
      
      p90_reference_mm =
        stats::quantile(
          rainfall_mm,
          probs = 0.90,
          na.rm = TRUE,
          type = quantile_type,
          names = FALSE
        ),
      
      p95_reference_mm =
        stats::quantile(
          rainfall_mm,
          probs = 0.95,
          na.rm = TRUE,
          type = quantile_type,
          names = FALSE
        ),
      
      p99_reference_mm =
        stats::quantile(
          rainfall_mm,
          probs = 0.99,
          na.rm = TRUE,
          type = quantile_type,
          names = FALSE
        )
    ),
    by = station_number
  ]
  
  thresholds[
    ,
    `:=`(
      reference_start_year =
        reference_start_year,
      
      reference_end_year =
        reference_end_year,
      
      wet_day_threshold_mm =
        wet_day_threshold_mm,
      
      quantile_type =
        quantile_type,
      
      minimum_reference_wet_days =
        minimum_reference_wet_days,
      
      reference_threshold_valid =
        n_reference_wet_days >=
        minimum_reference_wet_days
    )
  ]
  
  thresholds[
    reference_threshold_valid == FALSE,
    `:=`(
      p90_reference_mm = NA_real_,
      p95_reference_mm = NA_real_,
      p99_reference_mm = NA_real_
    )
  ]
  
  thresholds <- merge(
    network[
      ,
      .(
        station_number,
        station_name,
        latitude,
        longitude,
        elevation_m,
        state,
        rainfall_district
      )
    ],
    thresholds,
    by = "station_number",
    all.x = TRUE,
    sort = FALSE
  )
  
  thresholds[
    is.na(reference_threshold_valid),
    reference_threshold_valid := FALSE
  ]
  
  data.table::setorder(
    thresholds,
    station_number
  )
  
  rm(reference_wet_days)
  rm(wet_day_tables)
  
  invisible(
    gc()
  )
  
  thresholds
}


# ============================================================
# Summarise one daily chunk to station-year metrics
# ============================================================

summarise_extreme_rainfall_chunk <- function(
    rainfall,
    thresholds,
    wet_day_threshold_mm = 1) {
  
  rainfall <- add_rainfall_eligibility_flags(
    rainfall = rainfall,
    wet_day_threshold_mm = wet_day_threshold_mm
  )
  
  threshold_keep <- thresholds[
    ,
    .(
      station_number,
      p90_reference_mm,
      p95_reference_mm,
      p99_reference_mm,
      reference_threshold_valid
    )
  ]
  
  rainfall <- merge(
    rainfall,
    threshold_keep,
    by = "station_number",
    all.x = TRUE,
    sort = FALSE
  )
  
  rainfall[
    ,
    above_p90 :=
      valid_wet_day &
      reference_threshold_valid &
      !is.na(p90_reference_mm) &
      rainfall_mm > p90_reference_mm
  ]
  
  rainfall[
    ,
    above_p95 :=
      valid_wet_day &
      reference_threshold_valid &
      !is.na(p95_reference_mm) &
      rainfall_mm > p95_reference_mm
  ]
  
  rainfall[
    ,
    above_p99 :=
      valid_wet_day &
      reference_threshold_valid &
      !is.na(p99_reference_mm) &
      rainfall_mm > p99_reference_mm
  ]
  
  annual <- rainfall[
    ,
    .(
      first_record_date =
        min(date),
      
      last_record_date =
        max(date),
      
      n_rows =
        .N,
      
      expected_days =
        rainfall_days_in_year(
          data.table::first(year)
        ),
      
      n_observed =
        sum(rainfall_observed),
      
      n_missing =
        sum(!rainfall_observed),
      
      completeness =
        sum(rainfall_observed) /
        rainfall_days_in_year(
          data.table::first(year)
        ),
      
      n_valid_1day =
        sum(eligible_1day),
      
      n_invalid_1day =
        sum(
          rainfall_observed &
            !eligible_1day
        ),
      
      n_multiday_excluded =
        sum(multiday_record),
      
      n_nonstandard_period_excluded =
        sum(nonstandard_period),
      
      n_nonzero_quality_flag =
        sum(quality_flag_nonzero),
      
      n_aws_records =
        sum(
          aws_flag == 1,
          na.rm = TRUE
        ),
      
      # Complete recorded annual total. This retains
      # accumulations because they contribute to annual rainfall.
      annual_total_recorded_mm =
        if (
          all(is.na(rainfall_mm))
        ) {
          NA_real_
        } else {
          sum(
            rainfall_mm,
            na.rm = TRUE
          )
        },
      
      # Total from records eligible as one-day observations.
      annual_total_valid_1day_mm =
        if (
          !any(eligible_1day)
        ) {
          NA_real_
        } else {
          sum(
            rainfall_mm[
              eligible_1day
            ],
            na.rm = TRUE
          )
        },
      
      dry_days_valid =
        sum(
          eligible_1day &
            rainfall_mm <
            wet_day_threshold_mm,
          na.rm = TRUE
        ),
      
      wet_days_valid =
        sum(valid_wet_day),
      
      mean_wet_day_mm =
        if (
          !any(valid_wet_day)
        ) {
          NA_real_
        } else {
          mean(
            rainfall_mm[
              valid_wet_day
            ],
            na.rm = TRUE
          )
        },
      
      rx1day_mm =
        if (
          !any(eligible_1day)
        ) {
          NA_real_
        } else {
          max(
            rainfall_mm[
              eligible_1day
            ],
            na.rm = TRUE
          )
        },
      
      rx1day_date =
        if (
          !any(eligible_1day)
        ) {
          as.Date(NA)
        } else {
          date[
            eligible_1day
          ][
            which.max(
              rainfall_mm[
                eligible_1day
              ]
            )
          ]
        },
      
      days_ge_10mm =
        sum(
          eligible_1day &
            rainfall_mm >= 10,
          na.rm = TRUE
        ),
      
      days_ge_25mm =
        sum(
          eligible_1day &
            rainfall_mm >= 25,
          na.rm = TRUE
        ),
      
      days_ge_50mm =
        sum(
          eligible_1day &
            rainfall_mm >= 50,
          na.rm = TRUE
        ),
      
      days_ge_100mm =
        sum(
          eligible_1day &
            rainfall_mm >= 100,
          na.rm = TRUE
        ),
      
      annual_wet_day_q90_mm =
        if (
          sum(valid_wet_day) < 10L
        ) {
          NA_real_
        } else {
          stats::quantile(
            rainfall_mm[
              valid_wet_day
            ],
            probs = 0.90,
            na.rm = TRUE,
            type = 8,
            names = FALSE
          )
        },
      
      annual_wet_day_q95_mm =
        if (
          sum(valid_wet_day) < 20L
        ) {
          NA_real_
        } else {
          stats::quantile(
            rainfall_mm[
              valid_wet_day
            ],
            probs = 0.95,
            na.rm = TRUE,
            type = 8,
            names = FALSE
          )
        },
      
      annual_wet_day_q99_mm =
        if (
          sum(valid_wet_day) < 50L
        ) {
          NA_real_
        } else {
          stats::quantile(
            rainfall_mm[
              valid_wet_day
            ],
            probs = 0.99,
            na.rm = TRUE,
            type = 8,
            names = FALSE
          )
        },
      
      days_above_reference_p90 =
        sum(
          above_p90,
          na.rm = TRUE
        ),
      
      days_above_reference_p95 =
        sum(
          above_p95,
          na.rm = TRUE
        ),
      
      days_above_reference_p99 =
        sum(
          above_p99,
          na.rm = TRUE
        ),
      
      rainfall_above_reference_p90_mm =
        if (
          !any(above_p90, na.rm = TRUE)
        ) {
          0
        } else {
          sum(
            rainfall_mm[
              above_p90
            ],
            na.rm = TRUE
          )
        },
      
      rainfall_above_reference_p95_mm =
        if (
          !any(above_p95, na.rm = TRUE)
        ) {
          0
        } else {
          sum(
            rainfall_mm[
              above_p95
            ],
            na.rm = TRUE
          )
        },
      
      rainfall_above_reference_p99_mm =
        if (
          !any(above_p99, na.rm = TRUE)
        ) {
          0
        } else {
          sum(
            rainfall_mm[
              above_p99
            ],
            na.rm = TRUE
          )
        },
      
      p90_reference_mm =
        data.table::first(
          p90_reference_mm
        ),
      
      p95_reference_mm =
        data.table::first(
          p95_reference_mm
        ),
      
      p99_reference_mm =
        data.table::first(
          p99_reference_mm
        ),
      
      reference_threshold_valid =
        data.table::first(
          reference_threshold_valid
        )
    ),
    by = .(
      station_number,
      year
    )
  ]
  
  annual[
    ,
    valid_1day_fraction :=
      data.table::fifelse(
        n_observed > 0,
        n_valid_1day /
          n_observed,
        NA_real_
      )
  ]
  
  annual[
    ,
    wet_day_fraction :=
      data.table::fifelse(
        n_valid_1day > 0,
        wet_days_valid /
          n_valid_1day,
        NA_real_
      )
  ]
  
  annual[
    ,
    annual_complete_95 :=
      completeness >= 0.95
  ]
  
  data.table::setorder(
    annual,
    station_number,
    year
  )
  
  annual
}


# ============================================================
# Build annual extreme-rainfall metrics
# ============================================================

build_extreme_rainfall_metrics <- function(
    parquet_files,
    core_network,
    thresholds,
    analysis_start_year = 1950L,
    analysis_end_year = 2025L,
    wet_day_threshold_mm = 1) {
  
  check_extreme_rainfall_packages()
  
  network <- prepare_extreme_rainfall_network(
    core_network
  )
  
  station_numbers <- network$station_number
  
  analysis_files <- rainfall_files_for_period(
    parquet_files = parquet_files,
    start_year = analysis_start_year,
    end_year = analysis_end_year
  )
  
  if (length(analysis_files) == 0L) {
    
    stop(
      "No rainfall Parquet files overlap the analysis period.",
      call. = FALSE
    )
  }
  
  results <- vector(
    mode = "list",
    length = length(analysis_files)
  )
  
  for (i in seq_along(analysis_files)) {
    
    message(
      "\nAnalysis file ",
      i,
      " of ",
      length(analysis_files),
      ": ",
      basename(analysis_files[[i]])
    )
    
    rainfall <- read_extreme_rainfall_chunk(
      path = analysis_files[[i]],
      station_numbers = station_numbers,
      start_year = analysis_start_year,
      end_year = analysis_end_year
    )
    
    results[[i]] <- summarise_extreme_rainfall_chunk(
      rainfall = rainfall,
      thresholds = thresholds,
      wet_day_threshold_mm = wet_day_threshold_mm
    )
    
    rm(rainfall)
    
    invisible(
      gc()
    )
  }
  
  annual_metrics <- data.table::rbindlist(
    results,
    use.names = TRUE,
    fill = TRUE
  )
  
  duplicate_station_years <- annual_metrics[
    ,
    .N,
    by = .(
      station_number,
      year
    )
  ][
    N > 1L
  ]
  
  if (nrow(duplicate_station_years) > 0L) {
    
    stop(
      "Annual metrics contain duplicated station-year keys.",
      call. = FALSE
    )
  }
  
  network_metadata <- network[
    ,
    .(
      station_number,
      station_name,
      latitude,
      longitude,
      elevation_m,
      state,
      rainfall_district,
      area_name,
      station_open,
      open_date,
      close_date,
      qualifying_years,
      qualifying_fraction,
      selection_start_year,
      selection_end_year,
      selection_rule
    )
  ]
  
  annual_metrics <- merge(
    annual_metrics,
    network_metadata,
    by = "station_number",
    all.x = TRUE,
    sort = FALSE
  )
  
  annual_metrics[
    ,
    `:=`(
      analysis_start_year =
        analysis_start_year,
      
      analysis_end_year =
        analysis_end_year,
      
      wet_day_threshold_mm =
        wet_day_threshold_mm
    )
  ]
  
  data.table::setorder(
    annual_metrics,
    station_number,
    year
  )
  
  annual_metrics
}


# ============================================================
# Validate annual metrics
# ============================================================

validate_extreme_rainfall_metrics <- function(
    annual_metrics,
    expected_stations = NULL,
    start_year = 1950L,
    end_year = 2025L) {
  
  duplicate_station_years <- annual_metrics[
    ,
    .N,
    by = .(
      station_number,
      year
    )
  ][
    N > 1L
  ]
  
  expected_station_years <- if (
    is.null(expected_stations)
  ) {
    NA_integer_
  } else {
    length(
      unique(expected_stations)
    ) *
      (
        end_year -
          start_year +
          1L
      )
  }
  
  data.table::data.table(
    
    rows =
      nrow(annual_metrics),
    
    stations =
      data.table::uniqueN(
        annual_metrics$station_number
      ),
    
    first_year =
      min(
        annual_metrics$year
      ),
    
    last_year =
      max(
        annual_metrics$year
      ),
    
    expected_station_years =
      expected_station_years,
    
    missing_station_years =
      if (
        is.na(expected_station_years)
      ) {
        NA_integer_
      } else {
        expected_station_years -
          nrow(annual_metrics)
      },
    
    duplicate_station_years =
      nrow(
        duplicate_station_years
      ),
    
    complete_95_station_years =
      sum(
        annual_metrics$annual_complete_95,
        na.rm = TRUE
      ),
    
    missing_rx1day =
      sum(
        is.na(
          annual_metrics$rx1day_mm
        )
      ),
    
    invalid_reference_thresholds =
      data.table::uniqueN(
        annual_metrics[
          reference_threshold_valid == FALSE,
          station_number
        ]
      ),
    
    station_years_with_multiday_exclusions =
      sum(
        annual_metrics$
          n_multiday_excluded > 0,
        na.rm = TRUE
      ),
    
    station_years_with_nonstandard_exclusions =
      sum(
        annual_metrics$
          n_nonstandard_period_excluded > 0,
        na.rm = TRUE
      )
  )
}


# ============================================================
# Summarise station-level threshold availability
# ============================================================

summarise_rainfall_thresholds <- function(
    thresholds) {
  
  thresholds[
    ,
    .(
      stations =
        .N,
      
      valid_thresholds =
        sum(
          reference_threshold_valid,
          na.rm = TRUE
        ),
      
      invalid_thresholds =
        sum(
          !reference_threshold_valid |
            is.na(reference_threshold_valid),
          na.rm = TRUE
        ),
      
      minimum_reference_wet_days =
        min(
          n_reference_wet_days,
          na.rm = TRUE
        ),
      
      median_reference_wet_days =
        stats::median(
          n_reference_wet_days,
          na.rm = TRUE
        ),
      
      maximum_reference_wet_days =
        max(
          n_reference_wet_days,
          na.rm = TRUE
        ),
      
      median_p95_mm =
        stats::median(
          p95_reference_mm,
          na.rm = TRUE
        ),
      
      median_p99_mm =
        stats::median(
          p99_reference_mm,
          na.rm = TRUE
        )
    )
  ]
}