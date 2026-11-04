# ============================================================
# summarise_rainfall_availability.R
#
# Build station-year summaries from the national daily
# ADAM PRCP_DAYS Parquet archive.
#
# Designed for large datasets:
#   - reads one Parquet source chunk at a time
#   - summarises immediately to station-year level
#   - releases daily data before reading the next chunk
#
# Daily rainfall interpretation:
#   NA rainfall = observation unavailable
#   0 mm        = valid observed dry day
#   > 0 mm      = observed wet day
#
# No QC filtering is applied at this stage. The purpose is
# to describe the observational archive objectively before
# selecting the research network.
# ============================================================


# ============================================================
# Helper: expected days in a calendar year
# ============================================================

days_in_year <- function(year) {
  
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
# Summarise one Parquet chunk
# ============================================================

summarise_rainfall_parquet_file <- function(path) {
  
  if (!file.exists(path)) {
    
    stop(
      "Rainfall Parquet file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  
  message(
    "Summarising: ",
    basename(path)
  )
  
  
  # ----------------------------------------------------------
  # Read only variables needed for network diagnostics
  # ----------------------------------------------------------
  
  rainfall <- data.table::as.data.table(
    
    arrow::read_parquet(
      
      path,
      
      col_select = c(
        "station_number",
        "date",
        "rainfall_mm",
        "observation_quality_flag",
        "rainfall_period",
        "rainfall_day_count",
        "aws_flag"
      )
    )
  )
  
  
  if (nrow(rainfall) == 0L) {
    
    stop(
      "Rainfall Parquet file contains no rows: ",
      path,
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Calendar year
  # ----------------------------------------------------------
  
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
  
  
  # ==========================================================
  # Station-year summary
  # ==========================================================
  
  station_year <- rainfall[
    ,
    .(
      
      # ------------------------------------------------------
      # Actual rows in database
      # ------------------------------------------------------
      
      first_record_date =
        min(date),
      
      last_record_date =
        max(date),
      
      n_rows =
        .N,
      
      
      # ------------------------------------------------------
      # Rainfall availability
      # ------------------------------------------------------
      
      n_observed =
        sum(
          !is.na(rainfall_mm)
        ),
      
      n_missing =
        sum(
          is.na(rainfall_mm)
        ),
      
      
      # ------------------------------------------------------
      # Dry/wet observations
      # ------------------------------------------------------
      
      n_zero =
        sum(
          rainfall_mm == 0,
          na.rm = TRUE
        ),
      
      n_wet =
        sum(
          rainfall_mm > 0,
          na.rm = TRUE
        ),
      
      
      # ------------------------------------------------------
      # Annual rainfall
      #
      # Retained diagnostically even for incomplete years.
      # Do not interpret as annual total unless completeness
      # criteria have subsequently been met.
      # ------------------------------------------------------
      
      annual_rainfall_mm =
        if (
          all(
            is.na(rainfall_mm)
          )
        ) {
          NA_real_
        } else {
          sum(
            rainfall_mm,
            na.rm = TRUE
          )
        },
      
      
      # ------------------------------------------------------
      # Largest recorded daily value
      #
      # This is not yet our scientifically accepted daily
      # extreme because multi-day accumulation records have
      # not yet been excluded.
      # ------------------------------------------------------
      
      maximum_recorded_rainfall_mm =
        if (
          all(
            is.na(rainfall_mm)
          )
        ) {
          NA_real_
        } else {
          max(
            rainfall_mm,
            na.rm = TRUE
          )
        },
      
      
      # ------------------------------------------------------
      # Observation quality flags
      # ------------------------------------------------------
      
      n_quality_flag_missing =
        sum(
          is.na(
            observation_quality_flag
          )
        ),
      
      n_quality_flag_zero =
        sum(
          observation_quality_flag == 0,
          na.rm = TRUE
        ),
      
      n_quality_flag_nonzero =
        sum(
          observation_quality_flag != 0,
          na.rm = TRUE
        ),
      
      
      # ------------------------------------------------------
      # Rainfall-period diagnostics
      # ------------------------------------------------------
      
      n_period_missing =
        sum(
          is.na(
            rainfall_period
          )
        ),
      
      n_period_one =
        sum(
          rainfall_period == 1,
          na.rm = TRUE
        ),
      
      n_nonstandard_period =
        sum(
          !is.na(rainfall_period) &
            rainfall_period != 1
        ),
      
      
      # ------------------------------------------------------
      # Multi-day diagnostics
      # ------------------------------------------------------
      
      n_day_count_missing =
        sum(
          is.na(
            rainfall_day_count
          )
        ),
      
      n_day_count_one =
        sum(
          rainfall_day_count == 1,
          na.rm = TRUE
        ),
      
      n_multiday =
        sum(
          !is.na(rainfall_day_count) &
            rainfall_day_count > 1
        ),
      
      
      # ------------------------------------------------------
      # AWS representation
      # ------------------------------------------------------
      
      n_aws =
        sum(
          aws_flag == 1,
          na.rm = TRUE
        ),
      
      n_non_aws =
        sum(
          aws_flag == 0,
          na.rm = TRUE
        )
      
    ),
    by = .(
      station_number,
      year
    )
  ]
  
  
  # ==========================================================
  # Calendar completeness
  # ==========================================================
  
  station_year[
    ,
    expected_days :=
      days_in_year(year)
  ]
  
  
  station_year[
    ,
    completeness :=
      n_observed /
      expected_days
  ]
  
  
  # ----------------------------------------------------------
  # Row coverage measures
  #
  # row_coverage says whether PRCP_DAYS contains a database
  # row for each calendar day.
  #
  # completeness says whether rainfall itself is observed.
  # ----------------------------------------------------------
  
  station_year[
    ,
    row_coverage :=
      n_rows /
      expected_days
  ]
  
  
  # ----------------------------------------------------------
  # Wet-day fraction
  #
  # Denominator is observed rainfall days only.
  # ----------------------------------------------------------
  
  station_year[
    ,
    wet_day_fraction :=
      data.table::fifelse(
        n_observed > 0,
        n_wet / n_observed,
        NA_real_
      )
  ]
  
  
  # ----------------------------------------------------------
  # AWS fraction
  #
  # Based on database rows for which AWS_FLAG is available.
  # ----------------------------------------------------------
  
  station_year[
    ,
    aws_fraction :=
      data.table::fifelse(
        
        (n_aws + n_non_aws) > 0,
        
        n_aws /
          (n_aws + n_non_aws),
        
        NA_real_
      )
  ]
  
  
  # ----------------------------------------------------------
  # Provisional complete-year indicators
  #
  # These are diagnostic thresholds only.
  # No final scientific threshold has yet been selected.
  # ----------------------------------------------------------
  
  station_year[
    ,
    `:=`(
      
      complete_80 =
        completeness >= 0.80,
      
      complete_90 =
        completeness >= 0.90,
      
      complete_95 =
        completeness >= 0.95,
      
      complete_99 =
        completeness >= 0.99
    )
  ]
  
  
  # ----------------------------------------------------------
  # Source provenance
  # ----------------------------------------------------------
  
  station_year[
    ,
    source_parquet :=
      basename(path)
  ]
  
  
  data.table::setorder(
    station_year,
    station_number,
    year
  )
  
  
  # ----------------------------------------------------------
  # Release daily data
  # ----------------------------------------------------------
  
  rm(rainfall)
  
  invisible(
    gc()
  )
  
  
  station_year
}


# ============================================================
# Process every rainfall Parquet file
# ============================================================

summarise_all_rainfall_availability <- function(
    parquet_files) {
  
  if (length(parquet_files) == 0L) {
    
    stop(
      "No rainfall Parquet files supplied.",
      call. = FALSE
    )
  }
  
  
  missing_files <- parquet_files[
    !file.exists(
      parquet_files
    )
  ]
  
  
  if (length(missing_files) > 0L) {
    
    stop(
      "Missing rainfall Parquet files:\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  
  parquet_files <- sort(
    parquet_files
  )
  
  
  results <- vector(
    mode = "list",
    length = length(parquet_files)
  )
  
  
  for (i in seq_along(parquet_files)) {
    
    message(
      "\n[",
      i,
      "/",
      length(parquet_files),
      "] ",
      basename(
        parquet_files[[i]]
      )
    )
    
    
    results[[i]] <-
      summarise_rainfall_parquet_file(
        parquet_files[[i]]
      )
  }
  
  
  result <- data.table::rbindlist(
    results,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  # ----------------------------------------------------------
  # No station-year should appear in multiple source chunks
  # because extraction periods are non-overlapping.
  # ----------------------------------------------------------
  
  duplicated_station_year <- result[
    ,
    .N,
    by = .(
      station_number,
      year
    )
  ][
    N > 1L
  ]
  
  
  if (
    nrow(
      duplicated_station_year
    ) > 0L
  ) {
    
    stop(
      "Station-year summaries overlap across rainfall ",
      "Parquet chunks. Found ",
      nrow(
        duplicated_station_year
      ),
      " duplicated station-year keys.",
      call. = FALSE
    )
  }
  
  
  data.table::setorder(
    result,
    station_number,
    year
  )
  
  
  result
}


# ============================================================
# National network summary by year
# ============================================================

summarise_rainfall_network_by_year <- function(
    station_year) {
  
  required <- c(
    "station_number",
    "year",
    "n_observed",
    "completeness",
    "n_multiday",
    "n_nonstandard_period"
  )
  
  
  missing_columns <- setdiff(
    required,
    names(
      station_year
    )
  )
  
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "Station-year dataset is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  network <- station_year[
    ,
    .(
      
      stations_with_rows =
        uniqueN(
          station_number
        ),
      
      stations_with_observations =
        sum(
          n_observed > 0
        ),
      
      stations_80pct =
        sum(
          completeness >= 0.80,
          na.rm = TRUE
        ),
      
      stations_90pct =
        sum(
          completeness >= 0.90,
          na.rm = TRUE
        ),
      
      stations_95pct =
        sum(
          completeness >= 0.95,
          na.rm = TRUE
        ),
      
      stations_99pct =
        sum(
          completeness >= 0.99,
          na.rm = TRUE
        ),
      
      median_completeness =
        stats::median(
          completeness,
          na.rm = TRUE
        ),
      
      mean_completeness =
        mean(
          completeness,
          na.rm = TRUE
        ),
      
      station_years_with_multiday =
        sum(
          n_multiday > 0
        ),
      
      multiday_records =
        sum(
          n_multiday
        ),
      
      station_years_with_nonstandard_period =
        sum(
          n_nonstandard_period > 0
        ),
      
      nonstandard_period_records =
        sum(
          n_nonstandard_period
        )
    ),
    by = year
  ]
  
  
  data.table::setorder(
    network,
    year
  )
  
  
  network
}


# ============================================================
# Station-level record summary
# ============================================================

summarise_rainfall_station_records <- function(
    station_year) {
  
  
  station_summary <- station_year[
    ,
    .(
      
      first_year =
        min(year),
      
      last_year =
        max(year),
      
      years_with_rows =
        .N,
      
      years_with_observations =
        sum(
          n_observed > 0
        ),
      
      years_80pct =
        sum(
          completeness >= 0.80,
          na.rm = TRUE
        ),
      
      years_90pct =
        sum(
          completeness >= 0.90,
          na.rm = TRUE
        ),
      
      years_95pct =
        sum(
          completeness >= 0.95,
          na.rm = TRUE
        ),
      
      years_99pct =
        sum(
          completeness >= 0.99,
          na.rm = TRUE
        ),
      
      mean_completeness =
        mean(
          completeness,
          na.rm = TRUE
        ),
      
      median_completeness =
        stats::median(
          completeness,
          na.rm = TRUE
        ),
      
      total_observed_days =
        sum(
          n_observed
        ),
      
      total_missing_days =
        sum(
          n_missing
        ),
      
      total_wet_days =
        sum(
          n_wet
        ),
      
      total_multiday_records =
        sum(
          n_multiday
        ),
      
      total_nonstandard_period_records =
        sum(
          n_nonstandard_period
        ),
      
      maximum_recorded_rainfall_mm =
        if (
          all(
            is.na(
              maximum_recorded_rainfall_mm
            )
          )
        ) {
          NA_real_
        } else {
          max(
            maximum_recorded_rainfall_mm,
            na.rm = TRUE
          )
        }
    ),
    by = station_number
  ]
  
  
  station_summary[
    ,
    record_span_years :=
      last_year -
      first_year +
      1L
  ]
  
  
  data.table::setorder(
    station_summary,
    station_number
  )
  
  
  station_summary
}


# ============================================================
# Candidate-network table
#
# This does NOT choose the final network.
# It reports how many stations survive useful combinations of
# year-completeness and record-length criteria.
# ============================================================

summarise_candidate_rainfall_networks <- function(
    station_year,
    completeness_thresholds =
      c(
        0.80,
        0.90,
        0.95,
        0.99
      ),
    minimum_years =
      c(
        20L,
        30L,
        40L,
        50L,
        60L,
        70L,
        80L,
        100L
      )) {
  
  
  combinations <- data.table::CJ(
    
    completeness_threshold =
      completeness_thresholds,
    
    minimum_years =
      minimum_years
  )
  
  
  result <- combinations[
    ,
    {
      
      threshold <-
        completeness_threshold
      
      qualifying <- station_year[
        completeness >= threshold
      ]
      
      
      qualifying_station <- qualifying[
        ,
        .(
          qualifying_years =
            uniqueN(year)
        ),
        by = station_number
      ]
      
      
      .(
        stations =
          sum(
            qualifying_station$
              qualifying_years >=
              minimum_years
          )
      )
      
    },
    by = .(
      completeness_threshold,
      minimum_years
    )
  ]
  
  
  result
}