# ============================================================
# build_rainfall_pot_events.R
#
# Construct peaks-over-threshold rainfall event datasets for
# the frozen continental Australian rainfall network.
#
# Default configuration:
#
#   Analysis period:       1950-2025
#   Reference period:      1961-1990
#   Wet-day threshold:     1 mm
#   POT thresholds:        P95, P97.5 and P99
#   Primary run length:    3 non-exceedance days
#   Sensitivity runs:      2, 3 and 5 days
#
# Event definition:
#
#   Two exceedances belong to the same event if fewer than
#   run_length complete non-exceedance days occur between them.
#
# Example for run_length = 3:
#
#   exceedance Monday
#   non-exceedance Tuesday
#   non-exceedance Wednesday
#   exceedance Thursday
#
# belongs to one event.
#
# But:
#
#   exceedance Monday
#   non-exceedance Tuesday
#   non-exceedance Wednesday
#   non-exceedance Thursday
#   exceedance Friday
#
# begins a new event.
#
# Daily one-day eligibility:
#
#   rainfall is observed
#   AND rainfall_period is missing or equals 1
#   AND rainfall_day_count is missing or <= 1
#
# Multi-day and non-standard-period observations are excluded
# from threshold construction and POT event extraction.
# ============================================================


# ============================================================
# Package checks
# ============================================================

check_rainfall_pot_packages <- function() {
  required <- c("data.table", "arrow")
  
  missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_packages) > 0L) {
    stop(
      "Required packages are not installed: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Parse rainfall Parquet filename periods
# ============================================================

rainfall_parquet_inventory <- function(parquet_files) {
  if (length(parquet_files) == 0L) {
    stop("No rainfall Parquet files supplied.", call. = FALSE)
  }
  
  missing_files <- parquet_files[!file.exists(parquet_files)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Rainfall Parquet files do not exist:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }
  
  filenames <- basename(parquet_files)
  
  match <- regexec("^adam_prcp_days_([0-9]{4})_([0-9]{4})[.]parquet$",
                   filenames)
  
  parts <- regmatches(filenames, match)
  
  valid <- lengths(parts) == 3L
  
  if (!all(valid)) {
    stop(
      "Could not determine date range from rainfall files: ",
      paste(filenames[!valid], collapse = ", "),
      call. = FALSE
    )
  }
  
  result <- data.table::data.table(
    path =
      parquet_files,
    
    filename =
      filenames,
    
    start_year =
      as.integer(vapply(parts, `[[`, character(1), 2L)),
    
    end_year =
      as.integer(vapply(parts, `[[`, character(1), 3L))
  )
  
  data.table::setorder(result, start_year, end_year)
  
  result
}


rainfall_parquet_files_for_period <- function(parquet_files,
                                              requested_start_year,
                                              requested_end_year) {
  requested_start_year <- as.integer(requested_start_year)
  
  requested_end_year <- as.integer(requested_end_year)
  
  if (length(requested_start_year) != 1L ||
      length(requested_end_year) != 1L ||
      is.na(requested_start_year) ||
      is.na(requested_end_year)) {
    stop(
      "requested_start_year and requested_end_year ",
      "must each be one valid integer.",
      call. = FALSE
    )
  }
  
  if (requested_end_year <
      requested_start_year) {
    stop("requested_end_year cannot be earlier than ",
         "requested_start_year.",
         call. = FALSE)
  }
  
  inventory <- rainfall_parquet_inventory(parquet_files)
  
  keep <-
    inventory[["end_year"]] >=
    requested_start_year &
    inventory[["start_year"]] <=
    requested_end_year
  
  inventory[["path"]][keep]
}

# ============================================================
# Prepare core network
# ============================================================

prepare_pot_core_network <- function(core_network) {
  required <- c(
    "station_number",
    "station_name",
    "latitude",
    "longitude",
    "elevation_m",
    "state",
    "rainfall_district",
    "area_name"
  )
  
  missing_columns <- setdiff(required, names(core_network))
  
  if (length(missing_columns) > 0L) {
    stop("Core network is missing columns: ",
         paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  
  network <- data.table::copy(core_network)
  
  if ("continental_core" %in%
      names(network)) {
    network <- network[continental_core == TRUE]
  }
  
  duplicate_stations <- network[, .N, by = station_number][N > 1L]
  
  if (nrow(duplicate_stations) > 0L) {
    stop("Core network contains duplicated station numbers.", call. = FALSE)
  }
  
  if (any(is.na(network$station_number))) {
    stop("Core network contains missing station numbers.", call. = FALSE)
  }
  
  data.table::setorder(network, station_number)
  
  network
}


# ============================================================
# Read one daily rainfall chunk
# ============================================================

read_pot_rainfall_chunk <- function(path,
                                    station_numbers,
                                    start_year,
                                    end_year) {
  if (!file.exists(path)) {
    stop("Rainfall Parquet file does not exist: ", path, call. = FALSE)
  }
  
  rainfall <- data.table::as.data.table(arrow::read_parquet(
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
  ))
  
  rainfall <- rainfall[station_number %in%
                         station_numbers]
  
  rainfall[, year :=
             as.integer(format(date, "%Y"))]
  
  rainfall <- rainfall[year >= start_year &
                         year <= end_year]
  
  rainfall
}


# ============================================================
# One-day eligibility
# ============================================================

add_pot_eligibility <- function(rainfall, wet_day_threshold_mm = 1) {
  required <- c("station_number",
                "date",
                "rainfall_mm",
                "rainfall_period",
                "rainfall_day_count")
  
  missing_columns <- setdiff(required, names(rainfall))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Rainfall data are missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  result <- data.table::copy(rainfall)
  
  result[, rainfall_observed :=
           !is.na(rainfall_mm)]
  
  result[, multiday_record :=
           rainfall_observed &
           !is.na(rainfall_day_count) &
           rainfall_day_count > 1]
  
  result[, nonstandard_period :=
           rainfall_observed &
           !is.na(rainfall_period) &
           rainfall_period != 1]
  
  result[, eligible_1day :=
           rainfall_observed &
           (is.na(rainfall_period) |
              rainfall_period == 1) &
           (is.na(rainfall_day_count) |
              rainfall_day_count <= 1)]
  
  result[, eligible_wet_day :=
           eligible_1day &
           rainfall_mm >=
           wet_day_threshold_mm]
  
  result
}


# ============================================================
# Build station-specific reference thresholds
# ============================================================

build_pot_reference_thresholds <- function(parquet_files,
                                           core_network,
                                           reference_start_year = 1961L,
                                           reference_end_year = 1990L,
                                           threshold_probabilities = c(0.95, 0.975, 0.99),
                                           wet_day_threshold_mm = 1,
                                           minimum_reference_wet_days = 300L,
                                           quantile_type = 8L) {
  check_rainfall_pot_packages()
  
  network <- prepare_pot_core_network(core_network)
  
  station_numbers <- network$station_number
  
  reference_files <-
    rainfall_parquet_files_for_period(
      parquet_files = parquet_files,
      requested_start_year = reference_start_year,
      requested_end_year = reference_end_year
    )
  
  if (length(reference_files) == 0L) {
    stop("No rainfall Parquet files overlap the reference period.",
         call. = FALSE)
  }
  
  reference_tables <- vector(mode = "list", length = length(reference_files))
  
  for (i in seq_along(reference_files)) {
    message(
      "\nReference threshold file ",
      i,
      " of ",
      length(reference_files),
      ": ",
      basename(reference_files[[i]])
    )
    
    rainfall <- read_pot_rainfall_chunk(
      path =
        reference_files[[i]],
      
      station_numbers =
        station_numbers,
      
      start_year =
        reference_start_year,
      
      end_year =
        reference_end_year
    )
    
    rainfall <- add_pot_eligibility(rainfall =
                                      rainfall, wet_day_threshold_mm =
                                      wet_day_threshold_mm)
    
    reference_tables[[i]] <- rainfall[eligible_wet_day == TRUE, .(station_number, rainfall_mm)]
    
    rm(rainfall)
    
    invisible(gc())
  }
  
  reference_wet_days <- data.table::rbindlist(reference_tables, use.names = TRUE, fill = TRUE)
  
  if (nrow(reference_wet_days) == 0L) {
    stop("No eligible reference-period wet days were found.", call. = FALSE)
  }
  
  threshold_wide <- reference_wet_days[, {
    quantiles <- stats::quantile(
      rainfall_mm,
      probs =
        threshold_probabilities,
      na.rm = TRUE,
      type = quantile_type,
      names = FALSE
    )
    
    names(quantiles) <- paste0("q_", gsub(
      "[.]",
      "_",
      format(
        threshold_probabilities,
        trim = TRUE,
        scientific = FALSE
      )
    ))
    
    as.list(quantiles)
  }, by = station_number]
  
  reference_counts <- reference_wet_days[, .(n_reference_wet_days = .N), by = station_number]
  
  threshold_wide <- merge(
    threshold_wide,
    reference_counts,
    by = "station_number",
    all = TRUE,
    sort = FALSE
  )
  
  threshold_long <- data.table::melt(
    threshold_wide,
    id.vars = c("station_number", "n_reference_wet_days"),
    variable.name =
      "threshold_column",
    
    value.name =
      "threshold_mm"
  )
  
  probability_lookup <- data.table::data.table(threshold_column =
                                                 paste0("q_", gsub(
                                                   "[.]",
                                                   "_",
                                                   format(
                                                     threshold_probabilities,
                                                     trim = TRUE,
                                                     scientific = FALSE
                                                   )
                                                 )),
                                               
                                               threshold_probability =
                                                 threshold_probabilities)
  
  threshold_long <- merge(
    threshold_long,
    probability_lookup,
    by = "threshold_column",
    all.x = TRUE,
    sort = FALSE
  )
  
  threshold_long[, threshold_label :=
                   paste0("P",
                          format(
                            100 *
                              threshold_probability,
                            trim = TRUE,
                            scientific = FALSE
                          ))]
  
  threshold_long[, threshold_id :=
                   data.table::fcase(
                     abs(threshold_probability -
                           0.95) < 1e-10,
                     "p95",
                     
                     abs(threshold_probability -
                           0.975) < 1e-10,
                     "p975",
                     
                     abs(threshold_probability -
                           0.99) < 1e-10,
                     "p99",
                     
                     default =
                       paste0("p", gsub(
                         "[.]",
                         "",
                         format(
                           100 *
                             threshold_probability,
                           trim = TRUE,
                           scientific = FALSE
                         )
                       ))
                   )]
  
  threshold_long[, reference_threshold_valid :=
                   n_reference_wet_days >=
                   minimum_reference_wet_days]
  
  threshold_long[reference_threshold_valid == FALSE, threshold_mm := NA_real_]
  
  threshold_long[, `:=`(
    reference_start_year =
      as.integer(reference_start_year),
    
    reference_end_year =
      as.integer(reference_end_year),
    
    wet_day_threshold_mm =
      wet_day_threshold_mm,
    
    minimum_reference_wet_days =
      as.integer(minimum_reference_wet_days),
    
    quantile_type =
      as.integer(quantile_type)
  )]
  
  threshold_long <- merge(
    network[, .(
      station_number,
      station_name,
      latitude,
      longitude,
      elevation_m,
      state,
      rainfall_district,
      area_name
    )],
    threshold_long,
    by = "station_number",
    all.x = TRUE,
    allow.cartesian = TRUE,
    sort = FALSE
  )
  
  data.table::setorder(threshold_long, threshold_probability, station_number)
  
  rm(reference_wet_days)
  rm(reference_tables)
  
  invisible(gc())
  
  threshold_long
}


# ============================================================
# Extract daily threshold exceedances
# ============================================================

extract_pot_exceedances <- function(parquet_files,
                                    core_network,
                                    thresholds,
                                    analysis_start_year = 1950L,
                                    analysis_end_year = 2025L,
                                    wet_day_threshold_mm = 1) {
  check_rainfall_pot_packages()
  
  network <- prepare_pot_core_network(core_network)
  
  station_numbers <- network$station_number
  
  analysis_files <-
    rainfall_parquet_files_for_period(
      parquet_files = parquet_files,
      requested_start_year = analysis_start_year,
      requested_end_year = analysis_end_year
    )
  
  threshold_keep <- thresholds[reference_threshold_valid == TRUE &
                                 !is.na(threshold_mm), .(
                                   station_number,
                                   threshold_id,
                                   threshold_label,
                                   threshold_probability,
                                   threshold_mm
                                 )]
  
  if (nrow(threshold_keep) == 0L) {
    stop("No valid POT thresholds are available.", call. = FALSE)
  }
  
  exceedance_tables <- vector(mode = "list", length = length(analysis_files))
  
  daily_audits <- vector(mode = "list", length = length(analysis_files))
  
  for (i in seq_along(analysis_files)) {
    message(
      "\nPOT extraction file ",
      i,
      " of ",
      length(analysis_files),
      ": ",
      basename(analysis_files[[i]])
    )
    
    rainfall <- read_pot_rainfall_chunk(
      path =
        analysis_files[[i]],
      
      station_numbers =
        station_numbers,
      
      start_year =
        analysis_start_year,
      
      end_year =
        analysis_end_year
    )
    
    rainfall <- add_pot_eligibility(rainfall =
                                      rainfall, wet_day_threshold_mm =
                                      wet_day_threshold_mm)
    
    daily_audits[[i]] <- rainfall[, .(
      rows = .N,
      
      observed =
        sum(rainfall_observed),
      
      eligible_1day =
        sum(eligible_1day),
      
      eligible_wet_days =
        sum(eligible_wet_day),
      
      multiday_excluded =
        sum(multiday_record),
      
      nonstandard_period_excluded =
        sum(nonstandard_period)
    ), by = year][, source_parquet :=
                    basename(analysis_files[[i]])]
    
    eligible <- rainfall[eligible_wet_day == TRUE, .(
      station_number,
      date,
      year,
      rainfall_mm,
      observation_quality_flag,
      rainfall_quality,
      aws_flag,
      source_file
    )]
    
    if (nrow(eligible) == 0L) {
      exceedance_tables[[i]] <-
        data.table::data.table()
      
      rm(rainfall)
      rm(eligible)
      
      invisible(gc())
      
      next
    }
    
    candidates <- merge(
      eligible,
      threshold_keep,
      by = "station_number",
      allow.cartesian = TRUE,
      sort = FALSE
    )
    
    candidates <- candidates[rainfall_mm >
                               threshold_mm]
    
    candidates[, excess_mm :=
                 rainfall_mm -
                 threshold_mm]
    
    candidates[, source_parquet :=
                 basename(analysis_files[[i]])]
    
    exceedance_tables[[i]] <- candidates
    
    rm(rainfall)
    rm(eligible)
    rm(candidates)
    
    invisible(gc())
  }
  
  exceedances <- data.table::rbindlist(exceedance_tables, use.names = TRUE, fill = TRUE)
  
  daily_audit <- data.table::rbindlist(daily_audits, use.names = TRUE, fill = TRUE)
  
  duplicate_exceedances <- exceedances[, .N, by = .(station_number, threshold_id, date)][N > 1L]
  
  if (nrow(duplicate_exceedances) > 0L) {
    stop(
      "Duplicate station-threshold-date exceedances found: ",
      nrow(duplicate_exceedances),
      call. = FALSE
    )
  }
  
  data.table::setorder(exceedances, threshold_probability, station_number, date)
  
  list(exceedances =
         exceedances, daily_audit =
         daily_audit)
}


# ============================================================
# Decluster exceedances
# ============================================================

decluster_pot_exceedances <- function(exceedances, run_length_days = 3L) {
  run_length_days <- as.integer(run_length_days)
  
  if (length(run_length_days) != 1L ||
      is.na(run_length_days) ||
      run_length_days < 1L) {
    stop("run_length_days must be one positive integer.", call. = FALSE)
  }
  
  required <- c(
    "station_number",
    "date",
    "rainfall_mm",
    "threshold_id",
    "threshold_label",
    "threshold_probability",
    "threshold_mm",
    "excess_mm"
  )
  
  missing_columns <- setdiff(required, names(exceedances))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Exceedance data are missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  x <- data.table::copy(exceedances)
  
  data.table::setorder(x, station_number, threshold_probability, date)
  
  x[, previous_exceedance_date :=
      data.table::shift(date), by = .(station_number, threshold_id)]
  
  x[, days_since_previous_exceedance :=
      as.integer(date -
                   previous_exceedance_date)]
  
  # New event if there are at least run_length_days complete
  # non-exceedance days between successive exceedance dates.
  #
  # If date difference = run_length_days + 1, then exactly
  # run_length_days dates lie between the exceedances.
  
  x[, new_event :=
      is.na(previous_exceedance_date) |
      days_since_previous_exceedance >
      run_length_days]
  
  x[, event_sequence :=
      cumsum(new_event), by = .(station_number, threshold_id)]
  
  event_summary <- x[, {
    peak_index <- which.max(rainfall_mm)
    
    list(
      event_start =
        min(date),
      
      event_end =
        max(date),
      
      event_peak_date =
        date[peak_index],
      
      event_peak_mm =
        rainfall_mm[peak_index],
      
      threshold_mm =
        data.table::first(threshold_mm),
      
      excess_mm =
        rainfall_mm[peak_index] -
        data.table::first(threshold_mm),
      
      event_exceedance_total_mm =
        sum(rainfall_mm, na.rm = TRUE),
      
      event_excess_total_mm =
        sum(excess_mm, na.rm = TRUE),
      
      n_exceedance_days =
        .N,
      
      maximum_gap_between_exceedances_days =
        if (.N <= 1L) {
          0L
        } else {
          as.integer(max(diff(date)))
        },
      
      observation_quality_flag_peak =
        observation_quality_flag[peak_index],
      
      rainfall_quality_peak =
        rainfall_quality[peak_index],
      
      aws_flag_peak =
        aws_flag[peak_index]
    )
  }, by = .(
    station_number,
    threshold_id,
    threshold_label,
    threshold_probability,
    event_sequence
  )]
  
  event_summary[, event_duration_days :=
                  as.integer(event_end -
                               event_start) +
                  1L]
  
  event_summary[, run_length_days :=
                  run_length_days]
  
  event_summary[, event_id :=
                  sprintf(
                    "%s_%s_r%02d_%06d",
                    station_number,
                    threshold_id,
                    run_length_days,
                    event_sequence
                  )]
  
  event_summary[, `:=`(year =
                         as.integer(format(event_peak_date, "%Y")), month =
                         as.integer(format(event_peak_date, "%m")))]
  
  event_summary[, season :=
                  data.table::fcase(
                    month %in% c(12L, 1L, 2L),
                    "summer",
                    
                    month %in% c(3L, 4L, 5L),
                    "autumn",
                    
                    month %in% c(6L, 7L, 8L),
                    "winter",
                    
                    month %in% c(9L, 10L, 11L),
                    "spring",
                    
                    default =
                      NA_character_
                  )]
  
  event_summary[, water_year :=
                  data.table::fifelse(month >= 7L, year + 1L, year)]
  
  data.table::setorder(
    event_summary,
    threshold_probability,
    run_length_days,
    station_number,
    event_peak_date
  )
  
  event_summary
}


# ============================================================
# Add station metadata
# ============================================================

add_pot_station_metadata <- function(events, core_network) {
  network <- prepare_pot_core_network(core_network)
  
  metadata <- network[, .(
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
    close_date
  )]
  
  result <- merge(events,
                  metadata,
                  by = "station_number",
                  all.x = TRUE,
                  sort = FALSE)
  
  data.table::setorder(result,
                       threshold_probability,
                       run_length_days,
                       station_number,
                       event_peak_date)
  
  result
}


# ============================================================
# Add total eligible rainfall over each event span
#
# This makes a second sequential pass through the daily archive.
# It sums all eligible one-day rainfall between event_start and
# event_end, including eligible non-exceedance days between peaks.
# ============================================================

add_pot_event_span_totals <- function(events,
                                      parquet_files,
                                      core_network,
                                      analysis_start_year = 1950L,
                                      analysis_end_year = 2025L,
                                      wet_day_threshold_mm = 1) {
  if (nrow(events) == 0L) {
    events[, event_span_total_mm :=
             numeric()]
    
    return(events)
  }
  
  network <- prepare_pot_core_network(core_network)
  
  station_numbers <- network$station_number
  
  analysis_files <-
    rainfall_parquet_files_for_period(
      parquet_files = parquet_files,
      requested_start_year = analysis_start_year,
      requested_end_year = analysis_end_year
    )
  
  event_intervals <- events[, .(station_number, event_id, event_start, event_end)]
  
  event_intervals[, `:=`(interval_start =
                           event_start, interval_end =
                           event_end)]
  
  data.table::setkey(event_intervals,
                     station_number,
                     interval_start,
                     interval_end)
  
  span_tables <- vector(mode = "list", length = length(analysis_files))
  
  for (i in seq_along(analysis_files)) {
    message(
      "\nEvent-span rainfall file ",
      i,
      " of ",
      length(analysis_files),
      ": ",
      basename(analysis_files[[i]])
    )
    
    rainfall <- read_pot_rainfall_chunk(
      path =
        analysis_files[[i]],
      
      station_numbers =
        station_numbers,
      
      start_year =
        analysis_start_year,
      
      end_year =
        analysis_end_year
    )
    
    rainfall <- add_pot_eligibility(rainfall =
                                      rainfall, wet_day_threshold_mm =
                                      wet_day_threshold_mm)
    
    daily <- rainfall[eligible_1day == TRUE, .(station_number, date, rainfall_mm)]
    
    if (nrow(daily) == 0L) {
      span_tables[[i]] <-
        data.table::data.table()
      
      rm(rainfall)
      rm(daily)
      
      invisible(gc())
      
      next
    }
    
    daily[, `:=`(interval_start =
                   date, interval_end =
                   date)]
    
    data.table::setkey(daily, station_number, interval_start, interval_end)
    
    matched <- data.table::foverlaps(
      daily,
      event_intervals,
      by.x = c("station_number", "interval_start", "interval_end"),
      by.y = c("station_number", "interval_start", "interval_end"),
      type = "within",
      nomatch = 0L
    )
    
    span_tables[[i]] <- matched[, .(
      event_span_total_partial_mm =
        sum(rainfall_mm, na.rm = TRUE),
      
      event_span_observed_days_partial =
        .N
    ), by = event_id]
    
    rm(rainfall)
    rm(daily)
    rm(matched)
    
    invisible(gc())
  }
  
  span_totals <- data.table::rbindlist(span_tables, use.names = TRUE, fill = TRUE)
  
  span_totals <- span_totals[, .(
    event_span_total_mm =
      sum(event_span_total_partial_mm, na.rm = TRUE),
    
    event_span_observed_days =
      sum(event_span_observed_days_partial, na.rm = TRUE)
  ), by = event_id]
  
  result <- merge(events,
                  span_totals,
                  by = "event_id",
                  all.x = TRUE,
                  sort = FALSE)
  
  result[is.na(event_span_total_mm), event_span_total_mm := 0]
  
  result[is.na(event_span_observed_days), event_span_observed_days := 0L]
  
  data.table::setorder(result,
                       threshold_probability,
                       run_length_days,
                       station_number,
                       event_peak_date)
  
  result
}


# ============================================================
# Station-level POT summaries
# ============================================================

summarise_pot_events_by_station <- function(
    events,
    analysis_start_year = 1950L,
    analysis_end_year = 2025L) {
  
  analysis_years <-
    as.integer(
      analysis_end_year -
        analysis_start_year +
        1L
    )
  
  result <- events[
    ,
    .(
      events =
        as.integer(.N),
      
      events_per_year =
        as.numeric(
          .N / analysis_years
        ),
      
      first_event_date =
        min(
          event_peak_date
        ),
      
      last_event_date =
        max(
          event_peak_date
        ),
      
      median_peak_mm =
        as.numeric(
          stats::median(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      maximum_peak_mm =
        as.numeric(
          max(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      mean_excess_mm =
        as.numeric(
          mean(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      median_excess_mm =
        as.numeric(
          stats::median(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      maximum_excess_mm =
        as.numeric(
          max(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      mean_event_duration_days =
        as.numeric(
          mean(
            event_duration_days,
            na.rm = TRUE
          )
        ),
      
      maximum_event_duration_days =
        as.integer(
          max(
            event_duration_days,
            na.rm = TRUE
          )
        ),
      
      multiday_exceedance_events =
        as.integer(
          sum(
            n_exceedance_days > 1L
          )
        ),
      
      summer_events =
        as.integer(
          sum(
            season == "summer"
          )
        ),
      
      autumn_events =
        as.integer(
          sum(
            season == "autumn"
          )
        ),
      
      winter_events =
        as.integer(
          sum(
            season == "winter"
          )
        ),
      
      spring_events =
        as.integer(
          sum(
            season == "spring"
          )
        )
    ),
    by = .(
      station_number,
      station_name,
      latitude,
      longitude,
      elevation_m,
      state,
      rainfall_district,
      threshold_id,
      threshold_label,
      threshold_probability,
      threshold_mm,
      run_length_days
    )
  ]
  
  data.table::setorder(
    result,
    threshold_probability,
    run_length_days,
    station_number
  )
  
  result
}

# ============================================================
# Annual POT event summaries
# ============================================================

summarise_pot_events_by_year <- function(
    events) {
  
  result <- events[
    ,
    .(
      events =
        as.integer(.N),
      
      stations_with_events =
        as.integer(
          data.table::uniqueN(
            station_number
          )
        ),
      
      median_peak_mm =
        as.numeric(
          stats::median(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      maximum_peak_mm =
        as.numeric(
          max(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      median_excess_mm =
        as.numeric(
          stats::median(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      mean_excess_mm =
        as.numeric(
          mean(
            excess_mm,
            na.rm = TRUE
          )
        )
    ),
    by = .(
      threshold_id,
      threshold_label,
      threshold_probability,
      run_length_days,
      year
    )
  ]
  
  data.table::setorder(
    result,
    threshold_probability,
    run_length_days,
    year
  )
  
  result
}

# ============================================================
# Declustering sensitivity
# ============================================================

summarise_declustering_sensitivity <- function(events) {
  
  station_counts <- events[
    ,
    .(
      station_events =
        as.integer(.N)
    ),
    by = .(
      threshold_id,
      threshold_label,
      threshold_probability,
      run_length_days,
      station_number
    )
  ]
  
  event_summary <- events[
    ,
    .(
      events =
        as.integer(.N),
      
      stations_with_events =
        as.integer(
          data.table::uniqueN(
            station_number
          )
        ),
      
      mean_peak_mm =
        as.numeric(
          mean(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      median_peak_mm =
        as.numeric(
          stats::median(
            event_peak_mm,
            na.rm = TRUE
          )
        ),
      
      mean_excess_mm =
        as.numeric(
          mean(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      median_excess_mm =
        as.numeric(
          stats::median(
            excess_mm,
            na.rm = TRUE
          )
        ),
      
      median_duration_days =
        as.numeric(
          stats::median(
            event_duration_days,
            na.rm = TRUE
          )
        ),
      
      mean_duration_days =
        as.numeric(
          mean(
            event_duration_days,
            na.rm = TRUE
          )
        ),
      
      multi_exceedance_events =
        as.integer(
          sum(
            n_exceedance_days > 1L,
            na.rm = TRUE
          )
        ),
      
      multi_exceedance_event_fraction =
        as.numeric(
          mean(
            n_exceedance_days > 1L,
            na.rm = TRUE
          )
        )
    ),
    by = .(
      threshold_id,
      threshold_label,
      threshold_probability,
      run_length_days
    )
  ]
  
  station_summary <- station_counts[
    ,
    .(
      median_events_per_station =
        as.numeric(
          stats::median(
            station_events,
            na.rm = TRUE
          )
        ),
      
      mean_events_per_station =
        as.numeric(
          mean(
            station_events,
            na.rm = TRUE
          )
        ),
      
      q10_events_per_station =
        as.numeric(
          stats::quantile(
            station_events,
            probs = 0.10,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      q90_events_per_station =
        as.numeric(
          stats::quantile(
            station_events,
            probs = 0.90,
            na.rm = TRUE,
            names = FALSE
          )
        ),
      
      minimum_events_per_station =
        as.integer(
          min(
            station_events,
            na.rm = TRUE
          )
        ),
      
      maximum_events_per_station =
        as.integer(
          max(
            station_events,
            na.rm = TRUE
          )
        )
    ),
    by = .(
      threshold_id,
      threshold_label,
      threshold_probability,
      run_length_days
    )
  ]
  
  result <- merge(
    event_summary,
    station_summary,
    by = c(
      "threshold_id",
      "threshold_label",
      "threshold_probability",
      "run_length_days"
    ),
    all = TRUE,
    sort = FALSE
  )
  
  data.table::setorder(
    result,
    threshold_probability,
    run_length_days
  )
  
  result
}

# ============================================================
# Threshold diagnostics
# ============================================================

summarise_pot_thresholds <- function(thresholds, events) {
  threshold_summary <- thresholds[, .(
    stations =
      .N,
    
    valid_thresholds =
      sum(reference_threshold_valid, na.rm = TRUE),
    
    minimum_reference_wet_days =
      min(n_reference_wet_days, na.rm = TRUE),
    
    median_reference_wet_days =
      stats::median(n_reference_wet_days, na.rm = TRUE),
    
    median_threshold_mm =
      stats::median(threshold_mm, na.rm = TRUE),
    
    q10_threshold_mm =
      stats::quantile(threshold_mm, 0.10, na.rm = TRUE, names = FALSE),
    
    q90_threshold_mm =
      stats::quantile(threshold_mm, 0.90, na.rm = TRUE, names = FALSE)
  ), by = .(threshold_id, threshold_label, threshold_probability)]
  
  event_summary <- events[, .(
    declustered_events =
      .N,
    
    stations_with_events =
      data.table::uniqueN(station_number),
    
    median_excess_mm =
      stats::median(excess_mm, na.rm = TRUE),
    
    maximum_excess_mm =
      max(excess_mm, na.rm = TRUE)
  ), by = .(threshold_id,
            threshold_label,
            threshold_probability,
            run_length_days)]
  
  merge(
    threshold_summary,
    event_summary,
    by = c("threshold_id", "threshold_label", "threshold_probability"),
    all.x = TRUE,
    allow.cartesian = TRUE,
    sort = FALSE
  )[order(threshold_probability, run_length_days)]
}


# ============================================================
# Complete POT event build
# ============================================================

build_rainfall_pot_events <- function(parquet_files,
                                      core_network,
                                      thresholds = NULL,
                                      analysis_start_year = 1950L,
                                      analysis_end_year = 2025L,
                                      reference_start_year = 1961L,
                                      reference_end_year = 1990L,
                                      threshold_probabilities = c(0.95, 0.975, 0.99),
                                      run_lengths = c(2L, 3L, 5L),
                                      primary_run_length = 3L,
                                      wet_day_threshold_mm = 1,
                                      minimum_reference_wet_days = 300L,
                                      quantile_type = 8L,
                                      calculate_event_span_totals = TRUE) {
  check_rainfall_pot_packages()
  
  network <- prepare_pot_core_network(core_network)
  
  if (is.null(thresholds)) {
    thresholds <- build_pot_reference_thresholds(
      parquet_files =
        parquet_files,
      
      core_network =
        network,
      
      reference_start_year =
        reference_start_year,
      
      reference_end_year =
        reference_end_year,
      
      threshold_probabilities =
        threshold_probabilities,
      
      wet_day_threshold_mm =
        wet_day_threshold_mm,
      
      minimum_reference_wet_days =
        minimum_reference_wet_days,
      
      quantile_type =
        quantile_type
    )
  }
  
  extraction <- extract_pot_exceedances(
    parquet_files =
      parquet_files,
    
    core_network =
      network,
    
    thresholds =
      thresholds,
    
    analysis_start_year =
      analysis_start_year,
    
    analysis_end_year =
      analysis_end_year,
    
    wet_day_threshold_mm =
      wet_day_threshold_mm
  )
  
  exceedances <- extraction$exceedances
  
  event_tables <- lapply(as.integer(run_lengths), function(run_length) {
    message("\nDeclustering exceedances using ",
            run_length,
            "-day run length.")
    
    decluster_pot_exceedances(exceedances =
                                exceedances, run_length_days =
                                run_length)
  })
  
  events_all <- data.table::rbindlist(event_tables, use.names = TRUE, fill = TRUE)
  
  events_all <- add_pot_station_metadata(events =
                                           events_all, core_network =
                                           network)
  
  if (calculate_event_span_totals) {
    events_all <- add_pot_event_span_totals(
      events =
        events_all,
      
      parquet_files =
        parquet_files,
      
      core_network =
        network,
      
      analysis_start_year =
        analysis_start_year,
      
      analysis_end_year =
        analysis_end_year,
      
      wet_day_threshold_mm =
        wet_day_threshold_mm
    )
  }
  
  primary_events <- events_all[run_length_days ==
                                 as.integer(primary_run_length)]
  
  station_summary <-
    summarise_pot_events_by_station(
      events =
        events_all,
      
      analysis_start_year =
        analysis_start_year,
      
      analysis_end_year =
        analysis_end_year
    )
  
  annual_summary <-
    summarise_pot_events_by_year(events_all)
  
  declustering_sensitivity <-
    summarise_declustering_sensitivity(events_all)
  
  threshold_summary <-
    summarise_pot_thresholds(thresholds =
                               thresholds, events =
                               events_all)
  
  list(
    thresholds =
      thresholds,
    
    exceedances =
      exceedances,
    
    events_all =
      events_all,
    
    primary_events =
      primary_events,
    
    station_summary =
      station_summary,
    
    annual_summary =
      annual_summary,
    
    declustering_sensitivity =
      declustering_sensitivity,
    
    threshold_summary =
      threshold_summary,
    
    daily_audit =
      extraction$daily_audit
  )
}