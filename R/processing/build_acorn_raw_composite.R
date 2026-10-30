# ============================================================
# build_acorn_raw_composite.R
#
# Resolve duplicate raw observations and construct composite
# ACORN-location records using primarysites.txt.
#
# Policy:
#   1. Exact and complementary duplicates are consolidated.
#   2. Station-dates containing conflicting non-missing
#      temperatures are excluded completely.
#   3. All exclusions are retained in an audit table.
# ============================================================


resolve_acorn_duplicates_for_composite <- function(raw_data) {
  required <- c(
    "station_number",
    "date",
    "tmax",
    "tmin",
    "tmax_missing",
    "tmin_missing",
    "source_dataset",
    "source_version",
    "data_version",
    "source_member"
  )
  
  missing_columns <- setdiff(required, names(raw_data))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Raw data are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  x <- data.table::copy(raw_data)
  
  # Identify every row associated with a duplicated key.
  duplicate_indicator <-
    duplicated(x, by = c("station_number", "date")) |
    duplicated(x,
               by = c("station_number", "date"),
               fromLast = TRUE)
  
  unique_records <- x[!duplicate_indicator]
  
  duplicate_records <- x[duplicate_indicator]
  
  if (nrow(duplicate_records) == 0L) {
    return(
      list(
        data = x,
        conflict_audit = data.table::data.table(),
        duplicate_summary = data.table::data.table()
      )
    )
  }
  
  duplicate_summary <- duplicate_records[, .(
    n_source_rows = .N,
    
    n_tmax_nonmissing_values =
      data.table::uniqueN(tmax[!is.na(tmax)]),
    
    n_tmin_nonmissing_values =
      data.table::uniqueN(tmin[!is.na(tmin)])
  ), by = .(station_number, date)]
  
  conflict_keys <- duplicate_summary[n_tmax_nonmissing_values > 1L |
                                       n_tmin_nonmissing_values > 1L, .(station_number, date)]
  
  safe_keys <- duplicate_summary[n_tmax_nonmissing_values <= 1L &
                                   n_tmin_nonmissing_values <= 1L, .(station_number, date)]
  
  # Preserve every conflicting source row in the audit.
  conflict_audit <- duplicate_records[conflict_keys, on = .(station_number, date), nomatch = 0L, allow.cartesian = TRUE]
  
  if (nrow(conflict_audit) > 0L) {
    conflict_audit[, exclusion_reason :=
                     "conflicting_nonmissing_duplicate_values"]
  }
  
  # Retain only duplicates that can be safely consolidated.
  safe_records <- duplicate_records[safe_keys, on = .(station_number, date), nomatch = 0L, allow.cartesian = TRUE]
  
  resolved_safe <- safe_records[, {
    available_tmax <- tmax[!is.na(tmax)]
    
    available_tmin <- tmin[!is.na(tmin)]
    
    selected_tmax <- if (length(available_tmax) == 0L) {
      NA_real_
    } else {
      available_tmax[[1L]]
    }
    
    selected_tmin <- if (length(available_tmin) == 0L) {
      NA_real_
    } else {
      available_tmin[[1L]]
    }
    
    representative <- .SD[1L]
    
    representative[, `:=`(
      tmax = selected_tmax,
      tmin = selected_tmin,
      tmax_missing = is.na(selected_tmax),
      tmin_missing = is.na(selected_tmin)
    )]
    
    representative
  }, by = .(station_number, date)]
  
  # Conflicting keys are deliberately not included here.
  result <- data.table::rbindlist(list(unique_records, resolved_safe),
                                  use.names = TRUE,
                                  fill = TRUE)
  
  data.table::setorder(result, station_number, date)
  
  remaining_duplicates <- result[, .N, by = .(station_number, date)][N > 1L]
  
  if (nrow(remaining_duplicates) > 0L) {
    stop("Duplicate resolution failed: duplicate keys remain.",
         call. = FALSE)
  }
  
  message("Resolved ", nrow(safe_keys), " safe duplicate station-dates.")
  
  message(
    "Excluded ",
    nrow(conflict_keys),
    " conflicting station-dates represented by ",
    nrow(conflict_audit),
    " source rows."
  )
  
  list(
    data = result,
    conflict_audit = conflict_audit,
    duplicate_summary = duplicate_summary
  )
}



build_acorn_raw_composite <- function(raw_data, primary_sites_long) {
  raw_required <- c(
    "station_number",
    "date",
    "tmax",
    "tmin",
    "tmax_missing",
    "tmin_missing",
    "source_version",
    "source_member"
  )
  
  mapping_required <- c(
    "acorn_station_number",
    "source_station_number",
    "primary_site_order",
    "start_date",
    "end_date"
  )
  
  missing_raw <- setdiff(raw_required, names(raw_data))
  
  missing_mapping <- setdiff(mapping_required, names(primary_sites_long))
  
  if (length(missing_raw) > 0L) {
    stop("Raw data are missing columns: ",
         paste(missing_raw, collapse = ", "),
         call. = FALSE)
  }
  
  if (length(missing_mapping) > 0L) {
    stop(
      "Primary-site mapping is missing columns: ",
      paste(missing_mapping, collapse = ", "),
      call. = FALSE
    )
  }
  
  raw <- data.table::copy(raw_data)
  mapping <- data.table::copy(primary_sites_long)
  
  data.table::setkey(raw, station_number, date)
  
  segments <- vector("list", nrow(mapping))
  
  for (i in seq_len(nrow(mapping))) {
    acorn_id <-
      mapping$acorn_station_number[[i]]
    
    source_id <-
      mapping$source_station_number[[i]]
    
    site_order <-
      mapping$primary_site_order[[i]]
    
    period_start <-
      mapping$start_date[[i]]
    
    period_end <-
      mapping$end_date[[i]]
    
    segment <- raw[station_number == source_id &
                     date >= period_start &
                     date <= period_end]
    
    if (nrow(segment) == 0L) {
      warning(
        "No observations found for ACORN station ",
        acorn_id,
        ", source station ",
        source_id,
        ", period ",
        period_start,
        " to ",
        period_end,
        ".",
        call. = FALSE
      )
    }
    
    segment[, `:=`(
      acorn_station_number = acorn_id,
      source_station_number = station_number,
      primary_site_order = site_order,
      primary_period_start = period_start,
      primary_period_end = period_end,
      source_dataset = "acorn_sat_raw_composite",
      data_version =
        "quality_controlled_unadjusted_composite"
    )]
    
    segments[[i]] <- segment[, .(
      acorn_station_number,
      source_station_number,
      primary_site_order,
      primary_period_start,
      primary_period_end,
      date,
      tmax,
      tmin,
      tmax_missing,
      tmin_missing,
      source_dataset,
      source_version,
      data_version,
      source_member
    )]
  }
  
  result <- data.table::rbindlist(segments, use.names = TRUE, fill = TRUE)
  
  data.table::setorder(result, acorn_station_number, date)
  
  duplicate_dates <- result[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  if (nrow(duplicate_dates) > 0L) {
    stop(
      "Composite construction produced ",
      nrow(duplicate_dates),
      " duplicated ACORN station-dates.",
      call. = FALSE
    )
  }
  
  result
}



validate_acorn_raw_composite <- function(composite, primary_sites_long) {
  mapping_range <- primary_sites_long[, .(
    expected_first_date = min(start_date),
    expected_last_date = max(end_date)
  ), by = acorn_station_number]
  
  observed_range <- composite[, .(
    observed_first_date = min(date),
    observed_last_date = max(date),
    rows = .N,
    n_tmax = sum(!is.na(tmax)),
    n_tmin = sum(!is.na(tmin))
  ), by = acorn_station_number]
  
  audit <- merge(mapping_range,
                 observed_range,
                 by = "acorn_station_number",
                 all.x = TRUE)
  
  audit[, `:=`(
    first_date_matches =
      observed_first_date == expected_first_date,
    
    last_date_matches =
      observed_last_date == expected_last_date
  )]
  
  audit
}