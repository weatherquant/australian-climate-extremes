# ============================================================
# 01_validate_analysis_inputs.R
#
# Validate and freeze the matched ACORN-SAT v2.3 raw-versus-
# homogenised daily temperature analysis inputs.
#
# This script does NOT rebuild the raw composites or rematch the
# homogenised series. It validates the existing comparison product,
# creates one station-level audit row per ACORN location, and writes
# the paired daily and annual analysis datasets used downstream.
# ============================================================

library(data.table)
library(arrow)

# ============================================================
# Configuration
# ============================================================

comparison_file <- file.path(
  "data-derived", "acorn", "comparison",
  "acorn_raw_homogenised_comparison_v2.3.parquet"
)

primary_sites_file <- file.path(
  "data-derived", "acorn", "acorn_primary_sites_v2.3.parquet"
)

raw_composite_file <- file.path(
  "data-derived", "acorn", "acorn_raw_composite_v2.3.parquet"
)

homogenised_file <- file.path(
  "data-derived", "acorn", "homogenised",
  "acorn_homogenised_daily_v2.3.0.parquet"
)

output_dir <- file.path(
  "data-derived", "temperature_homogenisation"
)

table_output_dir <- file.path(
  "tables", "temperature_homogenisation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

expected_acorn_locations <- 112L
analysis_start_date <- as.Date("1910-01-01")
analysis_end_date <- as.Date("2021-12-31")
annual_completeness_threshold <- 0.90
adjustment_tolerance <- 1e-10

# ============================================================
# Input checks
# ============================================================

required_files <- c(
  comparison_file,
  primary_sites_file,
  raw_composite_file,
  homogenised_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

comparison <- as.data.table(read_parquet(comparison_file))
primary_sites <- as.data.table(read_parquet(primary_sites_file))

required_comparison_columns <- c(
  "acorn_station_number", "date", "source_station_number",
  "primary_site_order", "station_name",
  "tmax_raw", "tmax_hom", "tmin_raw", "tmin_hom",
  "tmax_adjustment", "tmin_adjustment",
  "matched_tmax", "matched_tmin",
  "tmax_match_status", "tmin_match_status"
)

missing_columns <- setdiff(required_comparison_columns, names(comparison))

if (length(missing_columns) > 0L) {
  stop(
    "Comparison product is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Restrict explicitly to the version's supported date range.
comparison <- comparison[
  date >= analysis_start_date &
    date <= analysis_end_date
]

# ============================================================
# Core key and location validation
# ============================================================

key_duplicates <- comparison[
  , .N,
  by = .(acorn_station_number, date)
][N > 1L]

if (nrow(key_duplicates) > 0L) {
  fwrite(
    key_duplicates,
    file.path(table_output_dir, "temperature_duplicate_station_dates.csv")
  )
  stop(
    "Duplicate ACORN location-date rows found in comparison product.",
    call. = FALSE
  )
}

location_count <- uniqueN(comparison$acorn_station_number)

if (location_count != expected_acorn_locations) {
  stop(
    "Expected ", expected_acorn_locations,
    " ACORN locations but found ", location_count, ".",
    call. = FALSE
  )
}

if (uniqueN(primary_sites$acorn_station_number) != expected_acorn_locations) {
  stop(
    "Primary-site table does not contain exactly ",
    expected_acorn_locations,
    " unique ACORN locations.",
    call. = FALSE
  )
}

# ============================================================
# Validate matched flags and adjustment arithmetic
# ============================================================

comparison[, expected_matched_tmax :=
             is.finite(tmax_raw) & is.finite(tmax_hom)
]

comparison[, expected_matched_tmin :=
             is.finite(tmin_raw) & is.finite(tmin_hom)
]

flag_mismatch_tmax <- comparison[matched_tmax != expected_matched_tmax]
flag_mismatch_tmin <- comparison[matched_tmin != expected_matched_tmin]

if (nrow(flag_mismatch_tmax) > 0L || nrow(flag_mismatch_tmin) > 0L) {
  fwrite(
    flag_mismatch_tmax,
    file.path(table_output_dir, "temperature_tmax_match_flag_mismatches.csv")
  )
  fwrite(
    flag_mismatch_tmin,
    file.path(table_output_dir, "temperature_tmin_match_flag_mismatches.csv")
  )
  stop("Matched-value flags are inconsistent with finite paired data.", call. = FALSE)
}

comparison[, calculated_tmax_adjustment := fifelse(
  matched_tmax,
  tmax_hom - tmax_raw,
  NA_real_
)]

comparison[, calculated_tmin_adjustment := fifelse(
  matched_tmin,
  tmin_hom - tmin_raw,
  NA_real_
)]

tmax_adjustment_error <- comparison[
  matched_tmax == TRUE,
  max(abs(tmax_adjustment - calculated_tmax_adjustment), na.rm = TRUE)
]

tmin_adjustment_error <- comparison[
  matched_tmin == TRUE,
  max(abs(tmin_adjustment - calculated_tmin_adjustment), na.rm = TRUE)
]

if (
  !is.finite(tmax_adjustment_error) ||
  !is.finite(tmin_adjustment_error) ||
  tmax_adjustment_error > adjustment_tolerance ||
  tmin_adjustment_error > adjustment_tolerance
) {
  stop(
    "Stored adjustment arithmetic failed validation. Maximum errors: Tmax=",
    signif(tmax_adjustment_error, 6),
    ", Tmin=", signif(tmin_adjustment_error, 6),
    call. = FALSE
  )
}

# ============================================================
# Validate primary-site membership by date
# ============================================================

primary_sites[
  ,
  `:=`(
    acorn_station_number =
      as.character(
        acorn_station_number
      ),
    
    source_station_number =
      as.character(
        source_station_number
      )
  )
]


# Each ACORN/source/order combination must identify one declared
# primary-site period.

duplicate_primary_site_keys <- primary_sites[
  ,
  .N,
  by = .(
    acorn_station_number,
    source_station_number,
    primary_site_order
  )
][
  N > 1L
]


if (nrow(duplicate_primary_site_keys) > 0L) {
  
  fwrite(
    duplicate_primary_site_keys,
    file.path(
      table_output_dir,
      "temperature_duplicate_primary_site_keys.csv"
    )
  )
  
  stop(
    "Duplicated ACORN/source/order keys found in the primary-site table.",
    call. = FALSE
  )
}


# Add the declared primary-site dates to matching comparison rows.
#
# Use an update join so the identifiers already present in comparison
# are preserved. Rows with no raw source assignment remain identifiable
# as no-source calendar rows.

comparison[
  ,
  `:=`(
    declared_source_start_date =
      as.Date(
        NA
      ),
    
    declared_source_end_date =
      as.Date(
        NA
      )
  )
]


comparison[
  primary_sites,
  on = .(
    acorn_station_number,
    source_station_number,
    primary_site_order
  ),
  `:=`(
    declared_source_start_date =
      i.start_date,
    
    declared_source_end_date =
      i.end_date
  )
]


# A source-period match is required when the comparison row carries
# a raw station assignment. Rows with no source assignment and no raw
# Tmax/Tmin value are retained as legitimate no-raw calendar rows.

comparison[
  ,
  raw_source_assigned :=
    !is.na(
      source_station_number
    ) &
    !is.na(
      primary_site_order
    )
]


comparison[
  ,
  raw_value_available :=
    is.finite(
      tmax_raw
    ) |
    is.finite(
      tmin_raw
    )
]


comparison[
  ,
  no_raw_calendar_row :=
    !raw_source_assigned &
    !raw_value_available
]


comparison[
  ,
  source_key_matched :=
    !is.na(
      declared_source_start_date
    ) &
    !is.na(
      declared_source_end_date
    )
]


comparison[
  ,
  source_period_valid :=
    data.table::fcase(
      
      raw_source_assigned &
        source_key_matched &
        date >=
        declared_source_start_date &
        date <=
        declared_source_end_date,
      TRUE,
      
      no_raw_calendar_row,
      TRUE,
      
      default =
        FALSE
    )
]


# These are genuine failures: a raw value exists but there is no
# associated source-station assignment.

raw_values_without_source <- comparison[
  raw_value_available == TRUE &
    raw_source_assigned == FALSE
]


if (nrow(raw_values_without_source) > 0L) {
  
  fwrite(
    raw_values_without_source,
    file.path(
      table_output_dir,
      "temperature_raw_values_without_source_assignment.csv"
    )
  )
  
  stop(
    nrow(
      raw_values_without_source
    ),
    " rows contain raw temperature values without a source-station assignment.",
    call. = FALSE
  )
}


# These rows carry a source assignment that does not match the declared
# primary-site table or falls outside its declared date interval.

invalid_source_periods <- comparison[
  raw_source_assigned == TRUE &
    source_period_valid == FALSE
]


if (nrow(invalid_source_periods) > 0L) {
  
  fwrite(
    invalid_source_periods,
    file.path(
      table_output_dir,
      "temperature_invalid_primary_site_periods.csv"
    )
  )
  
  stop(
    nrow(
      invalid_source_periods
    ),
    " source-assigned comparison rows fail primary-site period validation.",
    call. = FALSE
  )
}


# Retain a compact audit of rows with no raw source assignment. These
# are expected calendar rows where homogenised data may exist but no
# raw observation is available for that location and date.

no_raw_calendar_audit <- comparison[
  no_raw_calendar_row == TRUE,
  .(
    rows =
      .N,
    
    first_date =
      min(
        date
      ),
    
    last_date =
      max(
        date
      ),
    
    homogenised_tmax_days =
      sum(
        is.finite(
          tmax_hom
        )
      ),
    
    homogenised_tmin_days =
      sum(
        is.finite(
          tmin_hom
        )
      )
  ),
  by = acorn_station_number
]


fwrite(
  no_raw_calendar_audit,
  file.path(
    table_output_dir,
    "temperature_no_raw_calendar_rows_audit.csv"
  )
)

# ============================================================
# Build one canonical station name per ACORN location
# ============================================================

station_names <- comparison[
  !is.na(station_name) & nzchar(station_name),
  .N,
  by = .(acorn_station_number, station_name)
][
  order(acorn_station_number, -N, station_name)
][
  , .SD[1L],
  by = acorn_station_number
][
  , .(acorn_station_number, station_name)
]

if (nrow(station_names) != expected_acorn_locations) {
  stop(
    "Could not assign one canonical station name to every ACORN location.",
    call. = FALSE
  )
}

# Some existing summaries were split into one-day NA-name rows and the
# substantive named station row. Replace the row-level name with one
# canonical non-missing name per ACORN location before any grouping.
comparison[, station_name := NULL]
comparison[
  station_names,
  on = "acorn_station_number",
  station_name := i.station_name
]

# ============================================================
# Derive daily paired analysis fields
# ============================================================

comparison[, `:=`(
  year = as.integer(format(date, "%Y")),
  month = as.integer(format(date, "%m")),
  day_of_year = as.integer(format(date, "%j")),
  tmean_raw = fifelse(
    matched_tmax & matched_tmin,
    (tmax_raw + tmin_raw) / 2,
    NA_real_
  ),
  tmean_hom = fifelse(
    matched_tmax & matched_tmin,
    (tmax_hom + tmin_hom) / 2,
    NA_real_
  )
)]

comparison[, tmean_adjustment := fifelse(
  is.finite(tmean_raw) & is.finite(tmean_hom),
  tmean_hom - tmean_raw,
  NA_real_
)]

paired_daily <- comparison[
  matched_tmax == TRUE | matched_tmin == TRUE,
  .(
    acorn_station_number,
    station_name,
    date,
    year,
    month,
    day_of_year,
    source_station_number,
    primary_site_order,
    primary_period_start=declared_source_start_date,
    primary_period_end=declared_source_end_date,
    tmax_raw,
    tmax_hom,
    tmax_adjustment,
    matched_tmax,
    tmin_raw,
    tmin_hom,
    tmin_adjustment,
    matched_tmin,
    tmean_raw,
    tmean_hom,
    tmean_adjustment,
    matched_tmean = is.finite(tmean_raw) & is.finite(tmean_hom)
  )
]

setorder(paired_daily, acorn_station_number, date)

# ============================================================
# Station-level validation audit
# ============================================================

station_audit <- comparison[, .(
  first_date = min(date),
  last_date = max(date),
  calendar_days = .N,
  source_stations = uniqueN(source_station_number),
  first_source_station = source_station_number[which.min(date)],
  last_source_station = source_station_number[which.max(date)],
  
  raw_tmax_days = sum(is.finite(tmax_raw)),
  hom_tmax_days = sum(is.finite(tmax_hom)),
  matched_tmax_days = sum(matched_tmax),
  tmax_match_fraction_of_raw = sum(matched_tmax) / sum(is.finite(tmax_raw)),
  mean_tmax_adjustment = mean(tmax_adjustment[matched_tmax], na.rm = TRUE),
  median_tmax_adjustment = median(tmax_adjustment[matched_tmax], na.rm = TRUE),
  sd_tmax_adjustment = sd(tmax_adjustment[matched_tmax], na.rm = TRUE),
  min_tmax_adjustment = min(tmax_adjustment[matched_tmax], na.rm = TRUE),
  max_tmax_adjustment = max(tmax_adjustment[matched_tmax], na.rm = TRUE),
  proportion_tmax_changed = mean(
    abs(tmax_adjustment[matched_tmax]) > adjustment_tolerance,
    na.rm = TRUE
  ),
  
  raw_tmin_days = sum(is.finite(tmin_raw)),
  hom_tmin_days = sum(is.finite(tmin_hom)),
  matched_tmin_days = sum(matched_tmin),
  tmin_match_fraction_of_raw = sum(matched_tmin) / sum(is.finite(tmin_raw)),
  mean_tmin_adjustment = mean(tmin_adjustment[matched_tmin], na.rm = TRUE),
  median_tmin_adjustment = median(tmin_adjustment[matched_tmin], na.rm = TRUE),
  sd_tmin_adjustment = sd(tmin_adjustment[matched_tmin], na.rm = TRUE),
  min_tmin_adjustment = min(tmin_adjustment[matched_tmin], na.rm = TRUE),
  max_tmin_adjustment = max(tmin_adjustment[matched_tmin], na.rm = TRUE),
  proportion_tmin_changed = mean(
    abs(tmin_adjustment[matched_tmin]) > adjustment_tolerance,
    na.rm = TRUE
  ),
  
  matched_tmean_days = sum(
    is.finite(tmean_raw) & is.finite(tmean_hom)
  )
), by = acorn_station_number]

station_audit <- station_names[
  station_audit,
  on = "acorn_station_number"
]

setorder(station_audit, acorn_station_number)

if (nrow(station_audit) != expected_acorn_locations) {
  stop("Station audit does not contain exactly 112 locations.", call. = FALSE)
}

# ============================================================
# Annual paired series
# ============================================================

is_leap_year <- function(year) {
  (year %% 4L == 0L & year %% 100L != 0L) | year %% 400L == 0L
}

annual_paired <- paired_daily[, {
  expected_days <- if (is_leap_year(year[1L])) 366L else 365L
  
  n_tmax <- sum(matched_tmax)
  n_tmin <- sum(matched_tmin)
  n_tmean <- sum(matched_tmean)
  
  list(
    expected_days = expected_days,
    matched_tmax_days = n_tmax,
    matched_tmin_days = n_tmin,
    matched_tmean_days = n_tmean,
    tmax_completeness = n_tmax / expected_days,
    tmin_completeness = n_tmin / expected_days,
    tmean_completeness = n_tmean / expected_days,
    tmax_annual_raw = if (
      n_tmax / expected_days >= annual_completeness_threshold
    ) mean(tmax_raw[matched_tmax], na.rm = TRUE) else NA_real_,
    tmax_annual_hom = if (
      n_tmax / expected_days >= annual_completeness_threshold
    ) mean(tmax_hom[matched_tmax], na.rm = TRUE) else NA_real_,
    tmin_annual_raw = if (
      n_tmin / expected_days >= annual_completeness_threshold
    ) mean(tmin_raw[matched_tmin], na.rm = TRUE) else NA_real_,
    tmin_annual_hom = if (
      n_tmin / expected_days >= annual_completeness_threshold
    ) mean(tmin_hom[matched_tmin], na.rm = TRUE) else NA_real_,
    tmean_annual_raw = if (
      n_tmean / expected_days >= annual_completeness_threshold
    ) mean(tmean_raw[matched_tmean], na.rm = TRUE) else NA_real_,
    tmean_annual_hom = if (
      n_tmean / expected_days >= annual_completeness_threshold
    ) mean(tmean_hom[matched_tmean], na.rm = TRUE) else NA_real_
  )
}, by = .(acorn_station_number, station_name, year)]

annual_paired[, `:=`(
  tmax_annual_adjustment = tmax_annual_hom - tmax_annual_raw,
  tmin_annual_adjustment = tmin_annual_hom - tmin_annual_raw,
  tmean_annual_adjustment = tmean_annual_hom - tmean_annual_raw
)]

setorder(annual_paired, acorn_station_number, year)

# ============================================================
# Match-status and source-transition audits
# ============================================================

match_status_audit <- rbindlist(
  list(
    
    comparison[
      ,
      .(
        element =
          "tmax",
        
        match_status =
          tmax_match_status
      )
    ][
      ,
      .(
        rows =
          as.integer(.N)
      ),
      by = .(
        element,
        match_status
      )
    ],
    
    comparison[
      ,
      .(
        element =
          "tmin",
        
        match_status =
          tmin_match_status
      )
    ][
      ,
      .(
        rows =
          as.integer(.N)
      ),
      by = .(
        element,
        match_status
      )
    ]
  ),
  
  use.names = TRUE,
  fill = TRUE
)

setorder(
  match_status_audit,
  element,
  match_status
)

source_transition_audit <- primary_sites[
  order(acorn_station_number, primary_site_order),
  .(
    source_station_count = .N,
    first_date = min(start_date),
    last_date = max(end_date),
    source_station_sequence = paste(source_station_number, collapse = ";"),
    transition_dates = paste(
      format(start_date[-1L], "%Y-%m-%d"),
      collapse = ";"
    )
  ),
  by = acorn_station_number
]

source_transition_audit <- station_names[
  source_transition_audit,
  on = "acorn_station_number"
]

# ============================================================
# Global validation summary
# ============================================================

validation_summary <- data.table(
  metric = c(
    "ACORN locations",
    "Primary-site rows",
    "Comparison rows",
    "Paired daily rows",
    "First comparison date",
    "Last comparison date",
    "Matched Tmax days",
    "Matched Tmin days",
    "Matched Tmean days",
    "Maximum Tmax adjustment arithmetic error",
    "Maximum Tmin adjustment arithmetic error",
    "Annual completeness threshold"
  ),
  value = c(
    as.character(location_count),
    as.character(nrow(primary_sites)),
    as.character(nrow(comparison)),
    as.character(nrow(paired_daily)),
    as.character(min(comparison$date)),
    as.character(max(comparison$date)),
    as.character(sum(comparison$matched_tmax)),
    as.character(sum(comparison$matched_tmin)),
    as.character(sum(is.finite(comparison$tmean_raw))),
    format(tmax_adjustment_error, scientific = TRUE),
    format(tmin_adjustment_error, scientific = TRUE),
    as.character(annual_completeness_threshold)
  )
)

# ============================================================
# Save validated products
# ============================================================

write_parquet(
  paired_daily,
  file.path(output_dir, "temperature_raw_homogenised_paired_daily_v2.3.parquet"),
  compression = "zstd"
)

write_parquet(
  annual_paired,
  file.path(output_dir, "temperature_raw_homogenised_paired_annual_v2.3.parquet"),
  compression = "zstd"
)

write_parquet(
  station_audit,
  file.path(output_dir, "temperature_raw_homogenised_station_audit_v2.3.parquet"),
  compression = "zstd"
)

fwrite(
  station_audit,
  file.path(table_output_dir, "temperature_raw_homogenised_station_audit_v2.3.csv")
)

fwrite(
  validation_summary,
  file.path(table_output_dir, "temperature_input_validation_summary_v2.3.csv")
)

fwrite(
  match_status_audit,
  file.path(table_output_dir, "temperature_match_status_audit_v2.3.csv")
)

fwrite(
  source_transition_audit,
  file.path(table_output_dir, "temperature_source_transition_audit_v2.3.csv")
)

# ============================================================
# Console summary
# ============================================================

message("\n========================================")
message("TEMPERATURE ANALYSIS INPUTS VALIDATED")
message("========================================\n")

print(validation_summary)

message("\nMatch-status audit:")
print(match_status_audit)

message("\nStation audit summary:")
print(station_audit[, .(
  locations = .N,
  earliest_start = min(first_date),
  latest_start = max(first_date),
  minimum_matched_tmax_days = min(matched_tmax_days),
  median_matched_tmax_days = median(matched_tmax_days),
  minimum_matched_tmin_days = min(matched_tmin_days),
  median_matched_tmin_days = median(matched_tmin_days),
  median_tmax_adjustment = median(mean_tmax_adjustment, na.rm = TRUE),
  median_tmin_adjustment = median(mean_tmin_adjustment, na.rm = TRUE)
)])

message("\nAnnual paired-series coverage:")
print(annual_paired[, .(
  station_years = .N,
  valid_tmax_station_years = sum(is.finite(tmax_annual_raw)),
  valid_tmin_station_years = sum(is.finite(tmin_annual_raw)),
  valid_tmean_station_years = sum(is.finite(tmean_annual_raw))
)])
