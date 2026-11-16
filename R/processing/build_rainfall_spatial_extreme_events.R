# ============================================================
# build_rainfall_spatial_extreme_events.R
#
# Build daily station-level extreme-rainfall observations and
# daily national spatial-concurrence summaries for the frozen
# 742-station Australian rainfall network.
#
# Version 1 scope:
#   - P95, P97.5 and P99 station thresholds
#   - station-day eligibility
#   - unweighted station concurrence
#   - eligible-area-adjusted Voronoi concurrence
#   - daily magnitude/extremity summaries
#   - annual exposure/coverage audit
#
# Deliberately excluded from this version:
#   - distance-connected daily clusters
#   - multi-day episode matching
#   - trend models
#
# Required existing products:
#   data-derived/rainfall/network/
#     rainfall_continental_core_network_1950_2025.parquet
#
#   data-derived/rainfall/pot/
#     rainfall_pot_thresholds_1961_1990.parquet
#
#   data-derived/rainfall/spatial_extremes/
#     rainfall_station_voronoi_weights.parquet
#
#   data-derived/rainfall/daily/
#     one or more ADAM daily-rainfall parquet files
#
# IMPORTANT:
#   Adjust daily_input_dir and daily_file_pattern below if the
#   existing five-year ADAM parquet directory uses another path.
# ============================================================

# ============================================================
# Configuration
# ============================================================

analysis_start_date <- as.Date("1950-01-01")
analysis_end_date <- as.Date("2025-12-31")

network_file <- file.path(
  "data-derived",
  "rainfall",
  "network",
  "rainfall_continental_core_network_1950_2025.parquet"
)

threshold_file <- file.path("data-derived",
                            "rainfall",
                            "pot",
                            "rainfall_pot_thresholds_1961_1990.parquet")

voronoi_weight_file <- file.path(
  "data-derived",
  "rainfall",
  "spatial_extremes",
  "rainfall_station_voronoi_weights.parquet"
)

# EDIT THIS PATH ONLY IF YOUR EXISTING ADAM five-year daily parquet
# files live elsewhere. The expected filenames include strings such
# as adam_prcp_days_1950_1954.parquet.
daily_input_dir <- file.path("data-derived", "rainfall", "adam_prcp_days")

daily_file_pattern <- "^adam_prcp_days_[0-9]{4}_[0-9]{4}\\.parquet$"

output_dir <- file.path("data-derived", "rainfall", "spatial_extremes")

table_output_dir <- file.path("tables", "rainfall_spatial_extremes")

station_day_output_file <- file.path(output_dir,
                                     "rainfall_spatial_extreme_station_days_1950_2025.parquet")

daily_output_file <- file.path(output_dir,
                               "rainfall_spatial_extreme_days_1950_2025.parquet")

annual_output_file <- file.path(output_dir,
                                "rainfall_spatial_extreme_annual_coverage_1950_2025.parquet")

threshold_ids <- c("p95", "p975", "p99")
expected_station_count <- 742L

# A station-day is eligible only if it is a valid one-day rainfall
# observation and has a valid threshold. These values match the
# primary GPD quality rule.
accepted_rainfall_quality <- c(0L, 1L, 2L, 3L, 4L)

# Minimum daily coverage flags are diagnostic only. Days are not
# discarded from the saved product; later analyses can explicitly
# choose the required exposure standard.
minimum_station_coverage_fraction <- 0.90
minimum_area_coverage_fraction <- 0.90

# Save the full station-day product. This is useful for validation and
# later connected-cluster construction, but is larger than the daily
# summary. Set FALSE only if storage is constrained.
save_station_day_product <- TRUE

# ============================================================
# Packages
# ============================================================

required_packages <- c("data.table", "arrow")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

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
dir.create(table_output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

first_existing_name <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0L)
    NA_character_
  else
    hit[1L]
}

standardise_daily_columns <- function(x, source_file) {
  required_columns <- c(
    "station_number",
    "date",
    "rainfall_mm",
    "observation_quality_flag",
    "rainfall_period",
    "rainfall_type",
    "rainfall_quality",
    "rainfall_day_count"
  )
  
  missing_columns <- setdiff(required_columns, names(x))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns in ",
      basename(source_file),
      ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  
  x[, .(
    station_number =
      as.character(station_number),
    
    date =
      as.Date(date),
    
    rainfall_mm =
      as.numeric(rainfall_mm),
    
    observation_quality_flag =
      as.integer(observation_quality_flag),
    
    rainfall_period =
      as.numeric(rainfall_period),
    
    rainfall_type =
      as.character(rainfall_type),
    
    rainfall_quality =
      as.integer(rainfall_quality),
    
    rainfall_day_count =
      as.numeric(rainfall_day_count)
  )]
}

safe_numeric_summary <- function(x, fun, ...) {
  x <- x[is.finite(x)]
  if (length(x) == 0L)
    return(NA_real_)
  as.numeric(fun(x, ...))
}

# ============================================================
# Read and validate static products
# ============================================================

required_static_files <- c(network_file, threshold_file, voronoi_weight_file)

missing_static_files <- required_static_files[!file.exists(required_static_files)]

if (length(missing_static_files) > 0L) {
  stop(
    "Missing required static products:\n",
    paste(missing_static_files, collapse = "\n"),
    call. = FALSE
  )
}

network <- as.data.table(read_parquet(network_file))
thresholds <- as.data.table(read_parquet(threshold_file))
weights <- as.data.table(read_parquet(voronoi_weight_file))

if (nrow(network) != expected_station_count) {
  stop(
    "Expected ",
    expected_station_count,
    " stations in the frozen network, found ",
    nrow(network),
    ".",
    call. = FALSE
  )
}

if (uniqueN(network$station_number) != expected_station_count) {
  stop("Frozen network station numbers are not unique.", call. = FALSE)
}

if (nrow(weights) != expected_station_count) {
  stop(
    "Expected ",
    expected_station_count,
    " Voronoi weight rows, found ",
    nrow(weights),
    ".",
    call. = FALSE
  )
}

required_weight_columns <- c("station_number",
                             "voronoi_area_km2",
                             "continental_area_weight")

missing_weight_columns <- setdiff(required_weight_columns, names(weights))

if (length(missing_weight_columns) > 0L) {
  stop(
    "Voronoi weight product is missing: ",
    paste(missing_weight_columns, collapse = ", "),
    call. = FALSE
  )
}

if (abs(sum(weights$continental_area_weight) - 1) > 1e-8) {
  stop("Voronoi continental-area weights do not sum to one.", call. = FALSE)
}

required_threshold_columns <- c("station_number", "threshold_id", "threshold_mm")

missing_threshold_columns <- setdiff(required_threshold_columns, names(thresholds))

if (length(missing_threshold_columns) > 0L) {
  stop(
    "Threshold product is missing: ",
    paste(missing_threshold_columns, collapse = ", "),
    call. = FALSE
  )
}

thresholds <- thresholds[threshold_id %in% threshold_ids]

# Respect the validity flag if the threshold product includes it.
if ("reference_threshold_valid" %in% names(thresholds)) {
  thresholds <- thresholds[reference_threshold_valid == TRUE]
}

thresholds <- thresholds[is.finite(threshold_mm) & threshold_mm > 0]

if (nrow(thresholds) !=
    expected_station_count * length(threshold_ids)) {
  threshold_audit <- thresholds[, .(stations = uniqueN(station_number), rows = .N), by = threshold_id]
  
  print(threshold_audit)
  
  stop("Expected one valid P95/P97.5/P99 threshold for each frozen-network station.",
       call. = FALSE)
}

# Static station metadata used on the station-day product.
station_static <- merge(network[, .(station_number = as.character(station_number),
                                    station_name,
                                    state,
                                    latitude,
                                    longitude)],
                        weights[, .(station_number = as.character(station_number),
                                    voronoi_area_km2,
                                    continental_area_weight)],
                        by = "station_number",
                        all = TRUE,
                        sort = FALSE)

if (nrow(station_static) != expected_station_count ||
    anyNA(station_static$voronoi_area_km2)) {
  stop("Network-to-Voronoi join failed.", call. = FALSE)
}

thresholds[, station_number := as.character(station_number)]
thresholds <- thresholds[, .(
  station_number,
  threshold_id,
  threshold_probability = if ("threshold_probability" %in% names(thresholds))
    as.numeric(threshold_probability)
  else
    fcase(
      threshold_id == "p95",
      0.95,
      threshold_id == "p975",
      0.975,
      threshold_id == "p99",
      0.99,
      default = NA_real_
    ),
  threshold_mm = as.numeric(threshold_mm)
)]

setkey(station_static, station_number)
setkey(thresholds, station_number)

# ============================================================
# Find daily rainfall parquet files
# ============================================================

if (!dir.exists(daily_input_dir)) {
  stop(
    "Daily-rainfall directory not found: ",
    daily_input_dir,
    "\nUpdate daily_input_dir at the top of the script.",
    call. = FALSE
  )
}

daily_files <- list.files(daily_input_dir, pattern = daily_file_pattern, full.names = TRUE)

daily_files <- sort(daily_files)

if (length(daily_files) == 0L) {
  stop(
    "No ADAM daily-rainfall parquet files matched ",
    daily_file_pattern,
    " under ",
    daily_input_dir,
    ".",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Identify file coverage years
# ------------------------------------------------------------

file_years <- data.table(
  file = daily_files,
  
  file_start_year =
    as.integer(sub(
      ".*_([0-9]{4})_[0-9]{4}\\.parquet$", "\\1", daily_files
    )),
  
  file_end_year =
    as.integer(sub(
      ".*_[0-9]{4}_([0-9]{4})\\.parquet$", "\\1", daily_files
    ))
)


# ------------------------------------------------------------
# Retain only files overlapping the analysis period
# ------------------------------------------------------------

file_years <- file_years[file_end_year >=
                           as.integer(format(analysis_start_date, "%Y")) &
                           file_start_year <=
                           as.integer(format(analysis_end_date, "%Y"))]


daily_files <- file_years$file


if (length(daily_files) == 0L) {
  stop("No daily rainfall parquet files overlap the analysis period.",
       call. = FALSE)
}


message("Daily rainfall files overlapping analysis period: ",
        length(daily_files))


print(file_years)

# ============================================================
# Process daily rainfall files
# ============================================================

daily_summary_list <- vector("list", length(daily_files))
file_audit_list <- vector("list", length(daily_files))

station_day_chunk_dir <- file.path(
  output_dir,
  "station_day_chunks"
)

dir.create(
  station_day_chunk_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (i in seq_along(daily_files)) {
  file_path <- daily_files[i]
  
  message("Spatial concurrence file ",
          i,
          " of ",
          length(daily_files),
          ": ",
          basename(file_path))
  
  raw <- as.data.table(read_parquet(file_path))
  daily <- standardise_daily_columns(raw, file_path)
  
  raw_rows <- nrow(daily)
  
  daily <- daily[station_number %in%
                   station_static$station_number &
                   date >=
                   analysis_start_date &
                   date <=
                   analysis_end_date]
  
  selected_rows <- nrow(daily)
  
  # Prevent many-to-many expansion downstream. The prior archive audit
  # should already have resolved station-date duplicates; stop if any
  # remain in the selected series.
  duplicate_rows <- daily[, .N, by = .(station_number, date)][N > 1L]
  
  if (nrow(duplicate_rows) > 0L) {
    stop(
      "Duplicate station-date rows remain in ",
      basename(file_path),
      ". Inspect the archive duplicate audit before concurrence processing.",
      call. = FALSE
    )
  }
  
  daily <- station_static[daily, on = "station_number"]
  
  # Primary station-day eligibility. The rainfall-period fields are used
  # to exclude accumulated/non-one-day totals. If one of these fields was
  # absent in the source, standardise_daily_columns supplied the one-day
  # default; audit that assumption through source-column inspection before
  # accepting results for a differently structured archive.
  # ------------------------------------------------------------
  # Daily observation eligibility
  #
  # Valid zero-rainfall observations in ADAM normally have:
  #
  #   rainfall_mm        = 0
  #   rainfall_period    = NA
  #   rainfall_day_count = NA
  #   rainfall_type      = NA
  #
  # These are valid dry days and must remain in the exposure
  # denominator.
  #
  # Positive rainfall observations must represent a one-day total.
  # ------------------------------------------------------------
  
  daily[, valid_quality :=
          rainfall_quality %in%
          accepted_rainfall_quality]
  
  
  daily[, valid_dry_day :=
          is.finite(rainfall_mm) &
          rainfall_mm == 0 &
          valid_quality]
  
  
  daily[, valid_one_day_rain :=
          is.finite(rainfall_mm) &
          rainfall_mm > 0 &
          valid_quality &
          !is.na(rainfall_period) &
          rainfall_period == 1 &
          !is.na(rainfall_day_count) &
          rainfall_day_count == 1]
  
  
  daily[, eligible_observation :=
          valid_dry_day |
          valid_one_day_rain]
  
  # Cross with the three thresholds. This creates at most three rows per
  # selected station-day, not a complete date × 742 station grid.
  station_days <- thresholds[daily, on = "station_number", allow.cartesian = TRUE]
  
  station_days[, eligible :=
                 eligible_observation &
                 is.finite(threshold_mm) &
                 threshold_mm > 0]
  
  station_days[, exceeds :=
                 eligible &
                 rainfall_mm > threshold_mm]
  
  station_days[, threshold_ratio := fifelse(eligible, rainfall_mm / threshold_mm, NA_real_)]
  
  station_days[, excess_mm := fifelse(exceeds, rainfall_mm - threshold_mm, NA_real_)]
  
  # ------------------------------------------------------------
  # Calendar fields
  # ------------------------------------------------------------
  
  station_days[, `:=`(
    year =
      as.integer(format(date, "%Y")),
    
    month =
      as.integer(format(date, "%m")),
    
    day_of_year =
      as.integer(format(date, "%j"))
  )]
  
  
  # ------------------------------------------------------------
  # Season
  # ------------------------------------------------------------
  
  station_days[, season :=
                 data.table::fcase(
                   month %in%
                     c(12L, 1L, 2L),
                   "DJF",
                   
                   month %in%
                     3L:5L,
                   "MAM",
                   
                   month %in%
                     6L:8L,
                   "JJA",
                   
                   month %in%
                     9L:11L,
                   "SON",
                   
                   default =
                     NA_character_
                 )]
  
  
  # ------------------------------------------------------------
  # Represented-area fields
  # ------------------------------------------------------------
  
  station_days[, `:=`(
    eligible_area_km2 =
      data.table::fifelse(eligible, voronoi_area_km2, 0),
    
    exceeding_area_km2 =
      data.table::fifelse(exceeds, voronoi_area_km2, 0)
  )]
  
  
  station_day_chunk <- station_days[, .(
    station_number,
    station_name,
    state,
    latitude,
    longitude,
    
    date,
    year,
    month,
    day_of_year,
    season,
    
    rainfall_mm,
    observation_quality_flag,
    rainfall_period,
    rainfall_type,
    rainfall_quality,
    rainfall_day_count,
    
    valid_dry_day,
    valid_one_day_rain,
    
    threshold_id,
    threshold_probability,
    threshold_mm,
    
    eligible,
    exceeds,
    
    threshold_ratio,
    excess_mm,
    
    voronoi_area_km2,
    continental_area_weight,
    eligible_area_km2,
    exceeding_area_km2
  )]
  
  chunk_output_file <- file.path(
    station_day_chunk_dir,
    sub(
      "\\.parquet$",
      "_spatial_extremes.parquet",
      basename(file_path)
    )
  )
  
  arrow::write_parquet(station_day_chunk, chunk_output_file, compression = "zstd")
  
  daily_summary_list[[i]] <- station_day_chunk[, {
    eligible_index <- eligible == TRUE
    exceed_index <- exceeds == TRUE
    
    eligible_count <- sum(eligible_index)
    exceed_count <- sum(exceed_index)
    
    eligible_area <- sum(eligible_area_km2, na.rm = TRUE)
    
    exceed_area <- sum(exceeding_area_km2, na.rm = TRUE)
    
    list(
      eligible_stations =
        as.integer(eligible_count),
      
      exceeding_stations =
        as.integer(exceed_count),
      
      station_concurrence_fraction =
        if (eligible_count > 0L) {
          as.numeric(exceed_count /
                       eligible_count)
        } else {
          NA_real_
        },
      
      eligible_station_fraction =
        as.numeric(eligible_count /
                     expected_station_count),
      
      eligible_area_km2 =
        as.numeric(eligible_area),
      
      exceeding_area_km2 =
        as.numeric(exceed_area),
      
      area_concurrence_fraction =
        if (eligible_area > 0) {
          as.numeric(exceed_area /
                       eligible_area)
        } else {
          NA_real_
        },
      
      eligible_area_fraction =
        as.numeric(eligible_area /
                     sum(weights$voronoi_area_km2)),
      
      states_eligible =
        as.integer(uniqueN(state[eligible_index])),
      
      states_with_exceedance =
        as.integer(uniqueN(state[exceed_index])),
      
      maximum_rainfall_mm =
        safe_numeric_summary(rainfall_mm[eligible_index], max),
      
      median_rainfall_mm =
        safe_numeric_summary(rainfall_mm[eligible_index], median),
      
      maximum_exceeding_rainfall_mm =
        safe_numeric_summary(rainfall_mm[exceed_index], max),
      
      median_exceeding_rainfall_mm =
        safe_numeric_summary(rainfall_mm[exceed_index], median),
      
      maximum_threshold_ratio =
        safe_numeric_summary(threshold_ratio[eligible_index], max),
      
      median_threshold_ratio =
        safe_numeric_summary(threshold_ratio[eligible_index], median),
      
      maximum_exceeding_threshold_ratio =
        safe_numeric_summary(threshold_ratio[exceed_index], max),
      
      median_exceeding_threshold_ratio =
        safe_numeric_summary(threshold_ratio[exceed_index], median),
      
      maximum_excess_mm =
        safe_numeric_summary(excess_mm[exceed_index], max)
    )
  }, by = .(date,
            year,
            month,
            day_of_year,
            season,
            threshold_id,
            threshold_probability)]
  
  observation_audit <- daily[, .(
    selected_rows =
      .N,
    
    valid_dry_days =
      sum(valid_dry_day),
    
    valid_one_day_rain =
      sum(valid_one_day_rain),
    
    eligible_observations =
      sum(eligible_observation),
    
    positive_multiday_or_invalid =
      sum(is.finite(rainfall_mm) &
            rainfall_mm > 0 &
            !valid_one_day_rain),
    
    missing_rainfall =
      sum(!is.finite(rainfall_mm))
  )]
  
  file_audit_list[[i]] <- data.table(
    file =
      basename(file_path),
    
    raw_rows =
      as.integer(raw_rows),
    
    selected_network_period_rows =
      as.integer(selected_rows),
    
    valid_dry_days =
      as.integer(observation_audit$valid_dry_days),
    
    valid_one_day_rain =
      as.integer(observation_audit$valid_one_day_rain),
    
    eligible_observations =
      as.integer(observation_audit$eligible_observations),
    
    positive_multiday_or_invalid =
      as.integer(observation_audit$positive_multiday_or_invalid),
    
    missing_rainfall =
      as.integer(observation_audit$missing_rainfall),
    
    station_threshold_rows =
      as.integer(nrow(station_days)),
    
    eligible_station_threshold_rows =
      as.integer(sum(station_days$eligible)),
    
    exceedance_rows =
      as.integer(sum(station_days$exceeds))
  )
  
  rm(raw, daily, station_days, station_day_chunk)
  
  invisible(gc(full = TRUE))
}

# ============================================================
# Combine small daily summaries and file audits
# ============================================================

daily_spatial <- rbindlist(daily_summary_list, use.names = TRUE, fill = TRUE)

file_audit <- rbindlist(file_audit_list, use.names = TRUE, fill = TRUE)

rm(daily_summary_list, file_audit_list)

invisible(gc(full = TRUE))

setorder(daily_spatial, date, threshold_probability)

######################################################
# Diagnostic flags. Preserve all rows and let later analyses select the
# desired coverage standard explicitly.
daily_spatial[, station_coverage_adequate :=
                eligible_station_fraction >= minimum_station_coverage_fraction]

daily_spatial[, area_coverage_adequate :=
                eligible_area_fraction >= minimum_area_coverage_fraction]

daily_spatial[, primary_coverage_adequate :=
                station_coverage_adequate & area_coverage_adequate]

setorder(daily_spatial, date, threshold_probability)

# ============================================================
# Coverage tiers for spatial-extreme analyses
#
# Primary inference uses >=80% of both stations and represented
# area. >=90% is retained as a strict sensitivity analysis.
# ============================================================

daily_spatial[, coverage_strict :=
                eligible_station_fraction >= 0.90 &
                eligible_area_fraction >= 0.90]

daily_spatial[, coverage_primary :=
                eligible_station_fraction >= 0.80 &
                eligible_area_fraction >= 0.80]

daily_spatial[, coverage_tier :=
                data.table::fcase(coverage_strict,
                                  "high_90",
                                  
                                  coverage_primary,
                                  "acceptable_80",
                                  
                                  default = "low")]

# ============================================================
# Validate complete date × threshold coverage
# ============================================================

calendar <- data.table(date = seq(analysis_start_date, analysis_end_date, by = "day"))

calendar_grid <- CJ(date = calendar$date,
                    threshold_id = threshold_ids,
                    unique = TRUE)

coverage_check <- daily_spatial[calendar_grid, on = .(date, threshold_id)]

missing_daily_rows <- coverage_check[is.na(eligible_stations)]

if (nrow(missing_daily_rows) > 0L) {
  warning(
    nrow(missing_daily_rows),
    " date-threshold combinations are absent from the daily archive summary. ",
    "They remain identifiable in the annual audit but are not imputed.",
    call. = FALSE
  )
}

# ============================================================
# Annual coverage and concurrence audit
# ============================================================

annual_coverage <- daily_spatial[, .(
  days = as.integer(.N),
  days_primary_coverage_adequate = as.integer(sum(primary_coverage_adequate)),
  fraction_days_primary_coverage_adequate = as.numeric(mean(primary_coverage_adequate)),
  median_eligible_stations = as.numeric(median(eligible_stations, na.rm = TRUE)),
  minimum_eligible_stations = as.integer(min(eligible_stations, na.rm = TRUE)),
  q10_eligible_stations = as.numeric(quantile(
    eligible_stations, 0.10, na.rm = TRUE, names = FALSE
  )),
  median_eligible_station_fraction = as.numeric(median(eligible_station_fraction, na.rm = TRUE)),
  median_eligible_area_km2 = as.numeric(median(eligible_area_km2, na.rm = TRUE)),
  minimum_eligible_area_km2 = as.numeric(min(eligible_area_km2, na.rm = TRUE)),
  q10_eligible_area_km2 = as.numeric(quantile(
    eligible_area_km2, 0.10, na.rm = TRUE, names = FALSE
  )),
  median_eligible_area_fraction = as.numeric(median(eligible_area_fraction, na.rm = TRUE)),
  days_with_exceedance = as.integer(sum(exceeding_stations > 0L)),
  maximum_exceeding_stations = as.integer(max(exceeding_stations, na.rm = TRUE)),
  maximum_station_concurrence_fraction = as.numeric(max(station_concurrence_fraction, na.rm = TRUE)),
  maximum_exceeding_area_km2 = as.numeric(max(exceeding_area_km2, na.rm = TRUE)),
  maximum_area_concurrence_fraction = as.numeric(max(area_concurrence_fraction, na.rm = TRUE))
), by = .(year, threshold_id, threshold_probability)]

setorder(annual_coverage, year, threshold_probability)

# ============================================================
# Summary diagnostics
# ============================================================

threshold_summary <- daily_spatial[primary_coverage_adequate == TRUE, .(
  days = as.integer(.N),
  days_with_exceedance = as.integer(sum(exceeding_stations > 0L)),
  median_exceeding_stations_on_event_days = safe_numeric_summary(exceeding_stations[exceeding_stations > 0L], median),
  p90_exceeding_stations_on_event_days = safe_numeric_summary(
    exceeding_stations[exceeding_stations > 0L],
    quantile,
    probs = 0.90,
    names = FALSE
  ),
  p99_exceeding_stations_on_event_days = safe_numeric_summary(
    exceeding_stations[exceeding_stations > 0L],
    quantile,
    probs = 0.99,
    names = FALSE
  ),
  median_area_fraction_on_event_days = safe_numeric_summary(area_concurrence_fraction[exceeding_stations > 0L], median),
  p90_area_fraction_on_event_days = safe_numeric_summary(
    area_concurrence_fraction[exceeding_stations > 0L],
    quantile,
    probs = 0.90,
    names = FALSE
  ),
  p99_area_fraction_on_event_days = safe_numeric_summary(
    area_concurrence_fraction[exceeding_stations > 0L],
    quantile,
    probs = 0.99,
    names = FALSE
  ),
  correlation_station_area_fraction = as.numeric(
    cor(
      station_concurrence_fraction,
      area_concurrence_fraction,
      use = "complete.obs"
    )
  )
), by = .(threshold_id, threshold_probability)]

largest_days <- daily_spatial[primary_coverage_adequate == TRUE &
                                exceeding_stations > 0L][order(threshold_probability, -area_concurrence_fraction), head(.SD, 20L), by = threshold_id]

# ============================================================
# Save products
# ============================================================

# Detailed station-day data have already been written incrementally
# to:
#
#   data-derived/rainfall/spatial_extremes/station_day_chunks/
#
# Do not concatenate those files here. The directory can later be
# opened directly as an Arrow dataset.


write_parquet(daily_spatial, daily_output_file, compression = "zstd")


write_parquet(annual_coverage, annual_output_file, compression = "zstd")


fwrite(
  file_audit,
  file.path(table_output_dir, "rainfall_spatial_extreme_file_audit.csv")
)


fwrite(
  threshold_summary,
  file.path(
    table_output_dir,
    "rainfall_spatial_extreme_threshold_summary.csv"
  )
)


fwrite(
  annual_coverage,
  file.path(
    table_output_dir,
    "rainfall_spatial_extreme_annual_coverage.csv"
  )
)


fwrite(
  largest_days,
  file.path(
    table_output_dir,
    "rainfall_spatial_extreme_largest_days.csv"
  )
)

# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nDAILY RAINFALL SPATIAL EXTREMES COMPLETE",
  "\n========================================"
)

print(file_audit)
print(threshold_summary)

message("\nAnnual coverage summary:")

print(annual_coverage[, .(
  years =
    as.integer(.N),
  
  minimum_fraction_days_adequate =
    as.numeric(min(
      fraction_days_primary_coverage_adequate, na.rm = TRUE
    )),
  
  median_fraction_days_adequate =
    as.numeric(median(
      fraction_days_primary_coverage_adequate, na.rm = TRUE
    )),
  
  minimum_median_eligible_station_fraction =
    as.numeric(min(
      median_eligible_station_fraction, na.rm = TRUE
    )),
  
  minimum_median_eligible_area_fraction =
    as.numeric(min(median_eligible_area_fraction, na.rm = TRUE))
), by = threshold_id])

message("\nTop area-weighted concurrence days by threshold:")
print(largest_days[, .(
  threshold_id,
  date,
  eligible_stations,
  exceeding_stations,
  station_concurrence_fraction,
  eligible_area_km2,
  exceeding_area_km2,
  area_concurrence_fraction,
  states_with_exceedance,
  maximum_exceeding_rainfall_mm,
  maximum_exceeding_threshold_ratio
)])
