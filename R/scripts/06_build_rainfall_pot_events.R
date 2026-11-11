# ============================================================
# 06_build_rainfall_pot_events.R
# ============================================================

library(data.table)
library(arrow)

source(
  "R/processing/build_rainfall_pot_events.R"
)


# ============================================================
# Paths
# ============================================================

rainfall_parquet_dir <- file.path(
  "data-derived",
  "rainfall",
  "adam_prcp_days"
)

core_network_file <- file.path(
  "data-derived",
  "rainfall",
  "network",
  "rainfall_continental_core_network_1950_2025.parquet"
)

pot_output_dir <- file.path(
  "data-derived",
  "rainfall",
  "pot"
)

pot_table_dir <- file.path(
  "tables",
  "rainfall_pot"
)


dir.create(
  pot_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  pot_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Inputs
# ============================================================

parquet_files <- sort(
  list.files(
    rainfall_parquet_dir,
    pattern = "[.]parquet$",
    full.names = TRUE
  )
)

if (length(parquet_files) != 26L) {
  
  warning(
    "Expected 26 rainfall Parquet files but found ",
    length(parquet_files),
    ".",
    call. = FALSE
  )
}


if (!file.exists(core_network_file)) {
  
  stop(
    "Continental core network does not exist: ",
    core_network_file
  )
}


core_network <- data.table::as.data.table(
  arrow::read_parquet(
    core_network_file
  )
)


# ============================================================
# Build POT events
# ============================================================

rainfall_pot <- build_rainfall_pot_events(
  
  parquet_files =
    parquet_files,
  
  core_network =
    core_network,
  
  analysis_start_year =
    1950L,
  
  analysis_end_year =
    2025L,
  
  reference_start_year =
    1961L,
  
  reference_end_year =
    1990L,
  
  threshold_probabilities =
    c(
      0.95,
      0.975,
      0.99
    ),
  
  run_lengths =
    c(
      2L,
      3L,
      5L
    ),
  
  primary_run_length =
    3L,
  
  wet_day_threshold_mm =
    1,
  
  minimum_reference_wet_days =
    300L,
  
  quantile_type =
    8L,
  
  calculate_event_span_totals =
    TRUE
)


# ============================================================
# Validate
# ============================================================

print(
  rainfall_pot$threshold_summary
)

print(
  rainfall_pot$declustering_sensitivity
)


duplicate_events <- rainfall_pot$events_all[
  ,
  .N,
  by = .(
    station_number,
    threshold_id,
    run_length_days,
    event_id
  )
][
  N > 1L
]


if (nrow(duplicate_events) > 0L) {
  
  stop(
    "Duplicate POT event identifiers found.",
    call. = FALSE
  )
}


invalid_excess <- rainfall_pot$events_all[
  is.na(excess_mm) |
    excess_mm <= 0
]


if (nrow(invalid_excess) > 0L) {
  
  stop(
    "POT events contain missing or non-positive excesses.",
    call. = FALSE
  )
}


# ============================================================
# Save thresholds and raw exceedances
# ============================================================

arrow::write_parquet(
  
  rainfall_pot$thresholds,
  
  file.path(
    pot_output_dir,
    "rainfall_pot_thresholds_1961_1990.parquet"
  ),
  
  compression =
    "zstd"
)


arrow::write_parquet(
  
  rainfall_pot$exceedances,
  
  file.path(
    pot_output_dir,
    "rainfall_pot_daily_exceedances_1950_2025.parquet"
  ),
  
  compression =
    "zstd"
)


# ============================================================
# Save all declustering combinations
# ============================================================

arrow::write_parquet(
  
  rainfall_pot$events_all,
  
  file.path(
    pot_output_dir,
    "rainfall_pot_events_all_thresholds_runs_1950_2025.parquet"
  ),
  
  compression =
    "zstd"
)


# ============================================================
# Save primary 3-day event datasets separately
# ============================================================

for (
  threshold_name in
  c(
    "p95",
    "p975",
    "p99"
  )
) {
  
  threshold_events <- rainfall_pot$primary_events[
    threshold_id ==
      threshold_name
  ]
  
  arrow::write_parquet(
    
    threshold_events,
    
    file.path(
      pot_output_dir,
      paste0(
        "rainfall_pot_events_",
        threshold_name,
        "_run3_1950_2025.parquet"
      )
    ),
    
    compression =
      "zstd"
  )
}


# ============================================================
# Save summaries
# ============================================================

data.table::fwrite(
  
  rainfall_pot$station_summary,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_station_summary.csv"
  )
)


data.table::fwrite(
  
  rainfall_pot$annual_summary,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_annual_summary.csv"
  )
)


data.table::fwrite(
  
  rainfall_pot$declustering_sensitivity,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_declustering_sensitivity.csv"
  )
)


data.table::fwrite(
  
  rainfall_pot$threshold_summary,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_threshold_summary.csv"
  )
)


data.table::fwrite(
  
  rainfall_pot$daily_audit,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_daily_extraction_audit.csv"
  )
)


# ============================================================
# Final summary
# ============================================================

pot_build_summary <- data.table::data.table(
  
  stations =
    data.table::uniqueN(
      rainfall_pot$events_all$
        station_number
    ),
  
  thresholds =
    data.table::uniqueN(
      rainfall_pot$events_all$
        threshold_id
    ),
  
  run_lengths =
    data.table::uniqueN(
      rainfall_pot$events_all$
        run_length_days
    ),
  
  daily_exceedances =
    nrow(
      rainfall_pot$exceedances
    ),
  
  all_declustered_event_records =
    nrow(
      rainfall_pot$events_all
    ),
  
  primary_run3_events =
    nrow(
      rainfall_pot$primary_events
    ),
  
  first_event_date =
    min(
      rainfall_pot$events_all$
        event_peak_date
    ),
  
  last_event_date =
    max(
      rainfall_pot$events_all$
        event_peak_date
    )
)


data.table::fwrite(
  
  pot_build_summary,
  
  file.path(
    pot_table_dir,
    "rainfall_pot_build_summary.csv"
  )
)


message(
  "\n========================================",
  "\nRAINFALL POT EVENT BUILD COMPLETE",
  "\n========================================"
)


print(
  pot_build_summary
)