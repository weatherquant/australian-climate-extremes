# ============================================================
# 04_build_station_metadata.R
# ============================================================

library(data.table)
library(arrow)
library(yaml)

source(
  "R/ingest/read_adam_stations.R"
)

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

cfg <- yaml::read_yaml(
  "config.local.yml"
)

# ------------------------------------------------------------
# Source location
# ------------------------------------------------------------

station_metadata_file <- file.path(
  cfg$paths$raw_root,
  "stations",
  "metadata",
  "ADAM_STATIONS",
  "source",
  "adam_stations_metadata.csv"
)

if (!file.exists(station_metadata_file)) {
  
  stop(
    "ADAM station metadata file not found: ",
    station_metadata_file
  )
}

# ------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------

adam_stations <- read_adam_stations(
  station_metadata_file
)


# ------------------------------------------------------------
# Validate metadata
# ------------------------------------------------------------

station_validation <- validate_adam_stations(
  adam_stations
)


print(
  station_validation
)


# ------------------------------------------------------------
# Load rainfall station-year summary
# ------------------------------------------------------------

rainfall_station_year_file <- file.path(
  
  "data-derived",
  "rainfall",
  "summaries",
  "adam_rainfall_station_year_1900_2026.parquet"
)


if (!file.exists(rainfall_station_year_file)) {
  
  stop(
    "Rainfall station-year summary does not exist: ",
    rainfall_station_year_file
  )
}


rainfall_station_year <-
  data.table::as.data.table(
    
    arrow::read_parquet(
      rainfall_station_year_file
    )
  )


# ------------------------------------------------------------
# Rainfall-metadata coverage audit
# ------------------------------------------------------------

metadata_coverage <-
  audit_rainfall_station_metadata(
    
    rainfall_station_year =
      rainfall_station_year,
    
    stations =
      adam_stations
  )


print(
  metadata_coverage$summary
)


# ------------------------------------------------------------
# Join metadata
# ------------------------------------------------------------

rainfall_station_year_metadata <-
  join_rainfall_station_metadata(
    
    rainfall_station_year =
      rainfall_station_year,
    
    stations =
      adam_stations
  )


# ============================================================
# Output
# ============================================================

metadata_output_dir <- file.path(
  "data-derived",
  "station_metadata"
)


rainfall_summary_dir <- file.path(
  "data-derived",
  "rainfall",
  "summaries"
)


audit_dir <- file.path(
  "data-derived",
  "rainfall",
  "audit"
)


dir.create(
  metadata_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Canonical ADAM station metadata
# ------------------------------------------------------------

arrow::write_parquet(
  
  adam_stations,
  
  file.path(
    metadata_output_dir,
    "adam_stations.parquet"
  ),
  
  compression = "zstd"
)


# ------------------------------------------------------------
# Rainfall station-year plus metadata
# ------------------------------------------------------------

arrow::write_parquet(
  
  rainfall_station_year_metadata,
  
  file.path(
    rainfall_summary_dir,
    "adam_rainfall_station_year_metadata_1900_2026.parquet"
  ),
  
  compression = "zstd"
)


# ------------------------------------------------------------
# Audits
# ------------------------------------------------------------

data.table::fwrite(
  
  station_validation,
  
  file.path(
    audit_dir,
    "adam_stations_validation.csv"
  )
)


data.table::fwrite(
  
  metadata_coverage$summary,
  
  file.path(
    audit_dir,
    "adam_rainfall_metadata_coverage.csv"
  )
)


data.table::fwrite(
  
  metadata_coverage$unmatched_rainfall,
  
  file.path(
    audit_dir,
    "adam_rainfall_stations_without_metadata.csv"
  )
)


data.table::fwrite(
  
  metadata_coverage$metadata_without_rainfall,
  
  file.path(
    audit_dir,
    "adam_metadata_stations_without_rainfall.csv"
  )
)


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

message(
  "\n========================================",
  "\nSTATION METADATA BUILD COMPLETE",
  "\n========================================"
)


print(
  metadata_coverage$summary
)