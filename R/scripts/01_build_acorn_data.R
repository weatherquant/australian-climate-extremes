# ============================================================
# 01_build_acorn_data.R
#
# Build ACORN-SAT v2.3 raw and homogenised analytical datasets
#
# Produces:
#   - cleaned raw source-station data
#   - raw ACORN primary-site composite
#   - homogenised Tmax and Tmin
#   - raw-versus-homogenised comparison dataset
#   - QC and audit outputs
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(data.table)
library(arrow)
library(yaml)


# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

source("R/acquire/checksums.R")
source("R/acquire/download_file.R")
source("R/acquire/provenance.R")

source("R/acquire/acquire_acorn.R")
source("R/acquire/acquire_acorn_homogenised.R")

source("R/ingest/read_acorn_data.R")
source("R/ingest/read_acorn_primary_sites.R")
source("R/ingest/read_acorn_homogenised.R")

source("R/processing/build_acorn_raw_composite.R")
source("R/processing/match_acorn_raw_homogenised.R")


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

cfg <- yaml::read_yaml(
  "config.local.yml"
)


# ============================================================
# 1. ACORN RAW SOURCE
# ============================================================

message(
  "\n========================================",
  "\nACORN RAW DATA",
  "\n========================================"
)


# ------------------------------------------------------------
# Acquire raw ACORN bundle
# ------------------------------------------------------------

acorn_download <-
  acquire_acorn_raw(
    cfg
  )

acorn_zip <-
  acorn_download$destination[[1]]


if (!file.exists(acorn_zip)) {
  
  stop(
    "ACORN raw archive does not exist: ",
    acorn_zip
  )
}


# ------------------------------------------------------------
# Archive provenance
# ------------------------------------------------------------

acorn_provenance <-
  create_provenance_record(
    
    dataset_id =
      "acorn_raw",
    
    source_url =
      paste0(
        ACORN_BASE_URL,
        ACORN_RAW_BUNDLE
      ),
    
    source_version =
      ACORN_VERSION,
    
    path =
      acorn_zip,
    
    status =
      acorn_download$status[[1]]
  )


# ------------------------------------------------------------
# Read primary-site mapping
# ------------------------------------------------------------

acorn_primary_sites <-
  read_acorn_primary_sites(
    acorn_zip
  )

acorn_primary_sites_long_data <-
  acorn_primary_sites_long(
    acorn_primary_sites
  )


# ------------------------------------------------------------
# Read all raw station observations
# ------------------------------------------------------------

acorn_raw_all <-
  read_acorn_raw_all(
    zip_file = acorn_zip
  )


# ------------------------------------------------------------
# Validate original raw archive
# ------------------------------------------------------------

validation_raw <-
  validate_acorn_raw(
    acorn_raw_all
  )

print(
  validation_raw,
  width = Inf
)


# ------------------------------------------------------------
# Resolve duplicate source records
#
# Known v2.3 behaviour:
# conflicting raw duplicate station-dates are excluded and
# retained separately in the audit output.
# ------------------------------------------------------------

duplicate_resolution <-
  resolve_acorn_duplicates_for_composite(
    acorn_raw_all
  )

acorn_raw_clean <-
  duplicate_resolution$data

acorn_conflict_audit <-
  duplicate_resolution$conflict_audit

acorn_duplicate_summary <-
  duplicate_resolution$duplicate_summary


# ------------------------------------------------------------
# Validate cleaned raw data
# ------------------------------------------------------------

validation_clean <-
  validate_acorn_raw(
    acorn_raw_clean
  )

print(
  validation_clean,
  width = Inf
)


if (
  validation_clean$duplicate_station_dates[[1]] != 0L
) {
  
  stop(
    "Duplicate station-dates remain in cleaned ACORN raw data."
  )
}


# ------------------------------------------------------------
# Confirm all required primary stations are present
# ------------------------------------------------------------

raw_station_numbers <-
  unique(
    acorn_raw_clean$station_number
  )

required_primary_sites <-
  unique(
    acorn_primary_sites_long_data$
      source_station_number
  )

missing_primary_sites <-
  setdiff(
    required_primary_sites,
    raw_station_numbers
  )


if (length(missing_primary_sites) > 0L) {
  
  stop(
    "Required ACORN primary source stations are missing: ",
    paste(
      missing_primary_sites,
      collapse = ", "
    )
  )
}


# ============================================================
# 2. BUILD RAW ACORN COMPOSITE
# ============================================================

message(
  "\n========================================",
  "\nBUILD RAW ACORN COMPOSITE",
  "\n========================================"
)


acorn_raw_composite <-
  build_acorn_raw_composite(
    
    raw_data =
      acorn_raw_clean,
    
    primary_sites_long =
      acorn_primary_sites_long_data
  )


# ------------------------------------------------------------
# Composite duplicate check
# ------------------------------------------------------------

composite_duplicates <-
  acorn_raw_composite[
    ,
    .N,
    by = .(
      acorn_station_number,
      date
    )
  ][
    N > 1L
  ]


if (nrow(composite_duplicates) > 0L) {
  
  stop(
    "Raw ACORN composite contains duplicate station-dates."
  )
}


# ------------------------------------------------------------
# Composite boundary audit
# ------------------------------------------------------------

composite_audit <-
  validate_acorn_raw_composite(
    
    composite =
      acorn_raw_composite,
    
    primary_sites_long =
      acorn_primary_sites_long_data
  )


boundary_exceptions <-
  composite_audit[
    first_date_matches == FALSE |
      last_date_matches == FALSE
  ]


# ============================================================
# 3. SAVE RAW ACORN PRODUCTS
# ============================================================

derived_dir <- file.path(
  "data-derived",
  "acorn"
)

audit_dir <- file.path(
  derived_dir,
  "audit"
)


dir.create(
  derived_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


arrow::write_parquet(
  
  acorn_raw_clean,
  
  file.path(
    derived_dir,
    "acorn_raw_source_stations_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_primary_sites_long_data,
  
  file.path(
    derived_dir,
    "acorn_primary_sites_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_raw_composite,
  
  file.path(
    derived_dir,
    "acorn_raw_composite_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_conflict_audit,
  
  file.path(
    derived_dir,
    "acorn_raw_conflict_audit_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  composite_audit,
  
  file.path(
    derived_dir,
    "acorn_raw_composite_audit_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  boundary_exceptions,
  
  file.path(
    derived_dir,
    "acorn_raw_composite_boundary_exceptions_v2.3.parquet"
  ),
  
  compression = "zstd"
)


# ------------------------------------------------------------
# Audit CSVs
# ------------------------------------------------------------

data.table::fwrite(
  
  acorn_conflict_audit,
  
  file.path(
    audit_dir,
    "acorn_v2.3_excluded_conflicting_raw_records.csv"
  )
)


data.table::fwrite(
  
  boundary_exceptions,
  
  file.path(
    audit_dir,
    "acorn_v2.3_boundary_exceptions.csv"
  )
)


# ============================================================
# 4. HOMOGENISED ACORN
# ============================================================

message(
  "\n========================================",
  "\nACORN HOMOGENISED DATA",
  "\n========================================"
)


# ------------------------------------------------------------
# Acquire published homogenised Tmax/Tmin
# ------------------------------------------------------------

acorn_hom_download <-
  acquire_acorn_homogenised(
    cfg
  )


tmax_archive <-
  acorn_hom_download[
    variable == "tmax",
    destination
  ][[1]]


tmin_archive <-
  acorn_hom_download[
    variable == "tmin",
    destination
  ][[1]]


if (!file.exists(tmax_archive)) {
  stop(
    "Homogenised Tmax archive is missing."
  )
}


if (!file.exists(tmin_archive)) {
  stop(
    "Homogenised Tmin archive is missing."
  )
}


# ------------------------------------------------------------
# Read Tmax
# ------------------------------------------------------------

acorn_hom_tmax <-
  read_acorn_homogenised_archive(
    
    archive =
      tmax_archive,
    
    variable =
      "tmax",
    
    source_version =
      ACORN_HOMOGENISED_VERSION
  )


# ------------------------------------------------------------
# Read Tmin
# ------------------------------------------------------------

acorn_hom_tmin <-
  read_acorn_homogenised_archive(
    
    archive =
      tmin_archive,
    
    variable =
      "tmin",
    
    source_version =
      ACORN_HOMOGENISED_VERSION
  )


# ------------------------------------------------------------
# Combine Tmax and Tmin
# ------------------------------------------------------------

acorn_homogenised <-
  combine_acorn_homogenised(
    
    tmax_data =
      acorn_hom_tmax,
    
    tmin_data =
      acorn_hom_tmin
  )


# ------------------------------------------------------------
# Validate homogenised data
# ------------------------------------------------------------

validation_hom <-
  validate_acorn_homogenised(
    acorn_homogenised
  )


print(
  validation_hom,
  width = Inf
)


# ------------------------------------------------------------
# Confirm matching station networks
# ------------------------------------------------------------

raw_stations <-
  sort(
    unique(
      acorn_raw_composite$
        acorn_station_number
    )
  )


hom_stations <-
  sort(
    unique(
      acorn_homogenised$
        acorn_station_number
    )
  )


raw_only_stations <-
  setdiff(
    raw_stations,
    hom_stations
  )


hom_only_stations <-
  setdiff(
    hom_stations,
    raw_stations
  )


if (
  length(raw_only_stations) > 0L ||
  length(hom_only_stations) > 0L
) {
  
  stop(
    "Raw and homogenised ACORN station networks differ."
  )
}


# ============================================================
# 5. SAVE HOMOGENISED PRODUCTS
# ============================================================

hom_derived_dir <- file.path(
  derived_dir,
  "homogenised"
)


dir.create(
  hom_derived_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


arrow::write_parquet(
  
  acorn_hom_tmax,
  
  file.path(
    hom_derived_dir,
    "acorn_homogenised_tmax_v2.3.0.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_hom_tmin,
  
  file.path(
    hom_derived_dir,
    "acorn_homogenised_tmin_v2.3.0.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_homogenised,
  
  file.path(
    hom_derived_dir,
    "acorn_homogenised_daily_v2.3.0.parquet"
  ),
  
  compression = "zstd"
)


# ============================================================
# 6. RAW VERSUS HOMOGENISED COMPARISON
# ============================================================

message(
  "\n========================================",
  "\nMATCH RAW AND HOMOGENISED ACORN",
  "\n========================================"
)


acorn_temperature_comparison <-
  match_acorn_raw_homogenised(
    
    raw_composite =
      acorn_raw_composite,
    
    homogenised =
      acorn_homogenised
  )


matching_audit <-
  audit_acorn_matching(
    acorn_temperature_comparison
  )


acorn_adjustment_summary <-
  summarise_acorn_adjustments(
    acorn_temperature_comparison
  )


print(
  matching_audit
)


# ============================================================
# 7. SAVE COMPARISON PRODUCTS
# ============================================================

comparison_dir <- file.path(
  derived_dir,
  "comparison"
)


dir.create(
  comparison_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


arrow::write_parquet(
  
  acorn_temperature_comparison,
  
  file.path(
    comparison_dir,
    "acorn_raw_homogenised_comparison_v2.3.parquet"
  ),
  
  compression = "zstd"
)


arrow::write_parquet(
  
  acorn_adjustment_summary,
  
  file.path(
    comparison_dir,
    "acorn_adjustment_summary_v2.3.parquet"
  ),
  
  compression = "zstd"
)


data.table::fwrite(
  
  matching_audit,
  
  file.path(
    audit_dir,
    "acorn_raw_homogenised_matching_audit_v2.3.csv"
  )
)


# ============================================================
# 8. FINAL SUMMARY
# ============================================================

message(
  "\n========================================",
  "\nACORN BUILD COMPLETE",
  "\n========================================"
)


build_summary <- data.table::data.table(
  
  raw_source_rows =
    nrow(acorn_raw_all),
  
  clean_raw_rows =
    nrow(acorn_raw_clean),
  
  raw_source_stations =
    data.table::uniqueN(
      acorn_raw_clean$station_number
    ),
  
  raw_composite_rows =
    nrow(acorn_raw_composite),
  
  acorn_locations =
    data.table::uniqueN(
      acorn_raw_composite$
        acorn_station_number
    ),
  
  homogenised_rows =
    nrow(acorn_homogenised),
  
  comparison_rows =
    nrow(acorn_temperature_comparison),
  
  first_date =
    min(
      acorn_raw_composite$date
    ),
  
  last_date =
    max(
      acorn_raw_composite$date
    ),
  
  excluded_conflict_rows =
    nrow(acorn_conflict_audit),
  
  boundary_exceptions =
    nrow(boundary_exceptions)
)


print(
  build_summary
)


data.table::fwrite(
  
  build_summary,
  
  file.path(
    audit_dir,
    "acorn_v2.3_build_summary.csv"
  )
)
