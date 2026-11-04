# ============================================================
# ADAM PRCP_DAYS CSV TO PARQUET
# ============================================================

library(data.table)
library(arrow)
library(yaml)

source(
  "R/ingest/read_station_rainfall.R"
)

source(
  "R/processing/convert_rainfall_to_parquet.R"
)


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

cfg <- yaml::read_yaml(
  "config.local.yml"
)

rainfall_source_dir <- file.path(
  cfg$paths$raw_root,
  cfg$paths$station$rainfall
)

if (!dir.exists(rainfall_source_dir)) {
  
  stop(
    "Rainfall source directory does not exist: ",
    rainfall_source_dir
  )
}


# ------------------------------------------------------------
# Discover source files
# ------------------------------------------------------------

rainfall_files <- list.files(
  rainfall_source_dir,
  pattern = "^adam_prcp_days_.*[.]csv$",
  full.names = TRUE
)

if (length(rainfall_files) == 0L) {
  
  stop(
    "No ADAM rainfall CSV files found in: ",
    rainfall_source_dir
  )
}

rainfall_files <- sort(
  rainfall_files
)

message(
  "Discovered ",
  length(rainfall_files),
  " rainfall CSV files."
)

print(
  basename(rainfall_files)
)


# ------------------------------------------------------------
# Exclude pre-1900 until it is re-exported using an explicit
# four-digit ISO date.
# ------------------------------------------------------------

production_files <- rainfall_files[
  !grepl(
    "pre1900",
    basename(rainfall_files),
    ignore.case = TRUE
  )
]

message(
  "Processing ",
  length(production_files),
  " files covering 1900 onwards."
)


# ------------------------------------------------------------
# Output locations
# ------------------------------------------------------------

rainfall_parquet_dir <- file.path(
  "data-derived",
  "rainfall",
  "adam_prcp_days"
)

rainfall_audit_dir <- file.path(
  "data-derived",
  "rainfall",
  "audit"
)

dir.create(
  rainfall_parquet_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  rainfall_audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Source-file inventory and checksums
# ------------------------------------------------------------

rainfall_source_inventory <-
  data.table::rbindlist(
    lapply(
      rainfall_files,
      function(path) {
        
        info <- file.info(path)
        
        data.table::data.table(
          
          filename =
            basename(path),
          
          full_path =
            normalizePath(
              path,
              winslash = "/",
              mustWork = TRUE
            ),
          
          bytes =
            info$size,
          
          size_mb =
            info$size / 1024^2,
          
          modified =
            info$mtime,
          
          md5 =
            unname(
              tools::md5sum(path)
            )
        )
      }
    )
  )

data.table::fwrite(
  rainfall_source_inventory,
  file.path(
    rainfall_audit_dir,
    "adam_prcp_days_source_inventory.csv"
  )
)

print(
  rainfall_source_inventory[
    ,
    .(
      files = .N,
      total_gb =
        sum(bytes) / 1024^3
    )
  ]
)


# ------------------------------------------------------------
# Sequential conversion
#
# overwrite = FALSE means:
#   existing Parquet files are retained and validated;
#   missing Parquet files are created.
# ------------------------------------------------------------

rainfall_conversion <-
  convert_all_rainfall_files(
    
    source_files =
      production_files,
    
    destination_dir =
      rainfall_parquet_dir,
    
    audit_dir =
      rainfall_audit_dir,
    
    overwrite =
      FALSE
  )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print(
  rainfall_conversion$validation
)

print(
  rainfall_conversion$failures
)

conversion_summary <-
  rainfall_conversion$validation[
    ,
    .(
      source_files = .N,
      
      total_rows =
        sum(rows),
      
      source_gb =
        sum(source_bytes) /
        1024^3,
      
      parquet_gb =
        sum(parquet_bytes) /
        1024^3,
      
      compression_ratio =
        sum(parquet_bytes) /
        sum(source_bytes),
      
      duplicate_station_dates =
        sum(
          duplicate_station_dates
        ),
      
      negative_rainfall =
        sum(
          rainfall_negative
        ),
      
      missing_rainfall =
        sum(
          rainfall_missing
        ),
      
      zero_rainfall =
        sum(
          rainfall_zero
        ),
      
      positive_rainfall =
        sum(
          rainfall_positive
        )
    )
  ]

print(
  conversion_summary
)


# ------------------------------------------------------------
# Check Parquet inventory
# ------------------------------------------------------------

parquet_files <- list.files(
  rainfall_parquet_dir,
  pattern = "[.]parquet$",
  full.names = TRUE
)

message(
  "Parquet files present: ",
  length(parquet_files)
)

print(
  basename(parquet_files)
)


# ------------------------------------------------------------
# Open as an Arrow dataset
#
# This does not load all rows into R memory.
# ------------------------------------------------------------

rainfall_ds <- arrow::open_dataset(
  rainfall_parquet_dir,
  format = "parquet"
)

rainfall_ds

# ------------------------------------------------------------
# Save rainfall conversion summary
# ------------------------------------------------------------

data.table::fwrite(
  conversion_summary,
  file.path(
    rainfall_audit_dir,
    "adam_prcp_days_conversion_summary.csv"
  )
)


# ------------------------------------------------------------
# Get Git commit if Git is visible to the R process
# ------------------------------------------------------------

git_commit <- tryCatch(
  {
    result <- system2(
      "git",
      args = c(
        "rev-parse",
        "--short",
        "HEAD"
      ),
      stdout = TRUE,
      stderr = FALSE
    )
    
    if (length(result) == 0L) {
      NA_character_
    } else {
      result[[1]]
    }
  },
  error = function(e) {
    NA_character_
  }
)


# ------------------------------------------------------------
# Create pipeline provenance record
# ------------------------------------------------------------

rainfall_pipeline_record <- data.table::data.table(
  
  created_at = Sys.time(),
  
  git_commit = git_commit,
  
  source_files =
    conversion_summary$source_files,
  
  total_rows =
    conversion_summary$total_rows,
  
  source_gb =
    conversion_summary$source_gb,
  
  parquet_gb =
    conversion_summary$parquet_gb,
  
  compression_ratio =
    conversion_summary$compression_ratio,
  
  first_date =
    min(
      rainfall_conversion$validation$first_date
    ),
  
  last_date =
    max(
      rainfall_conversion$validation$last_date
    )
)


# ------------------------------------------------------------
# Save provenance record
# ------------------------------------------------------------

data.table::fwrite(
  rainfall_pipeline_record,
  file.path(
    rainfall_audit_dir,
    "adam_prcp_days_pipeline_record.csv"
  )
)


rainfall_pipeline_record
