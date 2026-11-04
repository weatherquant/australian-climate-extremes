# ============================================================
# convert_rainfall_to_parquet.R
#
# Sequentially convert ADAM PRCP_DAYS CSV files to compressed
# Parquet without loading the national archive simultaneously.
#
# Existing valid Parquet files can be skipped, making the
# conversion restartable after interruption.
# ============================================================


rainfall_parquet_name <- function(source_file) {
  
  sub(
    "[.]csv$",
    ".parquet",
    basename(source_file),
    ignore.case = TRUE
  )
}


validate_rainfall_chunk <- function(
    rainfall,
    source_file) {
  
  duplicate_keys <- rainfall[
    ,
    .N,
    by = .(
      station_number,
      date
    )
  ][
    N > 1L
  ]
  
  validation <- data.table::data.table(
    
    source_file =
      basename(source_file),
    
    rows =
      nrow(rainfall),
    
    stations =
      data.table::uniqueN(
        rainfall$station_number
      ),
    
    first_date =
      min(
        rainfall$date,
        na.rm = TRUE
      ),
    
    last_date =
      max(
        rainfall$date,
        na.rm = TRUE
      ),
    
    rainfall_missing =
      sum(
        is.na(rainfall$rainfall_mm)
      ),
    
    rainfall_zero =
      sum(
        rainfall$rainfall_mm == 0,
        na.rm = TRUE
      ),
    
    rainfall_positive =
      sum(
        rainfall$rainfall_mm > 0,
        na.rm = TRUE
      ),
    
    rainfall_negative =
      sum(
        rainfall$rainfall_mm < 0,
        na.rm = TRUE
      ),
    
    duplicate_station_dates =
      nrow(duplicate_keys),
    
    nonstandard_period =
      sum(
        !is.na(rainfall$rainfall_period) &
          rainfall$rainfall_period != 1
      ),
    
    multiday_count =
      sum(
        !is.na(rainfall$rainfall_day_count) &
          rainfall$rainfall_day_count > 1
      )
  )
  
  list(
    validation = validation,
    duplicates = duplicate_keys
  )
}


read_existing_parquet_validation <- function(
    parquet_file,
    source_file) {
  
  # Read only the fields required to rebuild the chunk audit.
  existing <- data.table::as.data.table(
    arrow::read_parquet(
      parquet_file,
      col_select = c(
        "station_number",
        "date",
        "rainfall_mm",
        "rainfall_period",
        "rainfall_day_count"
      )
    )
  )
  
  audit <- validate_rainfall_chunk(
    rainfall = existing,
    source_file = source_file
  )
  
  rm(existing)
  
  invisible(
    gc()
  )
  
  audit
}


convert_rainfall_file_to_parquet <- function(
    source_file,
    destination_dir,
    overwrite = FALSE) {
  
  if (length(source_file) != 1L) {
    
    stop(
      "Exactly one source file must be supplied.",
      call. = FALSE
    )
  }
  
  if (!file.exists(source_file)) {
    
    stop(
      "Source file does not exist: ",
      source_file,
      call. = FALSE
    )
  }
  
  dir.create(
    destination_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  output_name <- rainfall_parquet_name(
    source_file
  )
  
  output_file <- file.path(
    destination_dir,
    output_name
  )
  
  source_info <- file.info(
    source_file
  )
  
  if (
    file.exists(output_file) &&
    !overwrite
  ) {
    
    message(
      "Already converted, validating existing file: ",
      output_name
    )
    
    existing_audit <-
      read_existing_parquet_validation(
        parquet_file = output_file,
        source_file = source_file
      )
    
    output_info <- file.info(
      output_file
    )
    
    existing_audit$validation[
      ,
      `:=`(
        status = "existing",
        source_bytes = source_info$size,
        parquet_bytes = output_info$size,
        parquet_file = output_name
      )
    ]
    
    existing_audit$validation[
      ,
      compression_ratio :=
        parquet_bytes / source_bytes
    ]
    
    return(
      list(
        validation =
          existing_audit$validation,
        
        duplicates =
          existing_audit$duplicates,
        
        parquet_file =
          output_file,
        
        status =
          "existing"
      )
    )
  }
  
  message(
    "\n========================================"
  )
  
  message(
    "\nProcessing: ",
    basename(source_file)
  )
  
  message(
    "Source size: ",
    round(
      source_info$size / 1024^2,
      1
    ),
    " MB"
  )
  
  rainfall <- read_adam_rainfall_file(
    source_file
  )
  
  audit <- validate_rainfall_chunk(
    rainfall = rainfall,
    source_file = source_file
  )
  
  if (
    audit$validation$rainfall_negative[[1]] > 0L
  ) {
    
    warning(
      basename(source_file),
      " contains negative rainfall values. ",
      "They are retained in Parquet and flagged in validation.",
      call. = FALSE
    )
  }
  
  # Write to a temporary Parquet before replacing the final file.
  temporary_file <- paste0(
    output_file,
    ".part"
  )
  
  if (file.exists(temporary_file)) {
    unlink(temporary_file)
  }
  
  arrow::write_parquet(
    rainfall,
    temporary_file,
    compression = "zstd"
  )
  
  if (!file.exists(temporary_file)) {
    
    stop(
      "Temporary Parquet output was not created: ",
      temporary_file,
      call. = FALSE
    )
  }
  
  if (file.info(temporary_file)$size <= 0) {
    
    unlink(temporary_file)
    
    stop(
      "Temporary Parquet output is empty: ",
      temporary_file,
      call. = FALSE
    )
  }
  
  if (file.exists(output_file)) {
    unlink(output_file)
  }
  
  rename_success <- file.rename(
    temporary_file,
    output_file
  )
  
  if (!rename_success) {
    
    stop(
      "Could not rename temporary Parquet to: ",
      output_file,
      call. = FALSE
    )
  }
  
  output_info <- file.info(
    output_file
  )
  
  audit$validation[
    ,
    `:=`(
      status = "converted",
      source_bytes = source_info$size,
      parquet_bytes = output_info$size,
      parquet_file = output_name
    )
  ]
  
  audit$validation[
    ,
    compression_ratio :=
      parquet_bytes / source_bytes
  ]
  
  message(
    "Rows: ",
    format(
      nrow(rainfall),
      big.mark = ","
    )
  )
  
  message(
    "Stations: ",
    data.table::uniqueN(
      rainfall$station_number
    )
  )
  
  message(
    "Period: ",
    min(rainfall$date),
    " to ",
    max(rainfall$date)
  )
  
  message(
    "Parquet size: ",
    round(
      output_info$size / 1024^2,
      1
    ),
    " MB"
  )
  
  message(
    "Parquet/source ratio: ",
    round(
      output_info$size /
        source_info$size,
      3
    )
  )
  
  rm(rainfall)
  
  invisible(
    gc()
  )
  
  list(
    validation =
      audit$validation,
    
    duplicates =
      audit$duplicates,
    
    parquet_file =
      output_file,
    
    status =
      "converted"
  )
}


convert_all_rainfall_files <- function(
    source_files,
    destination_dir,
    audit_dir,
    overwrite = FALSE) {
  
  if (length(source_files) == 0L) {
    
    stop(
      "No rainfall source files were supplied.",
      call. = FALSE
    )
  }
  
  dir.create(
    destination_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    audit_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  results <- vector(
    mode = "list",
    length = length(source_files)
  )
  
  failures <- list()
  
  for (i in seq_along(source_files)) {
    
    message(
      "\n\n========================================",
      "\nFILE ",
      i,
      " OF ",
      length(source_files),
      "\n========================================"
    )
    
    result <- tryCatch(
      
      convert_rainfall_file_to_parquet(
        source_file =
          source_files[[i]],
        
        destination_dir =
          destination_dir,
        
        overwrite =
          overwrite
      ),
      
      error = function(e) {
        
        failure <- data.table::data.table(
          
          source_file =
            basename(
              source_files[[i]]
            ),
          
          error_message =
            conditionMessage(e),
          
          failed_at =
            Sys.time()
        )
        
        failures[[length(failures) + 1L]] <<-
          failure
        
        message(
          "FAILED: ",
          basename(source_files[[i]]),
          "\n",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    results[[i]] <- result
    
    # Save progress after every file.
    successful_results <- results[
      !vapply(
        results,
        is.null,
        logical(1)
      )
    ]
    
    if (length(successful_results) > 0L) {
      
      current_validation <-
        data.table::rbindlist(
          lapply(
            successful_results,
            `[[`,
            "validation"
          ),
          use.names = TRUE,
          fill = TRUE
        )
      
      data.table::fwrite(
        current_validation,
        file.path(
          audit_dir,
          "adam_prcp_days_conversion_validation.csv"
        )
      )
    }
    
    if (length(failures) > 0L) {
      
      current_failures <-
        data.table::rbindlist(
          failures,
          use.names = TRUE,
          fill = TRUE
        )
      
      data.table::fwrite(
        current_failures,
        file.path(
          audit_dir,
          "adam_prcp_days_conversion_failures.csv"
        )
      )
    }
    
    invisible(
      gc()
    )
  }
  
  successful_results <- results[
    !vapply(
      results,
      is.null,
      logical(1)
    )
  ]
  
  validation <- data.table::rbindlist(
    lapply(
      successful_results,
      `[[`,
      "validation"
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  duplicate_tables <- lapply(
    successful_results,
    `[[`,
    "duplicates"
  )
  
  duplicate_tables <- duplicate_tables[
    vapply(
      duplicate_tables,
      nrow,
      integer(1)
    ) > 0L
  ]
  
  if (length(duplicate_tables) > 0L) {
    
    duplicates <- data.table::rbindlist(
      duplicate_tables,
      use.names = TRUE,
      fill = TRUE
    )
    
  } else {
    
    duplicates <- data.table::data.table(
      station_number = character(),
      date = as.Date(character()),
      N = integer()
    )
  }
  
  data.table::fwrite(
    validation,
    file.path(
      audit_dir,
      "adam_prcp_days_conversion_validation.csv"
    )
  )
  
  data.table::fwrite(
    duplicates,
    file.path(
      audit_dir,
      "adam_prcp_days_duplicate_station_dates.csv"
    )
  )
  
  failure_table <- if (
    length(failures) > 0L
  ) {
    
    data.table::rbindlist(
      failures,
      use.names = TRUE,
      fill = TRUE
    )
    
  } else {
    
    data.table::data.table(
      source_file = character(),
      error_message = character(),
      failed_at = as.POSIXct(character())
    )
  }
  
  data.table::fwrite(
    failure_table,
    file.path(
      audit_dir,
      "adam_prcp_days_conversion_failures.csv"
    )
  )
  
  list(
    validation = validation,
    duplicates = duplicates,
    failures = failure_table,
    results = results
  )
}