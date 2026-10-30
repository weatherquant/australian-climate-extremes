# ============================================================
# read_acorn_data.R
#
# Read and standardise ACORN-SAT v2.3 raw temperature data.
#
# Source format:
#
#   NNNNNNYYYYMMDD  XXX  NNN
#
# NNNNNN   = station number
# YYYYMMDD = observation date
# XXX      = Tmax in tenths of degrees Celsius
# NNN      = Tmin in tenths of degrees Celsius
# -999     = missing value
#
# The raw files are quality controlled but not adjusted.
# ============================================================


acorn_raw_members <- function(zip_file) {
  if (!file.exists(zip_file)) {
    stop("ACORN ZIP file does not exist: ", zip_file, call. = FALSE)
  }
  
  inventory <- unzip(zip_file, list = TRUE)
  
  members <- inventory$Name[grepl("/raw-data/hqnew[0-9]{6}$", inventory$Name)]
  
  if (length(members) == 0L) {
    stop("No ACORN raw station members were found in: ", zip_file, call. = FALSE)
  }
  
  members
}


read_acorn_raw_member <- function(zip_file, member, source_version = "2.3") {
  if (!file.exists(zip_file)) {
    stop("ACORN ZIP file does not exist: ", zip_file, call. = FALSE)
  }
  
  connection <- unz(zip_file, member, open = "rt")
  
  on.exit(close(connection), add = TRUE)
  
  lines <- readLines(connection, warn = FALSE)
  
  lines <- trimws(lines)
  
  lines <- lines[nzchar(lines)]
  
  if (length(lines) == 0L) {
    return(
      data.table::data.table(
        station_number = character(),
        date = as.Date(character()),
        tmax = numeric(),
        tmin = numeric(),
        tmax_missing = logical(),
        tmin_missing = logical(),
        source_dataset = character(),
        source_version = character(),
        data_version = character(),
        source_member = character()
      )
    )
  }
  
  parsed <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    col.names = c("station_date", "tmax_tenths", "tmin_tenths"),
    colClasses = c("character", "integer", "integer"),
    showProgress = FALSE
  )
  
  if (ncol(parsed) != 3L) {
    stop("Unexpected number of columns in ACORN member: ",
         member,
         call. = FALSE)
  }
  
  parsed[, station_number :=
           substr(station_date, 1L, 6L)]
  
  parsed[, date_text :=
           substr(station_date, 7L, 14L)]
  
  parsed[, date :=
           as.Date(date_text, format = "%Y%m%d")]
  
  parsed[, tmax_missing :=
           is.na(tmax_tenths) |
           tmax_tenths == -999L]
  
  parsed[, tmin_missing :=
           is.na(tmin_tenths) |
           tmin_tenths == -999L]
  
  parsed[, tmax :=
           data.table::fifelse(tmax_missing, NA_real_, tmax_tenths / 10)]
  
  parsed[, tmin :=
           data.table::fifelse(tmin_missing, NA_real_, tmin_tenths / 10)]
  
  parsed[, `:=`(
    source_dataset = "acorn_sat_raw",
    source_version = source_version,
    data_version = "quality_controlled_unadjusted",
    source_member = member
  )]
  
  expected_station <- sub("^hqnew", "", basename(member))
  
  stations_found <- unique(parsed$station_number)
  
  stations_found <- stations_found[!is.na(stations_found)]
  
  if (length(stations_found) != 1L ||
      stations_found[[1]] != expected_station) {
    stop(
      "Station identifier mismatch in ",
      member,
      ". Expected ",
      expected_station,
      "; found ",
      paste(stations_found, collapse = ", "),
      call. = FALSE
    )
  }
  
  invalid_station_date <- nchar(parsed$station_date) != 14L
  
  if (any(invalid_station_date)) {
    stop("Invalid station/date identifier encountered in ",
         member,
         call. = FALSE)
  }
  
  if (any(is.na(parsed$date))) {
    stop("One or more dates could not be parsed in ", member, call. = FALSE)
  }
  
  data.table::setorder(parsed, station_number, date)
  
  parsed[, .(
    station_number,
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


read_acorn_raw_sample <- function(zip_file,
                                  n_stations = 3L,
                                  source_version = "2.3") {
  members <- acorn_raw_members(zip_file)
  
  n_stations <- min(as.integer(n_stations), length(members))
  
  result <- lapply(members[seq_len(n_stations)], function(member) {
    message("Reading ACORN raw member: ", basename(member))
    
    read_acorn_raw_member(zip_file = zip_file,
                          member = member,
                          source_version = source_version)
  })
  
  data.table::rbindlist(result, use.names = TRUE, fill = TRUE)
}


read_acorn_raw_all <- function(zip_file, source_version = "2.3") {
  members <- acorn_raw_members(zip_file)
  
  total <- length(members)
  
  result <- vector(mode = "list", length = total)
  
  for (i in seq_along(members)) {
    message("[", i, "/", total, "] Reading ", basename(members[[i]]))
    
    result[[i]] <- read_acorn_raw_member(zip_file = zip_file,
                                         member = members[[i]],
                                         source_version = source_version)
  }
  
  result <- data.table::rbindlist(result, use.names = TRUE, fill = TRUE)
  
  data.table::setorder(result, station_number, date)
  
  result
}


validate_acorn_raw <- function(data) {
  required <- c(
    "station_number",
    "date",
    "tmax",
    "tmin",
    "source_dataset",
    "source_version",
    "data_version",
    "source_member"
  )
  
  missing_columns <- setdiff(required, names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      "ACORN data are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  duplicate_rows <- data[, .N, by = .(station_number, date)][N > 1L]
  
  tibble::tibble(
    rows = nrow(data),
    stations = data.table::uniqueN(data$station_number),
    first_date = min(data$date, na.rm = TRUE),
    last_date = max(data$date, na.rm = TRUE),
    missing_tmax = sum(is.na(data$tmax)),
    missing_tmin = sum(is.na(data$tmin)),
    duplicate_station_dates = nrow(duplicate_rows),
    implausible_tmax = sum(data$tmax < -60 |
                             data$tmax > 60, na.rm = TRUE),
    implausible_tmin = sum(data$tmin < -60 |
                             data$tmin > 60, na.rm = TRUE),
    tmin_above_tmax = sum(data$tmin > data$tmax, na.rm = TRUE)
  )
}

resolve_acorn_raw_duplicates <- function(data) {
  required <- c("station_number",
                "date",
                "tmax",
                "tmin",
                "tmax_missing",
                "tmin_missing")
  
  missing_columns <- setdiff(required, names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Data are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  x <- data.table::copy(data)
  
  # Identify duplicated keys without aggregating every row.
  duplicate_indicator <- duplicated(x, by = c("station_number", "date")) |
    duplicated(x,
               by = c("station_number", "date"),
               fromLast = TRUE)
  
  unique_records <- x[!duplicate_indicator]
  
  duplicate_records <- x[duplicate_indicator]
  
  if (nrow(duplicate_records) == 0L) {
    message("No duplicated station-date records found.")
    return(x)
  }
  
  duplicate_audit <- duplicate_records[, .(
    n_rows = .N,
    
    n_tmax_nonmissing_values =
      data.table::uniqueN(tmax[!is.na(tmax)]),
    
    n_tmin_nonmissing_values =
      data.table::uniqueN(tmin[!is.na(tmin)])
  ), by = .(station_number, date)]
  
  genuine_conflicts <- duplicate_audit[n_tmax_nonmissing_values > 1L |
                                         n_tmin_nonmissing_values > 1L]
  
  if (nrow(genuine_conflicts) > 0L) {
    stop(
      "There are ",
      nrow(genuine_conflicts),
      " duplicated station-date keys containing genuinely ",
      "different non-missing temperature values. Inspect ",
      "true_conflict_records before continuing.",
      call. = FALSE
    )
  }
  
  # Consolidate exact duplicates and complementary records.
  resolved_duplicates <- duplicate_records[, {
    tmax_available <- tmax[!is.na(tmax)]
    tmin_available <- tmin[!is.na(tmin)]
    
    resolved_tmax <- if (length(tmax_available) == 0L) {
      NA_real_
    } else {
      tmax_available[[1]]
    }
    
    resolved_tmin <- if (length(tmin_available) == 0L) {
      NA_real_
    } else {
      tmin_available[[1]]
    }
    
    representative <- .SD[1L]
    
    representative[, `:=`(
      tmax = resolved_tmax,
      tmin = resolved_tmin,
      tmax_missing = is.na(resolved_tmax),
      tmin_missing = is.na(resolved_tmin)
    )]
    
    representative
  }, by = .(station_number, date)]
  
  result <- data.table::rbindlist(
    list(unique_records, resolved_duplicates),
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(result, station_number, date)
  
  remaining_duplicates <- result[, .N, by = .(station_number, date)][N > 1L]
  
  if (nrow(remaining_duplicates) > 0L) {
    stop("Duplicate resolution failed: duplicated keys remain.",
         call. = FALSE)
  }
  
  message(
    "Resolved ",
    nrow(duplicate_audit),
    " duplicated station-date keys and removed ",
    nrow(x) - nrow(result),
    " repeated rows."
  )
  
  result
}