# ============================================================
# read_acorn_homogenised.R
#
# Read published homogenised ACORN-SAT daily Tmax and Tmin
# station CSV files from the Bureau TAR.GZ archives.
#
# Each station CSV contains:
#   - row 1: station metadata
#   - rows 2 onward: daily observations
#
# Temperature values are already in degrees Celsius.
# ============================================================


acorn_homogenised_members <- function(archive, variable = c("tmax", "tmin")) {
  variable <- match.arg(variable)
  
  if (!file.exists(archive)) {
    stop("Archive does not exist: ", archive, call. = FALSE)
  }
  
  members <- utils::untar(archive, list = TRUE)
  
  pattern <- paste0("^", variable, "[.][0-9]{6}[.]daily[.]csv$")
  
  station_members <- members[grepl(pattern, members)]
  
  if (length(station_members) == 0L) {
    stop("No ",
         variable,
         " station CSV files found in archive: ",
         archive,
         call. = FALSE)
  }
  
  sort(station_members)
}


extract_acorn_station_number <- function(member) {
  filename <- basename(member)
  
  station_number <- substr(filename, 6L, 11L)
  
  if (!grepl("^[0-9]{6}$", station_number)) {
    stop(
      "Could not extract six-digit station number from: ",
      member,
      ". Extracted value was: ",
      station_number,
      call. = FALSE
    )
  }
  
  station_number
}


read_acorn_homogenised_member <- function(archive,
                                          member,
                                          variable = c("tmax", "tmin"),
                                          source_version = "2.3.0") {
  variable <- match.arg(variable)
  
  expected_station <-
    extract_acorn_station_number(member)
  
  extraction_dir <- file.path(tempdir(),
                              paste0("acorn_homogenised_", variable, "_", expected_station))
  
  dir.create(extraction_dir,
             recursive = TRUE,
             showWarnings = FALSE)
  
  utils::untar(archive, files = member, exdir = extraction_dir)
  
  extracted_file <- file.path(extraction_dir, member)
  
  if (!file.exists(extracted_file)) {
    stop("Failed to extract archive member: ", member, call. = FALSE)
  }
  
  x <- data.table::fread(extracted_file,
                         na.strings = c("", "NA", "NaN"),
                         showProgress = FALSE)
  
  required_common <- c("date", "site number", "site name")
  
  value_column <- switch(variable, tmax = "maximum temperature (degC)", tmin = "minimum temperature (degC)")
  
  required <- c(required_common, value_column)
  
  missing_columns <- setdiff(required, names(x))
  
  if (length(missing_columns) > 0L) {
    stop(
      "File ",
      member,
      " is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Row 1 contains station metadata.
  metadata_rows <- x[is.na(date) &
                       !is.na(`site number`)]
  
  if (nrow(metadata_rows) != 1L) {
    stop(
      "Expected exactly one station metadata row in ",
      member,
      "; found ",
      nrow(metadata_rows),
      ".",
      call. = FALSE
    )
  }
  
  station_number <- sprintf("%06d", as.integer(metadata_rows[["site number"]][[1]]))
  
  station_name <- trimws(as.character(metadata_rows[["site name"]][[1]]))
  
  if (station_number != expected_station) {
    stop(
      "Station identifier mismatch in ",
      member,
      ". Filename gives ",
      expected_station,
      "; metadata gives ",
      station_number,
      ".",
      call. = FALSE
    )
  }
  
  # Keep only dated observations.
  observations <- x[!is.na(date)]
  
  observations[, date := as.Date(date)]
  
  data.table::setnames(observations, value_column, variable)
  
  observations[, `:=`(
    acorn_station_number = station_number,
    station_name = station_name,
    source_dataset =
      paste0("acorn_sat_homogenised_", variable),
    source_version = source_version,
    data_version = "homogenised",
    source_member = member
  )]
  
  result <- observations[, c(
    "acorn_station_number",
    "station_name",
    "date",
    variable,
    "source_dataset",
    "source_version",
    "data_version",
    "source_member"
  ), with = FALSE]
  
  data.table::setorder(result, acorn_station_number, date)
  
  duplicate_dates <- result[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  if (nrow(duplicate_dates) > 0L) {
    stop("Duplicated dates found in ", member, ".", call. = FALSE)
  }
  
  result
}


read_acorn_homogenised_archive <- function(archive,
                                           variable = c("tmax", "tmin"),
                                           source_version = "2.3.0") {
  variable <- match.arg(variable)
  
  members <- acorn_homogenised_members(archive = archive, variable = variable)
  
  total <- length(members)
  
  result <- vector(mode = "list", length = total)
  
  for (i in seq_along(members)) {
    message("[",
            i,
            "/",
            total,
            "] Reading homogenised ",
            variable,
            ": ",
            members[[i]])
    
    result[[i]] <- read_acorn_homogenised_member(
      archive = archive,
      member = members[[i]],
      variable = variable,
      source_version = source_version
    )
  }
  
  result <- data.table::rbindlist(result, use.names = TRUE, fill = TRUE)
  
  data.table::setorder(result, acorn_station_number, date)
  
  result
}


combine_acorn_homogenised <- function(tmax_data, tmin_data) {
  tmax_required <- c("acorn_station_number", "station_name", "date", "tmax")
  
  tmin_required <- c("acorn_station_number", "station_name", "date", "tmin")
  
  missing_tmax <- setdiff(tmax_required, names(tmax_data))
  
  missing_tmin <- setdiff(tmin_required, names(tmin_data))
  
  if (length(missing_tmax) > 0L) {
    stop("Tmax data are missing columns: ",
         paste(missing_tmax, collapse = ", "),
         call. = FALSE)
  }
  
  if (length(missing_tmin) > 0L) {
    stop("Tmin data are missing columns: ",
         paste(missing_tmin, collapse = ", "),
         call. = FALSE)
  }
  
  tmax_keep <- tmax_data[, .(
    acorn_station_number,
    station_name_tmax = station_name,
    date,
    tmax_hom = tmax,
    tmax_source_member = source_member
  )]
  
  tmin_keep <- tmin_data[, .(
    acorn_station_number,
    station_name_tmin = station_name,
    date,
    tmin_hom = tmin,
    tmin_source_member = source_member
  )]
  
  combined <- merge(
    tmax_keep,
    tmin_keep,
    by = c("acorn_station_number", "date"),
    all = TRUE
  )
  
  station_name_mismatch <- combined[!is.na(station_name_tmax) &
                                      !is.na(station_name_tmin) &
                                      station_name_tmax != station_name_tmin]
  
  if (nrow(station_name_mismatch) > 0L) {
    stop("Tmax/Tmin station-name mismatches were found.", call. = FALSE)
  }
  
  combined[, station_name :=
             data.table::fcoalesce(station_name_tmax, station_name_tmin)]
  
  combined[, `:=`(
    source_dataset =
      "acorn_sat_homogenised",
    
    source_version =
      ACORN_HOMOGENISED_VERSION,
    
    data_version =
      "homogenised"
  )]
  
  combined <- combined[, .(
    acorn_station_number,
    station_name,
    date,
    tmax_hom,
    tmin_hom,
    source_dataset,
    source_version,
    data_version,
    tmax_source_member,
    tmin_source_member
  )]
  
  data.table::setorder(combined, acorn_station_number, date)
  
  duplicate_dates <- combined[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  if (nrow(duplicate_dates) > 0L) {
    stop("Combined homogenised dataset contains duplicate dates.",
         call. = FALSE)
  }
  
  combined
}


validate_acorn_homogenised <- function(data) {
  required <- c("acorn_station_number",
                "station_name",
                "date",
                "tmax_hom",
                "tmin_hom")
  
  missing_columns <- setdiff(required, names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Homogenised data are missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  duplicate_dates <- data[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  tibble::tibble(
    rows = nrow(data),
    
    stations = data.table::uniqueN(data$acorn_station_number),
    
    first_date = min(data$date, na.rm = TRUE),
    
    last_date = max(data$date, na.rm = TRUE),
    
    missing_tmax = sum(is.na(data$tmax_hom)),
    
    missing_tmin = sum(is.na(data$tmin_hom)),
    
    duplicate_station_dates =
      nrow(duplicate_dates),
    
    implausible_tmax = sum(data$tmax_hom < -60 |
                             data$tmax_hom > 60, na.rm = TRUE),
    
    implausible_tmin = sum(data$tmin_hom < -60 |
                             data$tmin_hom > 60, na.rm = TRUE),
    
    tmin_above_tmax = sum(data$tmin_hom > data$tmax_hom, na.rm = TRUE)
  )
}