# ============================================================
# read_adam_stations.R
#
# Import, standardise and validate ADAM STATIONS metadata.
# ============================================================


# ============================================================
# Read station metadata
# ============================================================

read_adam_stations <- function(path) {
  
  if (!file.exists(path)) {
    stop(
      "ADAM station metadata file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  message(
    "Reading station metadata: ",
    basename(path)
  )
  
  x <- data.table::fread(
    path,
    na.strings = c(
      "",
      "NA",
      "NULL",
      "(null)"
    ),
    encoding = "Latin-1",
    showProgress = TRUE
  )
  
  data.table::setnames(
    x,
    tolower(names(x))
  )
  
  required <- c(
    "stn_num",
    "pr_name",
    "lat",
    "lon",
    "stn_ht",
    "state",
    "dist_id",
    "area_name",
    "bgn_date",
    "curr_bgn_date",
    "end_date"
  )
  
  missing_columns <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "ADAM station metadata are missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Station number
  # ----------------------------------------------------------
  
  station_integer <- suppressWarnings(
    as.integer(x$stn_num)
  )
  
  if (any(is.na(station_integer))) {
    stop(
      sum(is.na(station_integer)),
      " rows contain invalid station numbers.",
      call. = FALSE
    )
  }
  
  station_number <- sprintf(
    "%06d",
    station_integer
  )
  
  
  # ----------------------------------------------------------
  # ISO date parser
  # ----------------------------------------------------------
  
  parse_iso_date <- function(value) {
    
    value <- trimws(
      as.character(value)
    )
    
    value[value == ""] <- NA_character_
    
    parsed <- as.Date(
      value,
      format = "%Y-%m-%d"
    )
    
    invalid <-
      !is.na(value) &
      is.na(parsed)
    
    if (any(invalid)) {
      
      examples <- unique(
        value[invalid]
      )
      
      stop(
        "Unparseable metadata dates. Examples: ",
        paste(
          head(examples, 10L),
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    parsed
  }
  
  
  open_date <- parse_iso_date(
    x$bgn_date
  )
  
  current_begin_date <- parse_iso_date(
    x$curr_bgn_date
  )
  
  close_date <- parse_iso_date(
    x$end_date
  )
  
  
  # ----------------------------------------------------------
  # Optional-column helper
  # ----------------------------------------------------------
  
  optional_character <- function(column) {
    
    if (column %in% names(x)) {
      trimws(
        as.character(
          x[[column]]
        )
      )
    } else {
      rep(
        NA_character_,
        nrow(x)
      )
    }
  }
  
  
  optional_numeric <- function(column) {
    
    if (column %in% names(x)) {
      suppressWarnings(
        as.numeric(
          x[[column]]
        )
      )
    } else {
      rep(
        NA_real_,
        nrow(x)
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Canonical metadata table
  # ----------------------------------------------------------
  
  result <- data.table::data.table(
    
    station_number =
      station_number,
    
    station_name =
      trimws(
        as.character(
          x$pr_name
        )
      ),
    
    latitude =
      suppressWarnings(
        as.numeric(
          x$lat
        )
      ),
    
    longitude =
      suppressWarnings(
        as.numeric(
          x$lon
        )
      ),
    
    elevation_m =
      suppressWarnings(
        as.numeric(
          x$stn_ht
        )
      ),
    
    elevation_quality =
      optional_character(
        "stn_ht_qual"
      ),
    
    barometer_height_m =
      optional_numeric(
        "bar_ht"
      ),
    
    barometer_height_quality =
      optional_character(
        "bar_ht_qual"
      ),
    
    aerodrome_height_m =
      optional_numeric(
        "aero_ht"
      ),
    
    aerodrome_height_quality =
      optional_character(
        "aero_ht_qual"
      ),
    
    aviation_id =
      optional_character(
        "aero_id"
      ),
    
    wmo_number =
      optional_character(
        "wmo_num"
      ),
    
    state =
      trimws(
        as.character(
          x$state
        )
      ),
    
    rainfall_district =
      trimws(
        as.character(
          x$dist_id
        )
      ),
    
    area_name =
      trimws(
        as.character(
          x$area_name
        )
      ),
    
    open_date =
      open_date,
    
    current_begin_date =
      current_begin_date,
    
    close_date =
      close_date,
    
    timezone =
      optional_character(
        "tm_zone"
      ),
    
    time_difference =
      optional_numeric(
        "tm_dif"
      ),
    
    comment =
      optional_character(
        "cmt"
      ),
    
    source_dataset =
      "adam_stations",
    
    source_file =
      basename(path)
  )
  
  
  # ----------------------------------------------------------
  # Convert empty strings to NA
  # ----------------------------------------------------------
  
  character_columns <- names(result)[
    vapply(
      result,
      is.character,
      logical(1)
    )
  ]
  
  for (column in character_columns) {
    
    values <- result[[column]]
    
    empty <- !is.na(values) &
      trimws(values) == ""
    
    if (any(empty)) {
      data.table::set(
        result,
        i = which(empty),
        j = column,
        value = NA_character_
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Operational status
  # ----------------------------------------------------------
  
  result[
    ,
    station_open :=
      is.na(close_date)
  ]
  
  
  # ----------------------------------------------------------
  # Coordinate QC
  # ----------------------------------------------------------
  
  result[
    ,
    valid_coordinates :=
      !is.na(latitude) &
      !is.na(longitude) &
      latitude >= -90 &
      latitude <= 90 &
      longitude >= -180 &
      longitude <= 180
  ]
  
  # Broad Australian region diagnostic only.
  # Do not use this as an automatic exclusion criterion.
  
  result[
    ,
    australia_coordinate_check :=
      valid_coordinates &
      latitude >= -45 &
      latitude <= -9 &
      longitude >= 110 &
      longitude <= 155
  ]
  
  
  # ----------------------------------------------------------
  # Date QC
  # ----------------------------------------------------------
  
  result[
    ,
    close_before_open :=
      !is.na(open_date) &
      !is.na(close_date) &
      close_date < open_date
  ]
  
  result[
    ,
    begin_date_disagreement_days :=
      data.table::fifelse(
        !is.na(open_date) &
          !is.na(current_begin_date),
        as.numeric(
          current_begin_date -
            open_date
        ),
        NA_real_
      )
  ]
  
  
  data.table::setorder(
    result,
    station_number
  )
  
  result
}


# ============================================================
# Validate station metadata
# ============================================================

validate_adam_stations <- function(stations) {
  
  required <- c(
    "station_number",
    "station_name",
    "latitude",
    "longitude",
    "elevation_m",
    "open_date",
    "close_date",
    "station_open",
    "valid_coordinates",
    "close_before_open"
  )
  
  missing_columns <- setdiff(
    required,
    names(stations)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Station metadata are missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  duplicate_station_numbers <- stations[
    ,
    .N,
    by = station_number
  ][
    N > 1L
  ]
  
  data.table::data.table(
    
    rows =
      nrow(stations),
    
    unique_station_numbers =
      data.table::uniqueN(
        stations$station_number
      ),
    
    duplicate_station_numbers =
      nrow(
        duplicate_station_numbers
      ),
    
    missing_station_name =
      sum(
        is.na(stations$station_name)
      ),
    
    missing_latitude =
      sum(
        is.na(stations$latitude)
      ),
    
    missing_longitude =
      sum(
        is.na(stations$longitude)
      ),
    
    missing_elevation =
      sum(
        is.na(stations$elevation_m)
      ),
    
    invalid_coordinates =
      sum(
        !stations$valid_coordinates,
        na.rm = TRUE
      ),
    
    outside_main_australian_bounds =
      sum(
        stations$valid_coordinates &
          !stations$australia_coordinate_check,
        na.rm = TRUE
      ),
    
    missing_open_date =
      sum(
        is.na(stations$open_date)
      ),
    
    currently_open =
      sum(
        stations$station_open,
        na.rm = TRUE
      ),
    
    closed =
      sum(
        !stations$station_open,
        na.rm = TRUE
      ),
    
    close_before_open =
      sum(
        stations$close_before_open,
        na.rm = TRUE
      )
  )
}


# ============================================================
# Audit rainfall station metadata coverage
# ============================================================

audit_rainfall_station_metadata <- function(
    rainfall_station_year,
    stations) {
  
  rainfall_stations <- sort(
    unique(
      rainfall_station_year$station_number
    )
  )
  
  metadata_stations <- sort(
    unique(
      stations$station_number
    )
  )
  
  matched <- intersect(
    rainfall_stations,
    metadata_stations
  )
  
  unmatched_rainfall <- setdiff(
    rainfall_stations,
    metadata_stations
  )
  
  metadata_without_rainfall <- setdiff(
    metadata_stations,
    rainfall_stations
  )
  
  summary <- data.table::data.table(
    
    rainfall_stations =
      length(rainfall_stations),
    
    metadata_stations =
      length(metadata_stations),
    
    matched_rainfall_stations =
      length(matched),
    
    unmatched_rainfall_stations =
      length(unmatched_rainfall),
    
    rainfall_match_fraction =
      length(matched) /
      length(rainfall_stations),
    
    metadata_without_rainfall =
      length(metadata_without_rainfall)
  )
  
  list(
    
    summary =
      summary,
    
    unmatched_rainfall =
      data.table::data.table(
        station_number =
          unmatched_rainfall
      ),
    
    metadata_without_rainfall =
      data.table::data.table(
        station_number =
          metadata_without_rainfall
      )
  )
}


# ============================================================
# Join rainfall availability to station metadata
# ============================================================

join_rainfall_station_metadata <- function(
    rainfall_station_year,
    stations) {
  
  required_rainfall <- c(
    "station_number",
    "year"
  )
  
  required_stations <- c(
    "station_number",
    "station_name",
    "latitude",
    "longitude"
  )
  
  missing_rainfall <- setdiff(
    required_rainfall,
    names(rainfall_station_year)
  )
  
  missing_stations <- setdiff(
    required_stations,
    names(stations)
  )
  
  if (length(missing_rainfall) > 0L) {
    stop(
      "Rainfall station-year data are missing columns: ",
      paste(
        missing_rainfall,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (length(missing_stations) > 0L) {
    stop(
      "Station metadata are missing columns: ",
      paste(
        missing_stations,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  duplicate_metadata <- stations[
    ,
    .N,
    by = station_number
  ][
    N > 1L
  ]
  
  if (nrow(duplicate_metadata) > 0L) {
    stop(
      "Station metadata contain ",
      nrow(duplicate_metadata),
      " duplicated station numbers. ",
      "Resolve these before joining.",
      call. = FALSE
    )
  }
  
  result <- merge(
    rainfall_station_year,
    stations,
    by = "station_number",
    all.x = TRUE,
    sort = FALSE
  )
  
  result[
    ,
    metadata_available :=
      !is.na(station_name) |
      (
        !is.na(latitude) &
          !is.na(longitude)
      )
  ]
  
  data.table::setorder(
    result,
    station_number,
    year
  )
  
  result
}