# ============================================================
# read_station_rainfall.R
#
# Read and standardise ADAM PRCP_DAYS CSV exports.
#
# Source filenames:
#   adam_prcp_days_1900_1904.csv
#   ...
#   adam_prcp_days_2025_2027.csv
#
# Important:
#   blank PRCP = missing observation
#   PRCP = 0   = observed zero rainfall
#
# The extraction filename is used to resolve ambiguous
# two-digit years in the Oracle CSV export.
# ============================================================

rainfall_period_from_filename <- function(path) {
  
  filename <- basename(path)
  
  match <- regexec(
    "^adam_prcp_days_([0-9]{4})_([0-9]{4})[.]csv$",
    filename
  )
  
  parts <- regmatches(
    filename,
    match
  )[[1]]
  
  if (length(parts) != 3L) {
    
    stop(
      "Could not determine extraction period from filename: ",
      filename,
      call. = FALSE
    )
  }
  
  start_year <- as.integer(
    parts[[2]]
  )
  
  end_year <- as.integer(
    parts[[3]]
  )
  
  if (
    is.na(start_year) ||
    is.na(end_year) ||
    end_year < start_year
  ) {
    
    stop(
      "Invalid extraction period in filename: ",
      filename,
      call. = FALSE
    )
  }
  
  list(
    start_year = start_year,
    end_year = end_year
  )
}

parse_adam_rainfall_dates <- function(
    date_text,
    start_year,
    end_year) {
  
  date_text <- trimws(
    as.character(date_text)
  )
  
  date_text[
    date_text == ""
  ] <- NA_character_
  
  if (any(is.na(date_text))) {
    
    stop(
      "Missing source dates were found.",
      call. = FALSE
    )
  }
  
  # Supports both:
  #   01/JAN/10
  #   1-Jan-10
  #
  # The final two characters are the Oracle two-digit year.
  
  two_digit_year <- suppressWarnings(
    as.integer(
      substr(
        date_text,
        nchar(date_text) - 1L,
        nchar(date_text)
      )
    )
  )
  
  if (any(is.na(two_digit_year))) {
    
    stop(
      "One or more two-digit years could not be extracted.",
      call. = FALSE
    )
  }
  
  start_century <-
    floor(start_year / 100) * 100L
  
  full_year <-
    start_century + two_digit_year
  
  # Allows a future chunk that crosses a century boundary.
  before_start <- full_year < start_year
  
  full_year[
    before_start
  ] <- full_year[
    before_start
  ] + 100L
  
  day_month <- substr(
    date_text,
    1L,
    nchar(date_text) - 2L
  )
  
  reconstructed <- paste0(
    day_month,
    sprintf(
      "%04d",
      full_year
    )
  )
  
  parsed_date <- as.Date(
    reconstructed,
    format = "%d/%b/%Y"
  )
  
  if (any(is.na(parsed_date))) {
    
    examples <- unique(
      date_text[
        is.na(parsed_date)
      ]
    )
    
    stop(
      "Could not parse ",
      sum(is.na(parsed_date)),
      " rainfall dates. Examples: ",
      paste(
        head(examples, 10L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  parsed_year <- as.integer(
    format(
      parsed_date,
      "%Y"
    )
  )
  
  outside_period <-
    parsed_year < start_year |
    parsed_year > end_year
  
  if (any(outside_period)) {
    
    stop(
      sum(outside_period),
      " dates fall outside the expected ",
      start_year,
      "-",
      end_year,
      " extraction period.",
      call. = FALSE
    )
  }
  
  parsed_date
}


read_adam_rainfall_file <- function(path) {
  
  if (length(path) != 1L) {
    
    stop(
      "read_adam_rainfall_file() requires exactly one file.",
      call. = FALSE
    )
  }
  
  if (!file.exists(path)) {
    
    stop(
      "Rainfall source file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  filename <- basename(path)
  
  if (identical(
    filename,
    "adam_prcp_days_pre1900.csv"
  )) {
    
    stop(
      "The pre-1900 source file must be re-exported with ",
      "a four-digit ISO date before processing.",
      call. = FALSE
    )
  }
  
  period <- rainfall_period_from_filename(
    path
  )
  
  message(
    "Reading: ",
    filename
  )
  
  x <- data.table::fread(
    path,
    na.strings = c(
      "",
      "NA",
      "NULL",
      "(null)"
    ),
    showProgress = TRUE,
    colClasses = list(
      character = c(
        "STN_NUM",
        "LSD",
        "PRCP_TYPE"
      )
    )
  )
  
  data.table::setnames(
    x,
    tolower(
      names(x)
    )
  )
  
  required <- c(
    "stn_num",
    "lsd",
    "ob_qual_flag",
    "prcp",
    "prcp_per",
    "prcp_type",
    "prcp_qual",
    "prcp_day_cnt",
    "aws_flag"
  )
  
  missing_columns <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing_columns) > 0L) {
    
    stop(
      filename,
      " is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  station_integer <- suppressWarnings(
    as.integer(
      x$stn_num
    )
  )
  
  if (any(is.na(station_integer))) {
    
    stop(
      filename,
      " contains invalid station numbers.",
      call. = FALSE
    )
  }
  
  rainfall_value <- suppressWarnings(
    as.numeric(
      x$prcp
    )
  )
  
  nonblank_prcp <-
    !is.na(x$prcp) &
    trimws(
      as.character(x$prcp)
    ) != ""
  
  invalid_rainfall <-
    nonblank_prcp &
    is.na(rainfall_value)
  
  if (any(invalid_rainfall)) {
    
    examples <- unique(
      x$prcp[
        invalid_rainfall
      ]
    )
    
    stop(
      filename,
      " contains non-numeric rainfall values. Examples: ",
      paste(
        head(examples, 10L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  result <- data.table::data.table(
    
    station_number = sprintf(
      "%06d",
      station_integer
    ),
    
    date = parse_adam_rainfall_dates(
      date_text = x$lsd,
      start_year = period$start_year,
      end_year = period$end_year
    ),
    
    rainfall_mm = rainfall_value,
    
    observation_quality_flag =
      suppressWarnings(
        as.integer(
          x$ob_qual_flag
        )
      ),
    
    rainfall_period =
      suppressWarnings(
        as.numeric(
          x$prcp_per
        )
      ),
    
    rainfall_type =
      trimws(
        as.character(
          x$prcp_type
        )
      ),
    
    rainfall_quality =
      suppressWarnings(
        as.integer(
          x$prcp_qual
        )
      ),
    
    rainfall_day_count =
      suppressWarnings(
        as.numeric(
          x$prcp_day_cnt
        )
      ),
    
    aws_flag =
      suppressWarnings(
        as.integer(
          x$aws_flag
        )
      ),
    
    source_dataset =
      "adam_prcp_days",
    
    source_file =
      filename
  )
  
  result[
    rainfall_type == "",
    rainfall_type := NA_character_
  ]
  
  negative_count <- result[
    rainfall_mm < 0 &
      !is.na(rainfall_mm),
    .N
  ]
  
  if (negative_count > 0L) {
    
    warning(
      negative_count,
      " negative rainfall values found in ",
      filename,
      ". Values were retained for QC review.",
      call. = FALSE
    )
  }
  
  data.table::setorder(
    result,
    station_number,
    date
  )
  
  result
}