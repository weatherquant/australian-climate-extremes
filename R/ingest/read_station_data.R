read_station_file <- function(path, source_dataset, data_version) {
  stopifnot(file.exists(path))
  raw <- data.table::fread(path, na.strings = c("", "NA", "NaN", "-999", "-9999"))
  attr(raw, "source_file") <- basename(path)
  attr(raw, "source_dataset") <- source_dataset
  attr(raw, "data_version") <- data_version
  raw
}

# Source-specific parsers should call read_station_file() and return the
# canonical station_daily schema. Add parsers once actual Bureau file layouts
# are known; do not guess source column positions here.
