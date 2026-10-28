read_station_metadata <- function(path) {
  stopifnot(file.exists(path))
  data.table::fread(path)
}
