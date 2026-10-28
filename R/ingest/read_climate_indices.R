read_climate_index <- function(path, index_name, source) {
  stopifnot(file.exists(path))
  x <- data.table::fread(path)
  attr(x, "index_name") <- index_name
  attr(x, "source") <- source
  x
}
