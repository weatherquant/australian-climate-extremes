read_barra_inventory <- function(path) {
  if (!dir.exists(path)) stop("BARRA2 directory does not exist: ", path, call. = FALSE)
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  data.frame(path = files, extension = tools::file_ext(files), stringsAsFactors = FALSE)
}
