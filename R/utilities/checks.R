assert_columns <- function(data, required, data_name = deparse(substitute(data))) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(data_name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(data)
}
