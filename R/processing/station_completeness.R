calculate_station_availability <- function(data) {
  assert_columns(data, c("station_id", "date", "variable", "value_clean", "data_version"))
  data |>
    dplyr::mutate(year = as.integer(format(.data$date, "%Y"))) |>
    dplyr::group_by(.data$station_id, .data$variable, .data$data_version, .data$year) |>
    dplyr::summarise(
      n_observed = dplyr::n(),
      n_valid = sum(!is.na(.data$value_clean)),
      n_missing = sum(is.na(.data$value_clean)),
      .groups = "drop"
    )
}
