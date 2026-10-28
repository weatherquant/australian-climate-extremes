standardise_station_daily <- function(data) {
  required <- c("station_id", "date", "variable", "value_original", "unit",
                "source_dataset", "data_version", "source_file")
  assert_columns(data, required)
  data |>
    dplyr::mutate(
      station_id = as.character(.data$station_id),
      date = as.Date(.data$date),
      variable = as.character(.data$variable),
      value_clean = as.numeric(.data$value_original),
      quality_status = dplyr::coalesce(as.character(.data$quality_status), "unknown")
    )
}
