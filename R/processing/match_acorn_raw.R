match_raw_homogenised_temperature <- function(raw, homogenised) {
  keys <- c("station_id", "date", "variable")
  assert_columns(raw, c(keys, "value_clean"), "raw")
  assert_columns(homogenised, c(keys, "value_clean"), "homogenised")
  dplyr::full_join(
    raw |> dplyr::select(dplyr::all_of(keys), raw_value = .data$value_clean),
    homogenised |> dplyr::select(dplyr::all_of(keys), homogenised_value = .data$value_clean),
    by = keys
  ) |>
    dplyr::mutate(adjustment = .data$homogenised_value - .data$raw_value)
}
