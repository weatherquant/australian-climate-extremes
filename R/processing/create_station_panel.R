create_station_panel <- function(data, metadata) {
  assert_columns(data, c("station_id"))
  assert_columns(metadata, c("station_id"))
  dplyr::left_join(data, metadata, by = "station_id")
}
