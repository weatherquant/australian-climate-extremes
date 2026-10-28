read_vector_layer <- function(path) sf::st_read(path, quiet = TRUE)
read_raster_layer <- function(path) terra::rast(path)
