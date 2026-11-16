# ============================================================
# build_rainfall_voronoi_network.R
#
# Build clipped Voronoi weights for the frozen Australian
# continental rainfall network.
#
# Purpose:
#   Provide station-density-adjusted spatial weights for the
#   rainfall concurrence / spatial-footprint analysis.
#
# Required inputs:
#   1. frozen station-network table containing at least:
#        station_number, station_name, latitude, longitude
#   2. mainland/continental Australian analysis-domain polygon
#
# IMPORTANT:
#   Set network_file and domain_file below to the existing project
#   products in this repository. The script deliberately does not
#   download or invent an Australia boundary.
#
# Outputs:
#   data-derived/rainfall/spatial_extremes/
#      rainfall_station_voronoi_weights.parquet
#      rainfall_station_voronoi.gpkg
#      rainfall_spatial_analysis_domain.gpkg
#
#   tables/rainfall_spatial_extremes/
#      rainfall_station_voronoi_audit.csv
#      rainfall_station_voronoi_summary.csv
#
# CRS:
#   EPSG:3577, GDA94 / Australian Albers, is used for area
#   calculations. Areas are converted from m^2 to km^2.
# ============================================================

# ============================================================
# Configuration: EDIT ONLY THESE INPUT PATHS IF REQUIRED
# ============================================================

network_file <- file.path(
  "data-derived",
  "rainfall",
  "network",
  "rainfall_continental_core_network_1950_2025.parquet"
)

domain_file <- file.path(
  "data-raw",
  "spatial",
  "abs_states_2021",
  "STE_2021_AUST_GDA2020.shp"
)

# If the domain GeoPackage has multiple layers, set this to the
# required layer name. Leave NULL for a single-layer file.
domain_layer <- NULL

output_dir <- file.path(
  "data-derived",
  "rainfall",
  "spatial_extremes"
)

table_output_dir <- file.path(
  "tables",
  "rainfall_spatial_extremes"
)

output_crs <- 3577L

expected_station_count <- 742L

# Small numerical tolerance for checking that clipped Voronoi areas
# sum to the analysis-domain area.
area_relative_tolerance <- 1e-6


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "data.table",
  "arrow",
  "sf"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(data.table)
library(arrow)
library(sf)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Input checks
# ============================================================

if (!file.exists(network_file)) {
  stop(
    "Frozen station-network file not found: ",
    network_file,
    "\nUpdate network_file at the top of the script to the existing 742-station network product.",
    call. = FALSE
  )
}

if (!file.exists(domain_file)) {
  stop(
    "Australian analysis-domain polygon not found: ",
    domain_file,
    "\nUpdate domain_file at the top of the script to the existing mainland/continental analysis-domain polygon.",
    call. = FALSE
  )
}


# ============================================================
# Read frozen network
# ============================================================

network <- as.data.table(
  arrow::read_parquet(
    network_file
  )
)

required_network_columns <- c(
  "station_number",
  "station_name",
  "latitude",
  "longitude"
)

missing_network_columns <- setdiff(
  required_network_columns,
  names(network)
)

if (length(missing_network_columns) > 0L) {
  stop(
    "Frozen network is missing columns: ",
    paste(missing_network_columns, collapse = ", "),
    call. = FALSE
  )
}

# If the supplied product contains more than the already-frozen core
# network and carries the established continental_core flag, use it.
if (
  "continental_core" %in% names(network) &&
  nrow(network) != expected_station_count
) {
  network <- network[continental_core == TRUE]
}

if (nrow(network) != expected_station_count) {
  stop(
    "Expected the frozen network to contain ",
    expected_station_count,
    " stations but found ",
    nrow(network),
    ". Do not proceed until network_file points to the established core network.",
    call. = FALSE
  )
}

if (anyNA(network$station_number)) {
  stop("Missing station_number values in frozen network.", call. = FALSE)
}

if (anyDuplicated(network$station_number)) {
  stop("Duplicated station_number values in frozen network.", call. = FALSE)
}

if (
  any(!is.finite(network$latitude)) ||
  any(!is.finite(network$longitude))
) {
  stop("Non-finite station coordinates found in frozen network.", call. = FALSE)
}

# Duplicate coordinates are fatal for an unambiguous station Voronoi
# tessellation. Stop rather than silently jittering stations.
coordinate_duplicates <- network[
  ,
  .N,
  by = .(longitude, latitude)
][N > 1L]

if (nrow(coordinate_duplicates) > 0L) {
  stop(
    "The frozen network contains duplicated station coordinates. ",
    "Inspect these before constructing Voronoi weights.",
    call. = FALSE
  )
}

setorder(network, station_number)


# ============================================================
# Read and prepare analysis domain
# ============================================================

if (is.null(domain_layer)) {
  domain <- sf::st_read(
    domain_file,
    quiet = TRUE
  )
} else {
  domain <- sf::st_read(
    domain_file,
    layer = domain_layer,
    quiet = TRUE
  )
}

# ============================================================
# Restrict domain to the eight standard Australian
# states and territories
# ============================================================

analysis_states <- c(
  "New South Wales",
  "Victoria",
  "Queensland",
  "South Australia",
  "Western Australia",
  "Tasmania",
  "Northern Territory",
  "Australian Capital Territory"
)

domain <- domain[
  domain$STE_NAME21 %in%
    analysis_states,
]

if (nrow(domain) != 8L) {
  stop(
    "Expected 8 Australian state/territory features after filtering, found ",
    nrow(domain),
    ".",
    call. = FALSE
  )
}

if (
  any(
    domain$STE_NAME21 %in%
    c(
      "Other Territories",
      "Outside Australia"
    )
  )
) {
  stop(
    "External or outside-Australia features remain in the analysis domain.",
    call. = FALSE
  )
}

message(
  "Analysis domain contains: ",
  paste(
    domain$STE_NAME21,
    collapse = ", "
  )
)

if (nrow(domain) == 0L) {
  stop("Analysis-domain polygon contains no features.", call. = FALSE)
}

if (is.na(sf::st_crs(domain))) {
  stop(
    "Analysis-domain polygon has no CRS. Assign its known CRS before running this script.",
    call. = FALSE
  )
}

# Make valid before unioning. st_make_valid is available via sf/lwgeom
# depending on installation, but sf exposes it directly.
domain <- sf::st_make_valid(domain)

domain <- sf::st_transform(
  domain,
  output_crs
)

# Dissolve all supplied mainland-analysis-domain pieces into one
# MULTIPOLYGON/POLYGON boundary. This preserves islands only if they
# are deliberately present in the supplied project domain.
domain_union <- sf::st_sf(
  domain_id = "continental_analysis_domain",
  geometry = sf::st_union(sf::st_geometry(domain))
)

sf::st_crs(domain_union) <- sf::st_crs(domain)

domain_union <- sf::st_make_valid(domain_union)

domain_area_km2 <- as.numeric(
  sf::st_area(domain_union)
) / 1e6

if (!is.finite(domain_area_km2) || domain_area_km2 <= 0) {
  stop("Invalid analysis-domain area produced.", call. = FALSE)
}


# ============================================================
# Convert stations to sf and project
# ============================================================

stations_sf <- sf::st_as_sf(
  network,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

stations_sf <- sf::st_transform(
  stations_sf,
  output_crs
)

inside_domain <- lengths(
  sf::st_intersects(
    stations_sf,
    domain_union
  )
) > 0L

if (!all(inside_domain)) {
  outside <- network$station_number[!inside_domain]

  stop(
    length(outside),
    " frozen-network stations fall outside the supplied analysis domain: ",
    paste(head(outside, 20L), collapse = ", "),
    if (length(outside) > 20L) " ..." else "",
    call. = FALSE
  )
}


# ============================================================
# Construct Voronoi cells
# ============================================================

# st_voronoi returns polygons in the order of the input point geometry
# for ordinary non-duplicated points. We nevertheless allocate station
# IDs via point-in-polygon intersection below instead of relying solely
# on ordering.
station_union <- sf::st_union(
  sf::st_geometry(stations_sf)
)

# Use the domain envelope as a finite tessellation envelope and clip
# the resulting cells to the actual analysis polygon immediately.
voronoi_raw <- sf::st_voronoi(
  station_union,
  envelope = sf::st_as_sfc(sf::st_bbox(domain_union))
)

voronoi_polygons <- sf::st_collection_extract(
  voronoi_raw,
  "POLYGON"
)

voronoi_sf <- sf::st_sf(
  voronoi_temp_id = seq_along(voronoi_polygons),
  geometry = voronoi_polygons
)

sf::st_crs(voronoi_sf) <- sf::st_crs(domain_union)

if (nrow(voronoi_sf) != nrow(stations_sf)) {
  stop(
    "Voronoi polygon count differs from station count: ",
    nrow(voronoi_sf),
    " polygons versus ",
    nrow(stations_sf),
    " stations.",
    call. = FALSE
  )
}


# ============================================================
# Attach each cell to its generating station
# ============================================================

# st_intersects handles points on polygon boundaries more robustly than
# st_within. A valid Voronoi cell should contain exactly one generator.
cell_station_hits <- sf::st_intersects(
  voronoi_sf,
  stations_sf
)

hit_counts <- lengths(cell_station_hits)

if (any(hit_counts != 1L)) {
  bad_cells <- which(hit_counts != 1L)

  stop(
    "Could not unambiguously associate ",
    length(bad_cells),
    " Voronoi cells with exactly one station. Bad temporary cell IDs: ",
    paste(head(bad_cells, 20L), collapse = ", "),
    call. = FALSE
  )
}

station_index <- vapply(
  cell_station_hits,
  `[[`,
  integer(1),
  1L
)

voronoi_sf$station_number <- stations_sf$station_number[station_index]

if (anyDuplicated(voronoi_sf$station_number)) {
  stop(
    "Voronoi-to-station assignment produced duplicated station IDs.",
    call. = FALSE
  )
}


# ============================================================
# Clip Voronoi cells to the analysis domain
# ============================================================

voronoi_clipped <- suppressWarnings(
  sf::st_intersection(
    voronoi_sf,
    domain_union
  )
)

voronoi_clipped <- sf::st_collection_extract(
  voronoi_clipped,
  "POLYGON",
  warn = FALSE
)

if (nrow(voronoi_clipped) == 0L) {
  stop("Voronoi/domain intersection produced zero polygons.", call. = FALSE)
}

# A clipped cell can be multipart after intersection. Dissolve all pieces
# back to one feature per station before calculating station weights.
voronoi_clipped <- aggregate(
  voronoi_clipped["geometry"],
  by = list(
    station_number = voronoi_clipped$station_number
  ),
  FUN = function(x) x[1]
)

# aggregate.sf unions geometries belonging to each group.
voronoi_clipped <- sf::st_make_valid(voronoi_clipped)

if (nrow(voronoi_clipped) != expected_station_count) {
  stop(
    "Expected ",
    expected_station_count,
    " clipped station cells but found ",
    nrow(voronoi_clipped),
    ".",
    call. = FALSE
  )
}


# ============================================================
# Reattach station metadata and compute weights
# ============================================================

station_metadata <- network[, .(
  station_number,
  station_name,
  latitude,
  longitude,
  state = if ("state" %in% names(network)) state else NA_character_,
  rainfall_district = if ("rainfall_district" %in% names(network)) rainfall_district else NA_character_,
  area_name = if ("area_name" %in% names(network)) area_name else NA_character_
)]

voronoi_clipped <- merge(
  voronoi_clipped,
  station_metadata,
  by = "station_number",
  all.x = TRUE,
  sort = FALSE
)

voronoi_clipped$voronoi_area_km2 <- as.numeric(
  sf::st_area(voronoi_clipped)
) / 1e6

if (
  any(!is.finite(voronoi_clipped$voronoi_area_km2)) ||
  any(voronoi_clipped$voronoi_area_km2 <= 0)
) {
  stop("At least one clipped Voronoi cell has invalid area.", call. = FALSE)
}

voronoi_clipped$continental_area_weight <-
  voronoi_clipped$voronoi_area_km2 /
  sum(voronoi_clipped$voronoi_area_km2)

voronoi_clipped$domain_fraction <-
  voronoi_clipped$voronoi_area_km2 /
  domain_area_km2


# ============================================================
# Validate tessellation area
# ============================================================

voronoi_total_area_km2 <- sum(
  voronoi_clipped$voronoi_area_km2
)

area_difference_km2 <-
  voronoi_total_area_km2 -
  domain_area_km2

area_relative_difference <-
  area_difference_km2 /
  domain_area_km2

if (
  abs(area_relative_difference) >
  area_relative_tolerance
) {
  stop(
    "Clipped Voronoi areas do not reconstruct the analysis-domain area within tolerance. ",
    "Relative difference = ",
    signif(area_relative_difference, 6),
    call. = FALSE
  )
}

if (
  abs(sum(voronoi_clipped$continental_area_weight) - 1) >
  area_relative_tolerance
) {
  stop(
    "Voronoi continental-area weights do not sum to one.",
    call. = FALSE
  )
}


# ============================================================
# Neighbour diagnostics
# ============================================================

# Queen-style polygon touching/intersection. This is diagnostic only;
# later 100/200/300-km rainfall-cluster connectivity will be calculated
# independently from station coordinates.
neighbours <- sf::st_touches(voronoi_clipped)

voronoi_clipped$voronoi_neighbours <- lengths(neighbours)


# ============================================================
# Non-spatial weight table
# ============================================================

weights_table <- sf::st_drop_geometry(voronoi_clipped)
weights_table <- as.data.table(weights_table)

setcolorder(
  weights_table,
  c(
    "station_number",
    "station_name",
    "latitude",
    "longitude",
    "state",
    "rainfall_district",
    "area_name",
    "voronoi_area_km2",
    "continental_area_weight",
    "domain_fraction",
    "voronoi_neighbours"
  )
)

setorder(weights_table, station_number)


# ============================================================
# Audit tables
# ============================================================

voronoi_audit <- data.table(
  stations = nrow(weights_table),
  crs_epsg = output_crs,
  domain_area_km2 = domain_area_km2,
  voronoi_total_area_km2 = voronoi_total_area_km2,
  area_difference_km2 = area_difference_km2,
  area_relative_difference = area_relative_difference,
  weight_sum = sum(weights_table$continental_area_weight),
  minimum_cell_area_km2 = min(weights_table$voronoi_area_km2),
  median_cell_area_km2 = median(weights_table$voronoi_area_km2),
  maximum_cell_area_km2 = max(weights_table$voronoi_area_km2),
  minimum_neighbours = min(weights_table$voronoi_neighbours),
  median_neighbours = median(weights_table$voronoi_neighbours),
  maximum_neighbours = max(weights_table$voronoi_neighbours)
)

voronoi_summary <- weights_table[, .(
  stations = .N,
  represented_area_km2 = sum(voronoi_area_km2),
  median_cell_area_km2 = median(voronoi_area_km2),
  q10_cell_area_km2 = quantile(voronoi_area_km2, 0.10, names = FALSE),
  q90_cell_area_km2 = quantile(voronoi_area_km2, 0.90, names = FALSE),
  minimum_cell_area_km2 = min(voronoi_area_km2),
  maximum_cell_area_km2 = max(voronoi_area_km2)
), by = state]

setorder(voronoi_summary, state)


# ============================================================
# Save outputs
# ============================================================

arrow::write_parquet(
  weights_table,
  file.path(
    output_dir,
    "rainfall_station_voronoi_weights.parquet"
  ),
  compression = "zstd"
)

voronoi_file <- file.path(
  output_dir,
  "rainfall_station_voronoi.gpkg"
)

if (file.exists(voronoi_file)) {
  file.remove(voronoi_file)
}

sf::st_write(
  voronoi_clipped,
  voronoi_file,
  layer = "station_voronoi",
  quiet = TRUE
)

domain_output_file <- file.path(
  output_dir,
  "rainfall_spatial_analysis_domain.gpkg"
)

if (file.exists(domain_output_file)) {
  file.remove(domain_output_file)
}

sf::st_write(
  domain_union,
  domain_output_file,
  layer = "analysis_domain",
  quiet = TRUE
)

fwrite(
  voronoi_audit,
  file.path(
    table_output_dir,
    "rainfall_station_voronoi_audit.csv"
  )
)

fwrite(
  voronoi_summary,
  file.path(
    table_output_dir,
    "rainfall_station_voronoi_summary.csv"
  )
)


# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nRAINFALL VORONOI NETWORK COMPLETE",
  "\n========================================"
)

print(voronoi_audit)
print(voronoi_summary)

message(
  "\nInspect the largest cells before proceeding to daily concurrence:"
)

print(
  weights_table[
    order(-voronoi_area_km2),
    .(
      station_number,
      station_name,
      state,
      latitude,
      longitude,
      voronoi_area_km2,
      continental_area_weight,
      voronoi_neighbours
    )
  ][1:min(20L, .N)]
)
