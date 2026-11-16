# ============================================================
# build_rainfall_spatial_extreme_clusters.R
#
# Build connected daily spatial clusters of P97.5 extreme rainfall
# for the frozen 742-station Australian rainfall network.
#
# Primary definition:
#   Threshold:                 station-specific P97.5
#   Connectivity distance:     200 km
#   Small connected cluster:   >= 2 stations
#   Coherent spatial cluster:  >= 3 stations
#   Primary coverage:          >=80% stations and >=80% represented area
#   Strict sensitivity:        >=90% stations and >=90% represented area
#
# Connectivity is based only on physical station-to-station distance.
# Voronoi area describes represented footprint but does not create links.
# Single remote stations therefore remain isolated extremes.
# ============================================================

# ============================================================
# Configuration
# ============================================================

station_day_chunk_dir <- file.path("data-derived",
                                   "rainfall",
                                   "spatial_extremes",
                                   "station_day_chunks")

daily_summary_file <- file.path(
  "data-derived",
  "rainfall",
  "spatial_extremes",
  "rainfall_spatial_extreme_days_1950_2025.parquet"
)

network_file <- file.path(
  "data-derived",
  "rainfall",
  "network",
  "rainfall_continental_core_network_1950_2025.parquet"
)

voronoi_weight_file <- file.path(
  "data-derived",
  "rainfall",
  "spatial_extremes",
  "rainfall_station_voronoi_weights.parquet"
)

output_dir <- file.path("data-derived", "rainfall", "spatial_extremes")

table_output_dir <- file.path("tables", "rainfall_spatial_extremes")

cluster_output_file <- file.path(output_dir,
                                 "rainfall_spatial_extreme_clusters_p975_200km_1950_2025.parquet")

daily_cluster_output_file <- file.path(
  output_dir,
  "rainfall_spatial_extreme_cluster_days_p975_200km_1950_2025.parquet"
)

annual_cluster_output_file <- file.path(
  output_dir,
  "rainfall_spatial_extreme_cluster_annual_p975_200km_1950_2025.parquet"
)

threshold_id_primary <- "p975"
connectivity_distance_km <- 200
minimum_connected_stations <- 2L
minimum_coherent_stations <- 3L
output_crs <- 3577L
expected_station_count <- 742L

# ============================================================
# Packages
# ============================================================

required_packages <- c("data.table", "arrow", "sf", "igraph")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

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
library(igraph)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# ============================================================
# Input checks
# ============================================================

required_files <- c(daily_summary_file, network_file, voronoi_weight_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Missing required files:\n",
       paste(missing_files, collapse = "\n"),
       call. = FALSE)
}

if (!dir.exists(station_day_chunk_dir)) {
  stop("Station-day chunk directory not found: ",
       station_day_chunk_dir,
       call. = FALSE)
}

chunk_files <- sort(
  list.files(
    station_day_chunk_dir,
    pattern = "_spatial_extremes\\.parquet$",
    full.names = TRUE
  )
)

if (length(chunk_files) != 16L) {
  stop("Expected 16 station-day chunks, found ",
       length(chunk_files),
       ".",
       call. = FALSE)
}

# ============================================================
# Static station geometry and weights
# ============================================================

network <- as.data.table(read_parquet(network_file))
weights <- as.data.table(read_parquet(voronoi_weight_file))

if (nrow(network) != expected_station_count) {
  stop("Expected 742 stations in frozen network.", call. = FALSE)
}

station_static <- merge(network[, .(station_number = as.character(station_number),
                                    station_name,
                                    state,
                                    latitude,
                                    longitude)],
                        weights[, .(station_number = as.character(station_number),
                                    voronoi_area_km2,
                                    continental_area_weight)],
                        by = "station_number",
                        all = TRUE,
                        sort = FALSE)

if (nrow(station_static) != expected_station_count ||
    anyNA(station_static$voronoi_area_km2)) {
  stop("Network-to-Voronoi join failed.", call. = FALSE)
}

station_sf <- st_as_sf(
  station_static,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

station_sf <- st_transform(station_sf, output_crs)
station_xy <- st_coordinates(station_sf)

station_static[, `:=`(x_m = station_xy[, 1], y_m = station_xy[, 2])]

setkey(station_static, station_number)

# Precompute the fixed 742-station 200-km adjacency matrix. Daily graphs
# are induced subgraphs of this fixed physical-distance network.
xy_matrix <- as.matrix(station_static[, .(x_m, y_m)])
distance_matrix_km <- as.matrix(dist(xy_matrix)) / 1000
adjacency_matrix <- distance_matrix_km <= connectivity_distance_km
diag(adjacency_matrix) <- FALSE

# ============================================================
# Daily coverage metadata
# ============================================================

daily_coverage <- as.data.table(read_parquet(
  daily_summary_file,
  col_select = c(
    "date",
    "year",
    "month",
    "season",
    "threshold_id",
    "eligible_stations",
    "eligible_station_fraction",
    "eligible_area_km2",
    "eligible_area_fraction",
    "coverage_primary",
    "coverage_strict",
    "coverage_tier"
  )
))

daily_coverage <- daily_coverage[threshold_id == threshold_id_primary]
setkey(daily_coverage, date)

# ============================================================
# Cluster helper
# ============================================================

cluster_one_day <- function(day_events) {
  n <- nrow(day_events)
  if (n == 0L)
    return(NULL)
  
  station_index <- match(day_events$station_number, station_static$station_number)
  
  if (anyNA(station_index)) {
    stop("An exceeding station is absent from station_static.", call. = FALSE)
  }
  
  if (n == 1L) {
    membership <- 1L
  } else {
    day_adjacency <- adjacency_matrix[station_index, station_index, drop = FALSE]
    day_graph <- graph_from_adjacency_matrix(day_adjacency, mode = "undirected", diag = FALSE)
    membership <- components(day_graph)$membership
  }
  
  day_events[, component_number := as.integer(membership)]
  
  cluster_table <- day_events[, {
    station_count_value <- .N
    area_value <- sum(voronoi_area_km2, na.rm = TRUE)
    weight_value <- sum(continental_area_weight, na.rm = TRUE)
    
    centroid_x <- weighted.mean(x_m,
                                w = pmax(voronoi_area_km2, 1e-9),
                                na.rm = TRUE)
    
    centroid_y <- weighted.mean(y_m,
                                w = pmax(voronoi_area_km2, 1e-9),
                                na.rm = TRUE)
    
    list(
      station_count = as.integer(station_count_value),
      connected_cluster = station_count_value >= minimum_connected_stations,
      coherent_cluster = station_count_value >= minimum_coherent_stations,
      represented_area_km2 = as.numeric(area_value),
      represented_area_fraction = as.numeric(weight_value),
      states_count = as.integer(uniqueN(state)),
      states = paste(sort(unique(state)), collapse = ";"),
      station_numbers = paste(sort(station_number), collapse = ";"),
      station_names = paste(station_name[order(station_number)], collapse = ";"),
      centroid_x_m = as.numeric(centroid_x),
      centroid_y_m = as.numeric(centroid_y),
      maximum_rainfall_mm = max(rainfall_mm, na.rm = TRUE),
      median_rainfall_mm = median(rainfall_mm, na.rm = TRUE),
      maximum_threshold_ratio = max(threshold_ratio, na.rm = TRUE),
      median_threshold_ratio = median(threshold_ratio, na.rm = TRUE),
      maximum_excess_mm = max(excess_mm, na.rm = TRUE),
      median_excess_mm = median(excess_mm, na.rm = TRUE)
    )
  }, by = component_number]
  
  centroid_sf <- st_as_sf(
    cluster_table,
    coords = c("centroid_x_m", "centroid_y_m"),
    crs = output_crs,
    remove = FALSE
  )
  
  centroid_sf <- st_transform(centroid_sf, 4326)
  centroid_lonlat <- st_coordinates(centroid_sf)
  
  cluster_table[, `:=`(centroid_longitude = centroid_lonlat[, 1],
                       centroid_latitude = centroid_lonlat[, 2])]
  
  setorder(cluster_table, -station_count, -represented_area_km2)
  cluster_table[, cluster_rank_station_count := seq_len(.N)]
  
  area_order <- order(-cluster_table$represented_area_km2,-cluster_table$station_count)
  cluster_table[, cluster_rank_area := NA_integer_]
  cluster_table$cluster_rank_area[area_order] <- seq_along(area_order)
  
  cluster_table[, cluster_id := sprintf(
    "%s_%s_C%03d",
    format(day_events$date[1], "%Y%m%d"),
    threshold_id_primary,
    cluster_rank_station_count
  )]
  
  cluster_table[, `:=`(
    date = day_events$date[1],
    year = day_events$year[1],
    month = day_events$month[1],
    season = day_events$season[1],
    threshold_id = threshold_id_primary,
    connectivity_distance_km = connectivity_distance_km
  )]
  
  setcolorder(
    cluster_table,
    c(
      "cluster_id",
      "date",
      "year",
      "month",
      "season",
      "threshold_id",
      "connectivity_distance_km",
      "component_number",
      "cluster_rank_station_count",
      "cluster_rank_area",
      "station_count",
      "connected_cluster",
      "coherent_cluster",
      "represented_area_km2",
      "represented_area_fraction",
      "states_count",
      "states",
      "centroid_latitude",
      "centroid_longitude",
      "maximum_rainfall_mm",
      "median_rainfall_mm",
      "maximum_threshold_ratio",
      "median_threshold_ratio",
      "maximum_excess_mm",
      "median_excess_mm",
      "station_numbers",
      "station_names"
    )
  )
  
  cluster_table
}

# ============================================================
# Process station-day chunks incrementally
# ============================================================

cluster_chunk_dir <- file.path(output_dir, "cluster_chunks_p975_200km")
dir.create(cluster_chunk_dir,
           recursive = TRUE,
           showWarnings = FALSE)

chunk_audit_list <- vector("list", length(chunk_files))
cluster_summary_list <- vector("list", length(chunk_files))

for (i in seq_along(chunk_files)) {
  file_path <- chunk_files[i]
  
  message("Cluster file ",
          i,
          " of ",
          length(chunk_files),
          ": ",
          basename(file_path))
  
  station_days <- as.data.table(read_parquet(
    file_path,
    col_select = c(
      "station_number",
      "station_name",
      "state",
      "latitude",
      "longitude",
      "date",
      "year",
      "month",
      "season",
      "rainfall_mm",
      "threshold_id",
      "eligible",
      "exceeds",
      "threshold_ratio",
      "excess_mm",
      "voronoi_area_km2",
      "continental_area_weight"
    )
  ))
  
  extremes <- station_days[threshold_id == threshold_id_primary &
                             eligible == TRUE &
                             exceeds == TRUE]
  
  extremes <- station_static[extremes, on = "station_number", nomatch = 0L, .(
    station_number = i.station_number,
    station_name = i.station_name,
    state = i.state,
    latitude = i.latitude,
    longitude = i.longitude,
    x_m,
    y_m,
    date = i.date,
    year = i.year,
    month = i.month,
    season = i.season,
    rainfall_mm = i.rainfall_mm,
    threshold_ratio = i.threshold_ratio,
    excess_mm = i.excess_mm,
    voronoi_area_km2 = i.voronoi_area_km2,
    continental_area_weight = i.continental_area_weight
  )]
  
  if (nrow(extremes) == 0L) {
    clusters_chunk <- data.table()
  } else {
    clusters_chunk <- rbindlist(lapply(split(
      extremes, by = "date", keep.by = TRUE
    ), cluster_one_day),
    use.names = TRUE,
    fill = TRUE)
  }
  
  cluster_chunk_file <- file.path(
    cluster_chunk_dir,
    sub(
      "_spatial_extremes\\.parquet$",
      "_clusters_p975_200km.parquet",
      basename(file_path)
    )
  )
  
  if (nrow(clusters_chunk) > 0L) {
    write_parquet(clusters_chunk, cluster_chunk_file, compression = "zstd")
  }
  
  cluster_summary_list[[i]] <- clusters_chunk
  
  chunk_audit_list[[i]] <- data.table(
    file = basename(file_path),
    exceeding_station_days = as.integer(nrow(extremes)),
    extreme_dates = as.integer(uniqueN(extremes$date)),
    clusters = as.integer(nrow(clusters_chunk)),
    connected_clusters = as.integer(sum(
      clusters_chunk$connected_cluster %in% TRUE
    )),
    coherent_clusters = as.integer(sum(
      clusters_chunk$coherent_cluster %in% TRUE
    )),
    isolated_clusters = as.integer(sum(clusters_chunk$station_count == 1L))
  )
  
  rm(station_days, extremes, clusters_chunk)
  invisible(gc(full = TRUE))
}

# The cluster table is much smaller than the station-day source and can
# safely be combined after incremental processing.
clusters <- rbindlist(cluster_summary_list, use.names = TRUE, fill = TRUE)
chunk_audit <- rbindlist(chunk_audit_list, use.names = TRUE, fill = TRUE)

rm(cluster_summary_list, chunk_audit_list)
invisible(gc(full = TRUE))

setorder(clusters, date, cluster_rank_station_count)

# ============================================================
# Attach daily coverage metadata
# ============================================================

clusters <- daily_coverage[clusters, on = "date", nomatch = 0L]

# ============================================================
# Daily cluster summaries
#
# Keep the largest-by-station and largest-by-area components
# separate. Do not combine attributes from different clusters.
# ============================================================

daily_clusters <- clusters[, {
  count_index <- which.max(station_count)
  
  area_index <- which.max(represented_area_km2)
  
  list(
    total_clusters =
      as.integer(.N),
    
    isolated_extremes =
      as.integer(sum(station_count == 1L)),
    
    connected_clusters =
      as.integer(sum(connected_cluster)),
    
    coherent_clusters =
      as.integer(sum(coherent_cluster)),
    
    stations_in_all_clusters =
      as.integer(sum(station_count)),
    
    stations_in_connected_clusters =
      as.integer(sum(station_count[connected_cluster])),
    
    stations_in_coherent_clusters =
      as.integer(sum(station_count[coherent_cluster])),
    
    # ------------------------------------------
    # Largest component by station count
    # ------------------------------------------
    
    largest_count_cluster_id =
      cluster_id[count_index],
    
    largest_count_cluster_stations =
      as.integer(station_count[count_index]),
    
    largest_count_cluster_area_km2 =
      as.numeric(represented_area_km2[count_index]),
    
    largest_count_cluster_area_fraction =
      as.numeric(represented_area_fraction[count_index]),
    
    largest_count_cluster_states =
      states[count_index],
    
    largest_count_cluster_maximum_rainfall_mm =
      as.numeric(maximum_rainfall_mm[count_index]),
    
    largest_count_cluster_maximum_threshold_ratio =
      as.numeric(maximum_threshold_ratio[count_index]),
    
    # ------------------------------------------
    # Largest component by represented area
    # ------------------------------------------
    
    largest_area_cluster_id =
      cluster_id[area_index],
    
    largest_area_cluster_stations =
      as.integer(station_count[area_index]),
    
    largest_area_cluster_area_km2 =
      as.numeric(represented_area_km2[area_index]),
    
    largest_area_cluster_area_fraction =
      as.numeric(represented_area_fraction[area_index]),
    
    largest_area_cluster_states =
      states[area_index],
    
    largest_area_cluster_maximum_rainfall_mm =
      as.numeric(maximum_rainfall_mm[area_index]),
    
    largest_area_cluster_maximum_threshold_ratio =
      as.numeric(maximum_threshold_ratio[area_index]),
    
    total_exceeding_represented_area_km2 =
      as.numeric(sum(represented_area_km2)),
    
    has_connected_cluster =
      any(connected_cluster),
    
    has_coherent_cluster =
      any(coherent_cluster)
  )
}, by = .(
  date,
  year,
  month,
  season,
  threshold_id,
  connectivity_distance_km,
  eligible_stations,
  eligible_station_fraction,
  eligible_area_km2,
  eligible_area_fraction,
  coverage_primary,
  coverage_strict,
  coverage_tier
)]

setorder(daily_clusters, date)

# ============================================================
# Annual cluster audit
# ============================================================

annual_clusters <- daily_clusters[, {
  primary_index <-
    coverage_primary == TRUE
  
  list(
    extreme_days =
      as.integer(.N),
    
    primary_coverage_extreme_days =
      as.integer(sum(primary_index)),
    
    strict_coverage_extreme_days =
      as.integer(sum(coverage_strict)),
    
    connected_cluster_days =
      as.integer(sum(has_connected_cluster &
                       primary_index)),
    
    coherent_cluster_days =
      as.integer(sum(has_coherent_cluster &
                       primary_index)),
    
    maximum_largest_count_cluster_stations =
      as.integer(max(
        largest_count_cluster_stations[primary_index], na.rm = TRUE
      )),
    
    maximum_largest_count_cluster_area_km2 =
      as.numeric(max(
        largest_count_cluster_area_km2[primary_index], na.rm = TRUE
      )),
    
    maximum_largest_area_cluster_area_km2 =
      as.numeric(max(
        largest_area_cluster_area_km2[primary_index], na.rm = TRUE
      )),
    
    median_largest_count_cluster_stations =
      as.numeric(median(
        largest_count_cluster_stations[primary_index], na.rm = TRUE
      )),
    
    p90_largest_count_cluster_stations =
      as.numeric(
        quantile(
          largest_count_cluster_stations[primary_index],
          0.90,
          na.rm = TRUE,
          names = FALSE
        )
      ),
    
    median_largest_count_cluster_area_km2 =
      as.numeric(median(
        largest_count_cluster_area_km2[primary_index], na.rm = TRUE
      )),
    
    p90_largest_count_cluster_area_km2 =
      as.numeric(
        quantile(
          largest_count_cluster_area_km2[primary_index],
          0.90,
          na.rm = TRUE,
          names = FALSE
        )
      )
  )
}, by = year]

setorder(annual_clusters, year)

# ============================================================
# Summary diagnostics
# ============================================================

cluster_size_summary <- clusters[, .(
  clusters =
    as.integer(.N),
  
  median_station_count =
    as.numeric(median(station_count)),
  
  p90_station_count =
    as.numeric(quantile(station_count, 0.90, names = FALSE)),
  
  p99_station_count =
    as.numeric(quantile(station_count, 0.99, names = FALSE)),
  
  median_area_km2 =
    as.numeric(median(represented_area_km2)),
  
  p90_area_km2 =
    as.numeric(quantile(represented_area_km2, 0.90, names = FALSE)),
  
  p99_area_km2 =
    as.numeric(quantile(represented_area_km2, 0.99, names = FALSE))
), by = .(
  cluster_class =
    data.table::fcase(
      station_count == 1L,
      "isolated_1",
      
      station_count == 2L,
      "connected_2",
      
      station_count >= 3L,
      "coherent_3plus",
      
      default = NA_character_
    )
)]

largest_coherent_days <- daily_clusters[coverage_primary == TRUE &
                                          has_coherent_cluster == TRUE][order(-largest_count_cluster_stations,-largest_count_cluster_area_km2)][1:min(50L, .N)]

# ============================================================
# Save outputs
# ============================================================

write_parquet(clusters, cluster_output_file, compression = "zstd")
write_parquet(daily_clusters, daily_cluster_output_file, compression = "zstd")
write_parquet(annual_clusters, annual_cluster_output_file, compression = "zstd")

fwrite(
  chunk_audit,
  file.path(
    table_output_dir,
    "rainfall_spatial_cluster_chunk_audit.csv"
  )
)

fwrite(
  cluster_size_summary,
  file.path(
    table_output_dir,
    "rainfall_spatial_cluster_size_summary.csv"
  )
)

fwrite(
  annual_clusters,
  file.path(
    table_output_dir,
    "rainfall_spatial_cluster_annual_summary.csv"
  )
)

fwrite(
  largest_coherent_days,
  file.path(
    table_output_dir,
    "rainfall_spatial_cluster_largest_days.csv"
  )
)

# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nDAILY RAINFALL SPATIAL CLUSTERS COMPLETE",
  "\n========================================"
)

print(chunk_audit)
print(cluster_size_summary)

message("\nCoverage and cluster-day summary:")
print(daily_clusters[, .(
  extreme_days = .N,
  primary_coverage_days = sum(coverage_primary),
  strict_coverage_days = sum(coverage_strict),
  connected_cluster_days = sum(has_connected_cluster),
  coherent_cluster_days = sum(has_coherent_cluster)
)])

message("\nLargest coherent P97.5 cluster days:")
print(largest_coherent_days[, .(
  date,
  coverage_tier,
  
  largest_count_cluster_id,
  largest_count_cluster_stations,
  largest_count_cluster_area_km2,
  largest_count_cluster_area_fraction,
  largest_count_cluster_states,
  largest_count_cluster_maximum_rainfall_mm,
  largest_count_cluster_maximum_threshold_ratio,
  
  largest_area_cluster_id,
  largest_area_cluster_stations,
  largest_area_cluster_area_km2,
  largest_area_cluster_states,
  
  coherent_clusters,
  isolated_extremes
)])
