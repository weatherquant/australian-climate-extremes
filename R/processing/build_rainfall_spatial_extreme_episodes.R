# ============================================================
# build_rainfall_spatial_extreme_episodes.R
#
# Join daily coherent P97.5 / 200-km rainfall clusters into
# multi-day spatial extreme-rainfall episodes.
#
# Primary daily cluster definition inherited from previous stage:
#   threshold:              P97.5
#   station connectivity:   <= 200 km
#   coherent cluster:       >= 3 stations
#
# Episode linkage between clusters on consecutive calendar days:
#   link if ANY of the following is true:
#     1. at least one station is shared; OR
#     2. station-set Jaccard overlap >= 0.10; OR
#     3. centroid displacement <= 200 km.
#
# Links are permitted only for clusters exactly one day apart.
# Connected components in this inter-day graph define episodes.
# This permits splitting/merging systems while preventing same-day
# independent clusters from being collapsed into one episode.
# ============================================================

library(data.table)
library(arrow)
library(sf)
library(igraph)

# ============================================================
# Configuration
# ============================================================

cluster_file <- file.path(
  "data-derived", "rainfall", "spatial_extremes",
  "rainfall_spatial_extreme_clusters_p975_200km_1950_2025.parquet"
)

daily_cluster_file <- file.path(
  "data-derived", "rainfall", "spatial_extremes",
  "rainfall_spatial_extreme_cluster_days_p975_200km_1950_2025.parquet"
)

output_dir <- file.path(
  "data-derived", "rainfall", "spatial_extremes"
)

table_output_dir <- file.path(
  "tables", "rainfall_spatial_extremes"
)

episode_output_file <- file.path(
  output_dir,
  "rainfall_spatial_extreme_episodes_p975_200km_1950_2025.parquet"
)

episode_day_output_file <- file.path(
  output_dir,
  "rainfall_spatial_extreme_episode_days_p975_200km_1950_2025.parquet"
)

annual_episode_output_file <- file.path(
  output_dir,
  "rainfall_spatial_extreme_episode_annual_p975_200km_1950_2025.parquet"
)

minimum_coherent_stations <- 3L
minimum_jaccard_overlap <- 0.10
maximum_centroid_displacement_km <- 200
maximum_day_gap <- 1L

# ============================================================
# Input checks
# ============================================================

for (f in c(cluster_file, daily_cluster_file)) {
  if (!file.exists(f)) {
    stop("Missing required input: ", f, call. = FALSE)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load coherent daily clusters only
# ============================================================

clusters <- as.data.table(read_parquet(cluster_file))

required_columns <- c(
  "cluster_id", "date", "year", "month", "season",
  "station_count", "coherent_cluster", "represented_area_km2",
  "represented_area_fraction", "states", "centroid_latitude",
  "centroid_longitude", "maximum_rainfall_mm",
  "maximum_threshold_ratio", "station_numbers"
)

missing_columns <- setdiff(required_columns, names(clusters))

if (length(missing_columns) > 0L) {
  stop(
    "Cluster file is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

clusters <- clusters[
  coherent_cluster == TRUE &
    station_count >= minimum_coherent_stations
]

if (nrow(clusters) == 0L) {
  stop("No coherent daily clusters found.", call. = FALSE)
}

setorder(clusters, date, cluster_id)
clusters[, node_id := seq_len(.N)]

# Parse station sets once.
station_sets <- strsplit(clusters$station_numbers, ";", fixed = TRUE)

# ============================================================
# Project centroids for kilometre distances
# ============================================================

centroid_sf <- st_as_sf(
  clusters,
  coords = c("centroid_longitude", "centroid_latitude"),
  crs = 4326,
  remove = FALSE
)
centroid_sf <- st_transform(centroid_sf, 3577)
xy <- st_coordinates(centroid_sf)
clusters[, `:=`(centroid_x_m = xy[, 1], centroid_y_m = xy[, 2])]

# ============================================================
# Build inter-day episode edges
# ============================================================

unique_dates <- sort(unique(clusters$date))
edge_list <- vector("list", max(0L, length(unique_dates) - 1L))
edge_index <- 0L

for (i in seq_len(length(unique_dates) - 1L)) {

  date_a <- unique_dates[i]
  date_b <- unique_dates[i + 1L]

  if (as.integer(date_b - date_a) != maximum_day_gap) {
    next
  }

  idx_a <- which(clusters$date == date_a)
  idx_b <- which(clusters$date == date_b)

  candidate_edges <- list()
  candidate_index <- 0L

  for (a in idx_a) {
    set_a <- station_sets[[a]]

    for (b in idx_b) {
      set_b <- station_sets[[b]]

      shared <- length(intersect(set_a, set_b))
      union_size <- length(union(set_a, set_b))
      jaccard <- if (union_size > 0L) shared / union_size else 0

      dx <- clusters$centroid_x_m[a] - clusters$centroid_x_m[b]
      dy <- clusters$centroid_y_m[a] - clusters$centroid_y_m[b]
      centroid_distance_km <- sqrt(dx^2 + dy^2) / 1000

      linked <-
        shared >= 1L ||
        jaccard >= minimum_jaccard_overlap ||
        centroid_distance_km <= maximum_centroid_displacement_km

      if (linked) {
        candidate_index <- candidate_index + 1L
        candidate_edges[[candidate_index]] <- data.table(
          from = clusters$node_id[a],
          to = clusters$node_id[b],
          shared_stations = as.integer(shared),
          jaccard_overlap = as.numeric(jaccard),
          centroid_distance_km = as.numeric(centroid_distance_km)
        )
      }
    }
  }

  if (candidate_index > 0L) {
    edge_index <- edge_index + 1L
    edge_list[[edge_index]] <- rbindlist(candidate_edges)
  }
}

if (edge_index == 0L) {
  edges <- data.table(
    from = integer(), to = integer(), shared_stations = integer(),
    jaccard_overlap = numeric(), centroid_distance_km = numeric()
  )
} else {
  edges <- rbindlist(edge_list[seq_len(edge_index)])
}

# ============================================================
# Connected components define multi-day episodes
# ============================================================

vertex_table <- data.frame(name = as.character(clusters$node_id))

if (nrow(edges) == 0L) {
  episode_membership <- seq_len(nrow(clusters))
} else {
  graph <- graph_from_data_frame(
    d = data.frame(
      from = as.character(edges$from),
      to = as.character(edges$to)
    ),
    directed = FALSE,
    vertices = vertex_table
  )

  comp <- components(graph)$membership
  episode_membership <- as.integer(comp[as.character(clusters$node_id)])
}

clusters[, episode_component := episode_membership]

# Stable episode IDs ordered by start date and component.
episode_order <- clusters[, .(
  start_date = min(date),
  end_date = max(date)
), by = episode_component]
setorder(episode_order, start_date, end_date, episode_component)
episode_order[, episode_id := sprintf("EP%05d", seq_len(.N))]

clusters <- episode_order[
  clusters,
  on = "episode_component"
]

# ============================================================
# Daily membership table
# ============================================================

episode_days <- clusters[, .(
  episode_id,
  cluster_id,
  date,
  year,
  month,
  season,
  station_count,
  represented_area_km2,
  represented_area_fraction,
  states,
  centroid_latitude,
  centroid_longitude,
  maximum_rainfall_mm,
  maximum_threshold_ratio,
  station_numbers
)]

setorder(episode_days, episode_id, date, -station_count)

# ============================================================
# Episode summaries
# ============================================================

episodes <- clusters[, {

  all_station_numbers <- unique(unlist(strsplit(station_numbers, ";", fixed = TRUE)))
  dates <- sort(unique(date))

  peak_count_index <- which.max(station_count)
  peak_area_index <- which.max(represented_area_km2)

  list(
    start_date = min(date),
    end_date = max(date),
    duration_days = as.integer(max(date) - min(date)) + 1L,
    active_days = as.integer(uniqueN(date)),
    daily_cluster_count = as.integer(.N),
    unique_stations_affected = as.integer(length(all_station_numbers)),

    maximum_simultaneous_stations = as.integer(max(station_count)),
    median_daily_cluster_stations = as.numeric(median(station_count)),

    peak_count_date = date[peak_count_index],
    peak_count_cluster_id = cluster_id[peak_count_index],
    peak_count_states = states[peak_count_index],
    peak_count_area_km2 = as.numeric(represented_area_km2[peak_count_index]),

    maximum_cluster_area_km2 = as.numeric(max(represented_area_km2)),
    peak_area_date = date[peak_area_index],
    peak_area_cluster_id = cluster_id[peak_area_index],
    peak_area_stations = as.integer(station_count[peak_area_index]),
    peak_area_states = states[peak_area_index],

    maximum_rainfall_mm = as.numeric(max(maximum_rainfall_mm, na.rm = TRUE)),
    maximum_threshold_ratio = as.numeric(max(maximum_threshold_ratio, na.rm = TRUE)),

    first_centroid_latitude = centroid_latitude[which.min(date)],
    first_centroid_longitude = centroid_longitude[which.min(date)],
    last_centroid_latitude = centroid_latitude[which.max(date)],
    last_centroid_longitude = centroid_longitude[which.max(date)],

    states_ever_affected = paste(
      sort(unique(unlist(strsplit(states, ";", fixed = TRUE)))),
      collapse = ";"
    ),

    station_numbers_ever = paste(sort(all_station_numbers), collapse = ";")
  )
}, by = episode_id]

setorder(episodes, start_date, episode_id)

# ============================================================
# Attach coverage information by episode day
# ============================================================

daily_coverage <- as.data.table(read_parquet(
  daily_cluster_file,
  col_select = c(
    "date", "coverage_primary", "coverage_strict", "coverage_tier"
  )
))

daily_coverage <- unique(daily_coverage, by = "date")

coverage_by_episode <- merge(
  unique(episode_days[, .(episode_id, date)]),
  daily_coverage,
  by = "date",
  all.x = TRUE
)

coverage_episode_summary <- coverage_by_episode[, .(
  primary_coverage_days = as.integer(sum(coverage_primary %in% TRUE)),
  strict_coverage_days = as.integer(sum(coverage_strict %in% TRUE)),
  low_coverage_days = as.integer(sum(!(coverage_primary %in% TRUE))),
  all_days_primary_coverage = all(coverage_primary %in% TRUE),
  all_days_strict_coverage = all(coverage_strict %in% TRUE)
), by = episode_id]

episodes <- coverage_episode_summary[
  episodes,
  on = "episode_id"
]

# Episodes used in primary inference must have all active cluster days at
# primary 80/80 coverage. Strict 90/90 is retained as sensitivity.

# ============================================================
# Annual episode summary
# ============================================================

episodes[, start_year := as.integer(format(start_date, "%Y"))]

annual_episodes <- episodes[, .(
  episodes_all = as.integer(.N),
  episodes_primary = as.integer(sum(all_days_primary_coverage)),
  episodes_strict = as.integer(sum(all_days_strict_coverage)),
  multi_day_episodes_primary = as.integer(sum(
    all_days_primary_coverage & duration_days >= 2L
  )),
  median_duration_days_primary = as.numeric(median(
    duration_days[all_days_primary_coverage], na.rm = TRUE
  )),
  maximum_duration_days_primary = as.integer(max(
    duration_days[all_days_primary_coverage], na.rm = TRUE
  )),
  median_maximum_simultaneous_stations_primary = as.numeric(median(
    maximum_simultaneous_stations[all_days_primary_coverage], na.rm = TRUE
  )),
  maximum_simultaneous_stations_primary = as.integer(max(
    maximum_simultaneous_stations[all_days_primary_coverage], na.rm = TRUE
  )),
  maximum_cluster_area_km2_primary = as.numeric(max(
    maximum_cluster_area_km2[all_days_primary_coverage], na.rm = TRUE
  ))
), by = start_year]

setorder(annual_episodes, start_year)

# ============================================================
# Diagnostics
# ============================================================

episode_summary <- episodes[, .(
  episodes = .N,
  primary_coverage_episodes = sum(all_days_primary_coverage),
  strict_coverage_episodes = sum(all_days_strict_coverage),
  single_day_episodes = sum(duration_days == 1L),
  multi_day_episodes = sum(duration_days >= 2L),
  median_duration_days = median(duration_days),
  p90_duration_days = quantile(duration_days, 0.90, names = FALSE),
  p99_duration_days = quantile(duration_days, 0.99, names = FALSE),
  median_unique_stations = median(unique_stations_affected),
  p90_unique_stations = quantile(unique_stations_affected, 0.90, names = FALSE),
  p99_unique_stations = quantile(unique_stations_affected, 0.99, names = FALSE),
  median_maximum_simultaneous_stations = median(maximum_simultaneous_stations),
  p90_maximum_simultaneous_stations = quantile(
    maximum_simultaneous_stations, 0.90, names = FALSE
  ),
  p99_maximum_simultaneous_stations = quantile(
    maximum_simultaneous_stations, 0.99, names = FALSE
  )
)]

largest_episodes <- episodes[
  all_days_primary_coverage == TRUE
][
  order(
    -maximum_simultaneous_stations,
    -unique_stations_affected,
    -maximum_cluster_area_km2
  )
][1:min(50L, .N)]

# ============================================================
# Save outputs
# ============================================================

write_parquet(episodes, episode_output_file, compression = "zstd")
write_parquet(episode_days, episode_day_output_file, compression = "zstd")
write_parquet(annual_episodes, annual_episode_output_file, compression = "zstd")

fwrite(
  episode_summary,
  file.path(table_output_dir, "rainfall_spatial_episode_summary.csv")
)

fwrite(
  annual_episodes,
  file.path(table_output_dir, "rainfall_spatial_episode_annual_summary.csv")
)

fwrite(
  largest_episodes,
  file.path(table_output_dir, "rainfall_spatial_episode_largest.csv")
)

fwrite(
  edges,
  file.path(table_output_dir, "rainfall_spatial_episode_link_audit.csv")
)

# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nRAINFALL SPATIAL EXTREME EPISODES COMPLETE",
  "\n========================================"
)

print(episode_summary)

message("\nLargest primary-coverage episodes:")
print(largest_episodes[, .(
  episode_id,
  start_date,
  end_date,
  duration_days,
  active_days,
  daily_cluster_count,
  unique_stations_affected,
  maximum_simultaneous_stations,
  peak_count_date,
  peak_count_states,
  peak_count_area_km2,
  maximum_cluster_area_km2,
  maximum_rainfall_mm,
  maximum_threshold_ratio,
  states_ever_affected
)])
