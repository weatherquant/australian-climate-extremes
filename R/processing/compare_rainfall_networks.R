# ============================================================
# compare_rainfall_networks.R
#
# Construct and compare candidate Australian rainfall networks
# for common periods beginning in 1930, 1950 and 1960.
#
# Default selection rule:
#
#   station-year completeness >= 0.95
#
#   and
#
#   qualifying station-years in >= 90% of all years
#   from the candidate start year through 2025.
#
# Spatial distances are calculated using Australian Albers:
# EPSG:3577.
# ============================================================


# ============================================================
# Validate required packages
# ============================================================

check_rainfall_network_packages <- function() {
  
  required_packages <- c(
    "data.table",
    "sf",
    "RANN"
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
      paste(
        missing_packages,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Build one candidate network
# ============================================================

build_candidate_rainfall_network <- function(
    station_year_metadata,
    start_year,
    end_year = 2025L,
    annual_completeness_threshold = 0.95,
    required_year_fraction = 0.90) {
  
  required_columns <- c(
    "station_number",
    "year",
    "completeness",
    "station_name",
    "latitude",
    "longitude",
    "elevation_m",
    "state",
    "rainfall_district",
    "area_name",
    "open_date",
    "close_date",
    "station_open",
    "valid_coordinates"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(station_year_metadata)
  )
  
  if (length(missing_columns) > 0L) {
    
    stop(
      "Station-year metadata are missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  start_year <- as.integer(start_year)
  end_year <- as.integer(end_year)
  
  if (end_year < start_year) {
    
    stop(
      "end_year must not be earlier than start_year.",
      call. = FALSE
    )
  }
  
  if (
    annual_completeness_threshold < 0 ||
    annual_completeness_threshold > 1
  ) {
    
    stop(
      "annual_completeness_threshold must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (
    required_year_fraction < 0 ||
    required_year_fraction > 1
  ) {
    
    stop(
      "required_year_fraction must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  x <- data.table::copy(
    station_year_metadata
  )
  
  period <- x[
    year >= start_year &
      year <= end_year
  ]
  
  period_years <-
    end_year -
    start_year +
    1L
  
  required_years <- ceiling(
    period_years *
      required_year_fraction
  )
  
  station_period <- period[
    ,
    .(
      years_in_table =
        data.table::uniqueN(year),
      
      qualifying_years =
        sum(
          completeness >=
            annual_completeness_threshold,
          na.rm = TRUE
        ),
      
      mean_completeness =
        mean(
          completeness,
          na.rm = TRUE
        ),
      
      median_completeness =
        stats::median(
          completeness,
          na.rm = TRUE
        ),
      
      minimum_completeness =
        min(
          completeness,
          na.rm = TRUE
        ),
      
      first_qualifying_year =
        if (
          any(
            completeness >=
            annual_completeness_threshold,
            na.rm = TRUE
          )
        ) {
          min(
            year[
              completeness >=
                annual_completeness_threshold
            ],
            na.rm = TRUE
          )
        } else {
          NA_integer_
        },
      
      last_qualifying_year =
        if (
          any(
            completeness >=
            annual_completeness_threshold,
            na.rm = TRUE
          )
        ) {
          max(
            year[
              completeness >=
                annual_completeness_threshold
            ],
            na.rm = TRUE
          )
        } else {
          NA_integer_
        },
      
      station_name =
        data.table::first(
          station_name
        ),
      
      latitude =
        data.table::first(
          latitude
        ),
      
      longitude =
        data.table::first(
          longitude
        ),
      
      elevation_m =
        data.table::first(
          elevation_m
        ),
      
      state =
        data.table::first(
          state
        ),
      
      rainfall_district =
        data.table::first(
          rainfall_district
        ),
      
      area_name =
        data.table::first(
          area_name
        ),
      
      open_date =
        data.table::first(
          open_date
        ),
      
      close_date =
        data.table::first(
          close_date
        ),
      
      station_open =
        data.table::first(
          station_open
        ),
      
      valid_coordinates =
        data.table::first(
          valid_coordinates
        )
    ),
    by = station_number
  ]
  
  station_period[
    ,
    qualifying_fraction :=
      qualifying_years /
      period_years
  ]
  
  station_period[
    ,
    `:=`(
      selection_start_year =
        start_year,
      
      selection_end_year =
        end_year,
      
      period_years =
        period_years,
      
      required_years =
        required_years,
      
      annual_completeness_threshold =
        annual_completeness_threshold,
      
      required_year_fraction =
        required_year_fraction
    )
  ]
  
  station_period[
    ,
    selected_temporally :=
      qualifying_years >=
      required_years
  ]
  
  station_period[
    ,
    selected_spatially :=
      selected_temporally &
      valid_coordinates &
      !is.na(latitude) &
      !is.na(longitude)
  ]
  
  station_period[
    ,
    selection_rule :=
      paste0(
        start_year,
        "-",
        end_year,
        "; annual completeness >= ",
        annual_completeness_threshold,
        "; qualifying in >= ",
        required_year_fraction,
        " of years"
      )
  ]
  
  data.table::setorder(
    station_period,
    station_number
  )
  
  station_period
}


# ============================================================
# Calculate nearest-neighbour distances
# ============================================================

add_nearest_neighbour_distance <- function(
    network,
    selected_column = "selected_spatially") {
  
  check_rainfall_network_packages()
  
  if (!selected_column %in% names(network)) {
    
    stop(
      "Selection column not found: ",
      selected_column,
      call. = FALSE
    )
  }
  
  selected <- data.table::copy(
    network[
      get(selected_column) == TRUE
    ]
  )
  
  if (nrow(selected) < 2L) {
    
    selected[
      ,
      nearest_neighbour_km :=
        NA_real_
    ]
    
    return(selected)
  }
  
  selected_sf <- sf::st_as_sf(
    selected,
    coords = c(
      "longitude",
      "latitude"
    ),
    crs = 4326,
    remove = FALSE
  )
  
  selected_albers <- sf::st_transform(
    selected_sf,
    3577
  )
  
  coordinates <- sf::st_coordinates(
    selected_albers
  )
  
  nearest <- RANN::nn2(
    data = coordinates,
    query = coordinates,
    k = 2
  )
  
  # First neighbour is the station itself.
  selected[
    ,
    nearest_neighbour_km :=
      nearest$nn.dists[, 2L] /
      1000
  ]
  
  selected
}


# ============================================================
# Summarise one candidate network
# ============================================================

summarise_candidate_rainfall_network <- function(
    network) {
  
  temporal <- network[
    selected_temporally == TRUE
  ]
  
  spatial <- network[
    selected_spatially == TRUE
  ]
  
  data.table::data.table(
    
    start_year =
      unique(
        network$selection_start_year
      ),
    
    end_year =
      unique(
        network$selection_end_year
      ),
    
    period_years =
      unique(
        network$period_years
      ),
    
    required_years =
      unique(
        network$required_years
      ),
    
    annual_completeness_threshold =
      unique(
        network$
          annual_completeness_threshold
      ),
    
    required_year_fraction =
      unique(
        network$
          required_year_fraction
      ),
    
    temporally_selected_stations =
      nrow(temporal),
    
    spatially_valid_stations =
      nrow(spatial),
    
    selected_without_valid_coordinates =
      nrow(temporal) -
      nrow(spatial),
    
    currently_open =
      sum(
        spatial$station_open,
        na.rm = TRUE
      ),
    
    closed =
      sum(
        !spatial$station_open,
        na.rm = TRUE
      ),
    
    missing_elevation =
      sum(
        is.na(
          spatial$elevation_m
        )
      ),
    
    minimum_latitude =
      min(
        spatial$latitude,
        na.rm = TRUE
      ),
    
    maximum_latitude =
      max(
        spatial$latitude,
        na.rm = TRUE
      ),
    
    minimum_longitude =
      min(
        spatial$longitude,
        na.rm = TRUE
      ),
    
    maximum_longitude =
      max(
        spatial$longitude,
        na.rm = TRUE
      ),
    
    median_elevation_m =
      stats::median(
        spatial$elevation_m,
        na.rm = TRUE
      ),
    
    q10_elevation_m =
      stats::quantile(
        spatial$elevation_m,
        0.10,
        na.rm = TRUE,
        names = FALSE
      ),
    
    q90_elevation_m =
      stats::quantile(
        spatial$elevation_m,
        0.90,
        na.rm = TRUE,
        names = FALSE
      )
  )
}


# ============================================================
# State-level summary
# ============================================================

summarise_candidate_network_by_state <- function(
    network) {
  
  network[
    selected_spatially == TRUE,
    .(
      stations =
        .N,
      
      currently_open =
        sum(
          station_open,
          na.rm = TRUE
        ),
      
      closed =
        sum(
          !station_open,
          na.rm = TRUE
        ),
      
      missing_elevation =
        sum(
          is.na(
            elevation_m
          )
        ),
      
      median_elevation_m =
        stats::median(
          elevation_m,
          na.rm = TRUE
        ),
      
      median_completeness =
        stats::median(
          mean_completeness,
          na.rm = TRUE
        ),
      
      median_qualifying_fraction =
        stats::median(
          qualifying_fraction,
          na.rm = TRUE
        )
    ),
    by = state
  ][
    order(state)
  ]
}


# ============================================================
# Rainfall-district summary
# ============================================================

summarise_candidate_network_by_district <- function(
    network) {
  
  network[
    selected_spatially == TRUE,
    .(
      stations =
        .N,
      
      currently_open =
        sum(
          station_open,
          na.rm = TRUE
        ),
      
      closed =
        sum(
          !station_open,
          na.rm = TRUE
        ),
      
      median_qualifying_fraction =
        stats::median(
          qualifying_fraction,
          na.rm = TRUE
        )
    ),
    by = .(
      state,
      rainfall_district
    )
  ][
    order(
      state,
      rainfall_district
    )
  ]
}


# ============================================================
# Compare multiple candidate start years
# ============================================================

compare_candidate_rainfall_networks <- function(
    station_year_metadata,
    start_years = c(
      1930L,
      1950L,
      1960L
    ),
    end_year = 2025L,
    annual_completeness_threshold = 0.95,
    required_year_fraction = 0.90) {
  
  check_rainfall_network_packages()
  
  networks <- lapply(
    start_years,
    function(start_year) {
      
      network <- build_candidate_rainfall_network(
        
        station_year_metadata =
          station_year_metadata,
        
        start_year =
          start_year,
        
        end_year =
          end_year,
        
        annual_completeness_threshold =
          annual_completeness_threshold,
        
        required_year_fraction =
          required_year_fraction
      )
      
      add_nearest_neighbour_distance(
        network
      )
    }
  )
  
  names(networks) <- paste0(
    "network_",
    start_years,
    "_",
    end_year
  )
  
  summaries <- data.table::rbindlist(
    lapply(
      networks,
      summarise_candidate_rainfall_network
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  nearest_neighbour_summary <-
    data.table::rbindlist(
      lapply(
        networks,
        function(network) {
          
          data.table::data.table(
            
            start_year =
              unique(
                network$
                  selection_start_year
              ),
            
            end_year =
              unique(
                network$
                  selection_end_year
              ),
            
            stations =
              nrow(network),
            
            median_nearest_neighbour_km =
              stats::median(
                network$
                  nearest_neighbour_km,
                na.rm = TRUE
              ),
            
            q10_nearest_neighbour_km =
              stats::quantile(
                network$
                  nearest_neighbour_km,
                0.10,
                na.rm = TRUE,
                names = FALSE
              ),
            
            q90_nearest_neighbour_km =
              stats::quantile(
                network$
                  nearest_neighbour_km,
                0.90,
                na.rm = TRUE,
                names = FALSE
              ),
            
            maximum_nearest_neighbour_km =
              max(
                network$
                  nearest_neighbour_km,
                na.rm = TRUE
              )
          )
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
  
  summaries <- merge(
    summaries,
    nearest_neighbour_summary,
    by = c(
      "start_year",
      "end_year"
    ),
    all.x = TRUE
  )
  
  state_summaries <-
    data.table::rbindlist(
      lapply(
        networks,
        function(network) {
          
          result <-
            summarise_candidate_network_by_state(
              network
            )
          
          result[
            ,
            start_year :=
              unique(
                network$
                  selection_start_year
              )
          ]
          
          result[
            ,
            end_year :=
              unique(
                network$
                  selection_end_year
              )
          ]
          
          result
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
  
  district_summaries <-
    data.table::rbindlist(
      lapply(
        networks,
        function(network) {
          
          result <-
            summarise_candidate_network_by_district(
              network
            )
          
          result[
            ,
            start_year :=
              unique(
                network$
                  selection_start_year
              )
          ]
          
          result[
            ,
            end_year :=
              unique(
                network$
                  selection_end_year
              )
          ]
          
          result
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
  
  list(
    networks =
      networks,
    
    summary =
      summaries,
    
    state_summary =
      state_summaries,
    
    district_summary =
      district_summaries
  )
}


# ============================================================
# Compare network membership between periods
# ============================================================

compare_candidate_network_membership <- function(
    comparison) {
  
  networks <- comparison$networks
  
  if (length(networks) < 2L) {
    
    stop(
      "At least two candidate networks are required.",
      call. = FALSE
    )
  }
  
  membership <- data.table::rbindlist(
    lapply(
      names(networks),
      function(network_name) {
        
        network <- networks[[network_name]]
        
        network[
          ,
          .(
            station_number,
            candidate =
              network_name,
            
            selected =
              selected_spatially
          )
        ]
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  membership_wide <- data.table::dcast(
    membership,
    station_number ~ candidate,
    value.var = "selected",
    fill = FALSE
  )
  
  selection_columns <- setdiff(
    names(membership_wide),
    "station_number"
  )
  
  membership_wide[
    ,
    networks_selected :=
      rowSums(
        .SD
      ),
    .SDcols =
      selection_columns
  ]
  
  membership_wide
}