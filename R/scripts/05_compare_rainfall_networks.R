# ============================================================
# 05_compare_rainfall_networks.R
#
# Compare 1930–2025, 1950–2025 and 1960–2025 candidate
# rainfall networks and save reproducible selection artifacts.
# ============================================================


library(data.table)
library(arrow)
library(sf)
library(RANN)

source(
  "R/processing/compare_rainfall_networks.R"
)

# ============================================================
# Input
# ============================================================

station_year_metadata_file <- file.path(
  "data-derived",
  "rainfall",
  "summaries",
  "adam_rainfall_station_year_metadata_1900_2026.parquet"
)

if (!file.exists(
  station_year_metadata_file
)) {
  
  stop(
    "Metadata-enhanced station-year file does not exist: ",
    station_year_metadata_file
  )
}

station_year_metadata <-
  data.table::as.data.table(
    arrow::read_parquet(
      station_year_metadata_file
    )
  )


# ============================================================
# Candidate network comparison
# ============================================================

rainfall_network_comparison <-
  compare_candidate_rainfall_networks(
    
    station_year_metadata =
      station_year_metadata,
    
    start_years =
      c(
        1930L,
        1950L,
        1960L
      ),
    
    end_year =
      2025L,
    
    annual_completeness_threshold =
      0.95,
    
    required_year_fraction =
      0.90
  )


print(
  rainfall_network_comparison$summary
)


print(
  rainfall_network_comparison$state_summary
)


# ============================================================
# Membership comparison
# ============================================================

rainfall_network_membership <-
  compare_candidate_network_membership(
    rainfall_network_comparison
  )


# ============================================================
# Output directories
# ============================================================

network_output_dir <- file.path(
  "data-derived",
  "rainfall",
  "network"
)


table_output_dir <- file.path(
  "tables",
  "rainfall_network"
)


dir.create(
  network_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Save each candidate network
# ============================================================

for (
  network_name in
  names(
    rainfall_network_comparison$networks
  )
) {
  
  network <-
    rainfall_network_comparison$
    networks[[network_name]]
  
  output_name <- paste0(
    network_name,
    ".parquet"
  )
  
  arrow::write_parquet(
    network,
    file.path(
      network_output_dir,
      output_name
    ),
    compression = "zstd"
  )
}


# ============================================================
# Save comparison tables
# ============================================================

data.table::fwrite(
  
  rainfall_network_comparison$summary,
  
  file.path(
    table_output_dir,
    "candidate_network_comparison.csv"
  )
)


data.table::fwrite(
  
  rainfall_network_comparison$state_summary,
  
  file.path(
    table_output_dir,
    "candidate_network_by_state.csv"
  )
)


data.table::fwrite(
  
  rainfall_network_comparison$district_summary,
  
  file.path(
    table_output_dir,
    "candidate_network_by_rainfall_district.csv"
  )
)


data.table::fwrite(
  
  rainfall_network_membership,
  
  file.path(
    table_output_dir,
    "candidate_network_membership.csv"
  )
)


# ============================================================
# Save provisional 1950–2025 core network
#
# This remains provisional until the maps are assessed.
# ============================================================

rainfall_core_network <-
  rainfall_network_comparison$networks[[
    "network_1950_2025"
  ]]

rainfall_core_network[
  ,
  continental_core :=
    state %in% c(
      "NSW",
      "NT",
      "QLD",
      "SA",
      "TAS",
      "VIC",
      "WA",
      "ACT"
    ) &
    latitude >= -44 &
    latitude <= -10 &
    longitude >= 112 &
    longitude <= 154
]

arrow::write_parquet(
  
  rainfall_core_network,
  
  file.path(
    network_output_dir,
    "rainfall_core_network_1950_2025.parquet"
  ),
  
  compression = "zstd"
)


# ============================================================
# Final summary
# ============================================================

message(
  "\n========================================",
  "\nRAINFALL NETWORK COMPARISON COMPLETE",
  "\n========================================"
)


print(
  rainfall_network_comparison$summary
)