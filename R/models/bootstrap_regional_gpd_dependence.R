# ============================================================
# bootstrap_regional_gpd_dependence.R
#
# Dependence-aware block-bootstrap pilot for regional rainfall
# GPD temporal signals.
#
# Primary event definition:
#   P97.5 threshold
#   3-day runs declustering
#   PRCP_QUAL 0-4
#   1950-2025
#   NSW, QLD, VIC
#
# Bootstrap designs:
#   M1: circular moving blocks over the full 1950-2025 series.
#       Source-year blocks are mapped to ordered target years so
#       the linear temporal estimand remains defined.
#
#   M2: circular moving blocks sampled separately within each
#       fixed climate period; original period labels are retained.
#
# Default is a 50-replicate, 3-year-block pilot. After confirming
# runtime and >=95% successful fits, set n_bootstrap to 200.
#
# The script saves each successful or failed replicate so the
# Quarto report never reruns the bootstrap.
# ============================================================

library(data.table)
library(arrow)

# ============================================================
# Configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L
centre_year <- mean(c(analysis_start_year, analysis_end_year))

regions <- c("NSW", "QLD", "VIC")
accepted_rainfall_quality <- c(0L, 1L, 2L, 3L, 4L)

n_bootstrap <- 200L
block_length_years <- 3L
random_seed <- 20261112L

# Save interim outputs every this many replicates per region.
checkpoint_every <- 10L

pot_file <- file.path(
  "data-derived", "rainfall", "pot",
  "rainfall_pot_events_p975_run3_1950_2025.parquet"
)

helper_file <- file.path(
  "R", "models", "gpd_helpers.R"
)

model_output_dir <- file.path(
  "data-derived", "rainfall", "models", "gpd"
)

table_output_dir <- file.path(
  "tables", "rainfall_models", "gpd"
)

dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_file <- file.path(
  model_output_dir,
  paste0(
    "regional_gpd_dependence_bootstrap_checkpoint_b",
    block_length_years,
    "_n",
    n_bootstrap,
    ".rds"
  )
)

final_rds_file <- file.path(
  model_output_dir,
  paste0(
    "regional_gpd_dependence_bootstrap_b",
    block_length_years,
    "_n",
    n_bootstrap,
    ".rds"
  )
)

if (!file.exists(pot_file)) {
  stop("Missing primary POT file: ", pot_file, call. = FALSE)
}

if (!file.exists(helper_file)) {
  stop("Missing GPD helper file: ", helper_file, call. = FALSE)
}

source(helper_file)

# ============================================================
# Bootstrap sampling functions
# ============================================================

sample_circular_blocks <- function(
    source_years,
    target_years,
    block_length) {

  source_years <- as.integer(source_years)
  target_years <- as.integer(target_years)
  block_length <- as.integer(block_length)

  if (length(source_years) != length(target_years)) {
    stop("source_years and target_years must have equal length.", call. = FALSE)
  }

  n_years <- length(source_years)

  if (block_length < 1L || block_length > n_years) {
    stop("Invalid block length.", call. = FALSE)
  }

  sampled_source <- integer(0L)

  while (length(sampled_source) < n_years) {
    start_index <- sample.int(n_years, size = 1L)
    block_index <- ((start_index - 1L + 0:(block_length - 1L)) %% n_years) + 1L
    sampled_source <- c(sampled_source, source_years[block_index])
  }

  sampled_source <- sampled_source[seq_len(n_years)]

  data.table(
    target_year = target_years,
    source_year = sampled_source,
    draw_position = seq_len(n_years)
  )
}

make_m1_year_map <- function(
    years,
    block_length) {

  years <- sort(unique(as.integer(years)))

  sample_circular_blocks(
    source_years = years,
    target_years = years,
    block_length = block_length
  )
}

make_m2_year_map <- function(
    block_length) {

  periods <- list(
    `1950-1975` = 1950:1975,
    `1976-2000` = 1976:2000,
    `2001-2025` = 2001:2025
  )

  result <- rbindlist(
    lapply(
      names(periods),
      function(period_name) {
        years <- periods[[period_name]]
        out <- sample_circular_blocks(
          source_years = years,
          target_years = years,
          block_length = min(block_length, length(years))
        )
        out[, climate_period := period_name]
        out
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  setorder(result, target_year)
  result
}

expand_bootstrap_events <- function(
    events_region,
    year_map,
    model_type = c("M1", "M2")) {

  model_type <- match.arg(model_type)

  # Each source year may appear multiple times in the bootstrap.
  # The merge therefore duplicates all same-year station events together.
  out <- merge(
    year_map,
    events_region,
    by.x = "source_year",
    by.y = "year",
    allow.cartesian = TRUE,
    sort = FALSE
  )

  # Retain the sampled source year as the M1 covariate. This is a
  # moving-block pairs bootstrap: whole region-year clusters, including
  # their original time coordinates, are sampled together. target_year
  # records bootstrap ordering only.
  out[, bootstrap_year := as.integer(target_year)]
  out[, year_decade_boot := (as.integer(source_year) - centre_year) / 10]

  if (model_type == "M2") {
    out[, climate_period_boot := climate_period]
  }

  out
}

# ============================================================
# Safe fit wrappers
# ============================================================

safe_fit_m1 <- function(x) {

  started <- proc.time()[["elapsed"]]

  fit <- tryCatch(
    fit_gpd_m1(
      excess = x$normalized_excess,
      year_decade = x$year_decade_boot
    ),
    error = function(e) e
  )

  elapsed <- proc.time()[["elapsed"]] - started

  if (inherits(fit, "error")) {
    return(list(
      success = FALSE,
      error = conditionMessage(fit),
      elapsed_seconds = elapsed,
      fit = NULL
    ))
  }

  list(
    success = isTRUE(fit$fit_valid),
    error = if (isTRUE(fit$fit_valid)) "" else "fit_valid was FALSE",
    elapsed_seconds = elapsed,
    fit = fit
  )
}

safe_fit_m2 <- function(x) {

  started <- proc.time()[["elapsed"]]

  fit <- tryCatch(
    fit_gpd_m2(
      excess = x$normalized_excess,
      climate_period = x$climate_period_boot
    ),
    error = function(e) e
  )

  elapsed <- proc.time()[["elapsed"]] - started

  if (inherits(fit, "error")) {
    return(list(
      success = FALSE,
      error = conditionMessage(fit),
      elapsed_seconds = elapsed,
      fit = NULL
    ))
  }

  list(
    success = isTRUE(fit$fit_valid),
    error = if (isTRUE(fit$fit_valid)) "" else "fit_valid was FALSE",
    elapsed_seconds = elapsed,
    fit = fit
  )
}

extract_m1_bootstrap <- function(region, replicate_id, fit_result, n_events) {

  if (!fit_result$success || is.null(fit_result$fit)) {
    return(data.table(
      region = region,
      replicate = as.integer(replicate_id),
      model = "M1",
      success = FALSE,
      events = as.integer(n_events),
      elapsed_seconds = as.numeric(fit_result$elapsed_seconds),
      error = fit_result$error,
      beta_log_scale_decade = NA_real_,
      percent_scale_change_decade = NA_real_,
      shape = NA_real_,
      convergence = NA_integer_,
      hessian_positive_definite = FALSE,
      maximum_gradient = NA_real_
    ))
  }

  fit <- fit_result$fit
  beta <- as.numeric(fit$parameters_raw["scale_time_per_decade"])
  shape <- as.numeric(shape_from_raw(fit$parameters_raw["raw_shape"]))

  data.table(
    region = region,
    replicate = as.integer(replicate_id),
    model = "M1",
    success = TRUE,
    events = as.integer(n_events),
    elapsed_seconds = as.numeric(fit_result$elapsed_seconds),
    error = "",
    beta_log_scale_decade = beta,
    percent_scale_change_decade = 100 * (exp(beta) - 1),
    shape = shape,
    convergence = as.integer(fit$convergence),
    hessian_positive_definite = fit$hessian_positive_definite,
    maximum_gradient = as.numeric(fit$maximum_gradient)
  )
}

extract_m2_bootstrap <- function(region, replicate_id, fit_result, n_events) {

  if (!fit_result$success || is.null(fit_result$fit)) {
    return(data.table(
      region = region,
      replicate = as.integer(replicate_id),
      model = "M2",
      success = FALSE,
      events = as.integer(n_events),
      elapsed_seconds = as.numeric(fit_result$elapsed_seconds),
      error = fit_result$error,
      middle_log_scale_difference = NA_real_,
      recent_log_scale_difference = NA_real_,
      middle_vs_early_pct = NA_real_,
      recent_vs_early_pct = NA_real_,
      shape = NA_real_,
      convergence = NA_integer_,
      hessian_positive_definite = FALSE,
      maximum_gradient = NA_real_
    ))
  }

  fit <- fit_result$fit
  middle <- as.numeric(fit$parameters_raw["period_1976_2000"])
  recent <- as.numeric(fit$parameters_raw["period_2001_2025"])
  shape <- as.numeric(shape_from_raw(fit$parameters_raw["raw_shape"]))

  data.table(
    region = region,
    replicate = as.integer(replicate_id),
    model = "M2",
    success = TRUE,
    events = as.integer(n_events),
    elapsed_seconds = as.numeric(fit_result$elapsed_seconds),
    error = "",
    middle_log_scale_difference = middle,
    recent_log_scale_difference = recent,
    middle_vs_early_pct = 100 * (exp(middle) - 1),
    recent_vs_early_pct = 100 * (exp(recent) - 1),
    shape = shape,
    convergence = as.integer(fit$convergence),
    hessian_positive_definite = fit$hessian_positive_definite,
    maximum_gradient = as.numeric(fit$maximum_gradient)
  )
}

# ============================================================
# Load primary event data
# ============================================================

events <- as.data.table(arrow::read_parquet(pot_file))

required_columns <- c(
  "station_number", "state", "event_peak_date",
  "threshold_mm", "excess_mm", "rainfall_quality_peak"
)

missing_columns <- setdiff(required_columns, names(events))

if (length(missing_columns) > 0L) {
  stop(
    "Primary POT file is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

events <- events[
  state %in% regions &
    rainfall_quality_peak %in% accepted_rainfall_quality &
    is.finite(excess_mm) & excess_mm > 0 &
    is.finite(threshold_mm) & threshold_mm > 0
]

events[, year := as.integer(format(event_peak_date, "%Y"))]
events <- events[year >= analysis_start_year & year <= analysis_end_year]
events[, normalized_excess := excess_mm / threshold_mm]
events[, climate_period := assign_climate_period(year)]

if (nrow(events) == 0L) {
  stop("No eligible regional events were found.", call. = FALSE)
}

# Keep only columns needed repeatedly in bootstrap merges.
events_small <- events[, .(
  station_number,
  state,
  year,
  normalized_excess
)]

# ============================================================
# Original pooled regional fits for reference
# ============================================================

message("Fitting original pooled regional M1/M2 reference models")

pooled_reference_list <- lapply(
  regions,
  function(region_value) {

    x <- events[state == region_value]
    x[, year_decade := (year - centre_year) / 10]

    m1 <- fit_gpd_m1(
      excess = x$normalized_excess,
      year_decade = x$year_decade
    )

    m2 <- fit_gpd_m2(
      excess = x$normalized_excess,
      climate_period = x$climate_period
    )

    m1_par <- extract_m1_parameters(m1)
    m1_time <- m1_par[term == "scale_time_per_decade"]

    middle <- as.numeric(m2$parameters_raw["period_1976_2000"])
    recent <- as.numeric(m2$parameters_raw["period_2001_2025"])

    data.table(
      region = region_value,
      events = as.integer(nrow(x)),
      stations = as.integer(uniqueN(x$station_number)),
      pooled_m1_beta = as.numeric(m1_time$estimate),
      pooled_m1_se = as.numeric(m1_time$std_error),
      pooled_m1_ci_low_pct = 100 * (exp(m1_time$conf_low) - 1),
      pooled_m1_ci_high_pct = 100 * (exp(m1_time$conf_high) - 1),
      pooled_m1_change_pct = 100 * (exp(m1_time$estimate) - 1),
      pooled_m2_middle_pct = 100 * (exp(middle) - 1),
      pooled_m2_recent_pct = 100 * (exp(recent) - 1),
      pooled_m1_valid = m1$fit_valid,
      pooled_m2_valid = m2$fit_valid
    )
  }
)

pooled_reference <- rbindlist(
  pooled_reference_list,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# Bootstrap
# ============================================================

set.seed(random_seed)

m1_results <- vector("list", length(regions) * n_bootstrap)
m2_results <- vector("list", length(regions) * n_bootstrap)

result_index <- 0L
run_started <- proc.time()[["elapsed"]]

for (region_index in seq_along(regions)) {

  region_value <- regions[region_index]
  region_events <- events_small[state == region_value]

  message(
    "\n========================================",
    "\nDependence bootstrap region ", region_index, " of ", length(regions),
    ": ", region_value,
    "\n========================================"
  )

  for (b in seq_len(n_bootstrap)) {

    result_index <- result_index + 1L

    # M1 full-series moving blocks.
    m1_map <- make_m1_year_map(
      years = analysis_start_year:analysis_end_year,
      block_length = block_length_years
    )

    m1_events <- expand_bootstrap_events(
      events_region = region_events,
      year_map = m1_map,
      model_type = "M1"
    )

    m1_fit <- safe_fit_m1(m1_events)

    m1_results[[result_index]] <- extract_m1_bootstrap(
      region = region_value,
      replicate_id = b,
      fit_result = m1_fit,
      n_events = nrow(m1_events)
    )

    # M2 period-stratified moving blocks.
    m2_map <- make_m2_year_map(
      block_length = block_length_years
    )

    m2_events <- expand_bootstrap_events(
      events_region = region_events,
      year_map = m2_map,
      model_type = "M2"
    )

    m2_fit <- safe_fit_m2(m2_events)

    m2_results[[result_index]] <- extract_m2_bootstrap(
      region = region_value,
      replicate_id = b,
      fit_result = m2_fit,
      n_events = nrow(m2_events)
    )

    if (
      b %% checkpoint_every == 0L ||
      b == n_bootstrap
    ) {

      completed_m1 <- rbindlist(
        m1_results[seq_len(result_index)],
        use.names = TRUE,
        fill = TRUE
      )

      completed_m2 <- rbindlist(
        m2_results[seq_len(result_index)],
        use.names = TRUE,
        fill = TRUE
      )

      saveRDS(
        list(
          configuration = list(
            n_bootstrap = n_bootstrap,
            block_length_years = block_length_years,
            random_seed = random_seed,
            regions = regions,
            completed_region = region_value,
            completed_replicate = b
          ),
          pooled_reference = pooled_reference,
          M1 = completed_m1,
          M2 = completed_m2
        ),
        checkpoint_file,
        compress = "xz"
      )

      message(
        "  Completed ",
        region_value,
        " replicate ",
        b,
        " of ",
        n_bootstrap
      )
    }

    rm(m1_map, m1_events, m1_fit, m2_map, m2_events, m2_fit)
    invisible(gc(FALSE))
  }
}

run_elapsed_seconds <- proc.time()[["elapsed"]] - run_started

m1_bootstrap <- rbindlist(m1_results, use.names = TRUE, fill = TRUE)
m2_bootstrap <- rbindlist(m2_results, use.names = TRUE, fill = TRUE)

# ============================================================
# Summaries
# ============================================================

m1_summary <- m1_bootstrap[, .(
  requested_replicates = as.integer(n_bootstrap),
  successful_replicates = as.integer(sum(success)),
  failed_replicates = as.integer(sum(!success)),
  success_fraction = as.numeric(mean(success)),
  median_events = as.numeric(median(events)),
  median_fit_seconds = as.numeric(median(elapsed_seconds)),
  bootstrap_mean_pct = as.numeric(mean(percent_scale_change_decade[success], na.rm = TRUE)),
  bootstrap_median_pct = as.numeric(median(percent_scale_change_decade[success], na.rm = TRUE)),
  bootstrap_sd_log_beta = as.numeric(sd(beta_log_scale_decade[success], na.rm = TRUE)),
  bootstrap_sd_pct = as.numeric(sd(percent_scale_change_decade[success], na.rm = TRUE)),
  bootstrap_ci_low_pct = as.numeric(quantile(
    percent_scale_change_decade[success], 0.025,
    na.rm = TRUE, names = FALSE
  )),
  bootstrap_ci_high_pct = as.numeric(quantile(
    percent_scale_change_decade[success], 0.975,
    na.rm = TRUE, names = FALSE
  )),
  proportion_negative = as.numeric(mean(
    percent_scale_change_decade[success] < 0,
    na.rm = TRUE
  )),
  proportion_positive = as.numeric(mean(
    percent_scale_change_decade[success] > 0,
    na.rm = TRUE
  ))
), by = region]

m1_summary <- merge(
  m1_summary,
  pooled_reference[, .(
    region,
    pooled_change_pct = pooled_m1_change_pct,
    pooled_ci_low_pct = pooled_m1_ci_low_pct,
    pooled_ci_high_pct = pooled_m1_ci_high_pct,
    pooled_se_log_beta = pooled_m1_se
  )],
  by = "region",
  all.x = TRUE,
  sort = FALSE
)

m1_summary[, variance_inflation_ratio :=
  bootstrap_sd_log_beta / pooled_se_log_beta
]

m1_summary[, sign_stable_95 :=
  fifelse(
    region == "NSW",
    bootstrap_ci_high_pct < 0,
    fifelse(
      region == "VIC",
      bootstrap_ci_low_pct > 0,
      bootstrap_ci_low_pct <= 0 & bootstrap_ci_high_pct >= 0
    )
  )
]

m2_summary_long <- rbindlist(
  list(
    m2_bootstrap[success == TRUE, .(
      comparison = "1976-2000 vs 1950-1975",
      estimate_pct = middle_vs_early_pct
    ), by = .(region, replicate)],
    m2_bootstrap[success == TRUE, .(
      comparison = "2001-2025 vs 1950-1975",
      estimate_pct = recent_vs_early_pct
    ), by = .(region, replicate)]
  ),
  use.names = TRUE,
  fill = TRUE
)

m2_summary <- m2_summary_long[, .(
  successful_replicates = as.integer(.N),
  bootstrap_mean_pct = as.numeric(mean(estimate_pct)),
  bootstrap_median_pct = as.numeric(median(estimate_pct)),
  bootstrap_sd_pct = as.numeric(sd(estimate_pct)),
  bootstrap_ci_low_pct = as.numeric(quantile(
    estimate_pct, 0.025,
    names = FALSE
  )),
  bootstrap_ci_high_pct = as.numeric(quantile(
    estimate_pct, 0.975,
    names = FALSE
  )),
  proportion_negative = as.numeric(mean(estimate_pct < 0)),
  proportion_positive = as.numeric(mean(estimate_pct > 0))
), by = .(region, comparison)]

m2_success_summary <- m2_bootstrap[, .(
  requested_replicates = as.integer(n_bootstrap),
  successful_replicates = as.integer(sum(success)),
  failed_replicates = as.integer(sum(!success)),
  success_fraction = as.numeric(mean(success)),
  median_fit_seconds = as.numeric(median(elapsed_seconds))
), by = region]

run_summary <- data.table(
  n_bootstrap = n_bootstrap,
  block_length_years = block_length_years,
  random_seed = random_seed,
  regions = paste(regions, collapse = ","),
  total_elapsed_seconds = as.numeric(run_elapsed_seconds),
  total_elapsed_minutes = as.numeric(run_elapsed_seconds / 60),
  m1_success_fraction = mean(m1_bootstrap$success),
  m2_success_fraction = mean(m2_bootstrap$success),
  completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

# ============================================================
# Save final outputs
# ============================================================

saveRDS(
  list(
    configuration = list(
      analysis_start_year = analysis_start_year,
      analysis_end_year = analysis_end_year,
      centre_year = centre_year,
      regions = regions,
      accepted_rainfall_quality = accepted_rainfall_quality,
      n_bootstrap = n_bootstrap,
      block_length_years = block_length_years,
      random_seed = random_seed
    ),
    pooled_reference = pooled_reference,
    M1 = m1_bootstrap,
    M2 = m2_bootstrap,
    M1_summary = m1_summary,
    M2_summary = m2_summary,
    M2_success_summary = m2_success_summary,
    run_summary = run_summary
  ),
  final_rds_file,
  compress = "xz"
)

prefix <- paste0(
  "regional_gpd_dependence_bootstrap_b",
  block_length_years,
  "_n",
  n_bootstrap
)

fwrite(
  pooled_reference,
  file.path(table_output_dir, paste0(prefix, "_pooled_reference.csv"))
)

fwrite(
  m1_bootstrap,
  file.path(table_output_dir, paste0(prefix, "_m1_replicates.csv"))
)

fwrite(
  m2_bootstrap,
  file.path(table_output_dir, paste0(prefix, "_m2_replicates.csv"))
)

fwrite(
  m1_summary,
  file.path(table_output_dir, paste0(prefix, "_m1_summary.csv"))
)

fwrite(
  m2_summary,
  file.path(table_output_dir, paste0(prefix, "_m2_summary.csv"))
)

fwrite(
  m2_success_summary,
  file.path(table_output_dir, paste0(prefix, "_m2_success.csv"))
)

fwrite(
  run_summary,
  file.path(table_output_dir, paste0(prefix, "_run_summary.csv"))
)

message(
  "\n========================================",
  "\nREGIONAL GPD DEPENDENCE BOOTSTRAP COMPLETE",
  "\n========================================"
)

print(run_summary)
print(m1_summary)
print(m2_success_summary)
print(m2_summary)
