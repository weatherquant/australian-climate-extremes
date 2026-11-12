# ============================================================
# fit_regional_gpd_threshold_sensitivity.R
#
# Fast regional threshold sensitivity for the three regions
# showing material temporal structure in the primary P97.5 run3
# analysis: NSW, QLD, VIC.
#
# Threshold sensitivity only:
#   P95 / 3-day
#   P97.5 / 3-day
#   P99 / 3-day
#
# Fits in each region/specification:
#   M0: stationary scale, constant shape
#   M1: linear scale trend per decade, constant shape
#   M2: three-period scale, constant shape
#
# Quality rule: PRCP_QUAL 0-4
# Analysis period: 1950-2025
# ============================================================

library(data.table)
library(arrow)

analysis_start_year <- 1950L
analysis_end_year <- 2025L
centre_year <- mean(c(analysis_start_year, analysis_end_year))
accepted_rainfall_quality <- c(0L, 1L, 2L, 3L, 4L)

regions <- c("NSW", "QLD", "VIC")
thresholds <- c("p95", "p975", "p99")
run_length_days <- 3L

all_events_file <- file.path(
  "data-derived", "rainfall", "pot",
  "rainfall_pot_events_all_thresholds_runs_1950_2025.parquet"
)

helper_file <- file.path(
  "R",
  "models",
  "gpd_helpers.R"
)

model_output_dir <- file.path(
  "data-derived", "rainfall", "models", "gpd"
)

table_output_dir <- file.path(
  "tables", "rainfall_models", "gpd"
)

dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(all_events_file)) {
  stop("Missing POT event file: ", all_events_file, call. = FALSE)
}

if (!file.exists(helper_file)) {
  stop("Missing GPD helper script: ", helper_file, call. = FALSE)
}

# ------------------------------------------------------------
# Load reusable GPD functions only
# ------------------------------------------------------------

source(
  helper_file
)

required_functions <- c(
  "fit_gpd_m0", "fit_gpd_m1", "fit_gpd_m2",
  "extract_m1_parameters", "shape_from_raw",
  "assign_climate_period"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Missing helper functions: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

pot <- as.data.table(arrow::read_parquet(all_events_file))

pot <- pot[
  rainfall_quality_peak %in% accepted_rainfall_quality &
    run_length_days == 3L &
    threshold_id %in% thresholds &
    state %in% regions &
    is.finite(excess_mm) & excess_mm > 0 &
    is.finite(threshold_mm) & threshold_mm > 0
]

pot[, year := as.integer(format(event_peak_date, "%Y"))]
pot <- pot[year >= analysis_start_year & year <= analysis_end_year]
pot[, year_decade := (year - centre_year) / 10]
pot[, climate_period := assign_climate_period(year)]
pot[, normalized_excess := excess_mm / threshold_mm]

# ------------------------------------------------------------
# Design and audit
# ------------------------------------------------------------

design <- CJ(
  region = regions,
  threshold_id = thresholds,
  unique = TRUE
)

design[, threshold_order := match(threshold_id, thresholds)]
design[, region_order := match(region, regions)]
setorder(design, region_order, threshold_order)
design[, c("threshold_order", "region_order") := NULL]

sample_audit <- pot[, .(
  events = as.integer(.N),
  stations = as.integer(uniqueN(station_number)),
  events_1950_1975 = as.integer(sum(climate_period == "1950-1975")),
  events_1976_2000 = as.integer(sum(climate_period == "1976-2000")),
  events_2001_2025 = as.integer(sum(climate_period == "2001-2025"))
), by = .(state, threshold_id)]

# ------------------------------------------------------------
# Result extractor
# ------------------------------------------------------------

extract_result <- function(region, threshold_id, x, m0, m1, m2) {

  m1_par <- extract_m1_parameters(m1)
  m1_time <- m1_par[term == "scale_time_per_decade"]

  beta <- as.numeric(m1_time$estimate)
  beta_low <- as.numeric(m1_time$conf_low)
  beta_high <- as.numeric(m1_time$conf_high)

  middle <- as.numeric(m2$parameters_raw["period_1976_2000"])
  recent <- as.numeric(m2$parameters_raw["period_2001_2025"])

  lr1 <- 2 * (m1$log_likelihood - m0$log_likelihood)
  lr2 <- 2 * (m2$log_likelihood - m0$log_likelihood)

  shape0 <- as.numeric(shape_from_raw(m0$parameters_raw["raw_shape"]))

  data.table(
    region = region,
    threshold_id = threshold_id,
    events = as.integer(nrow(x)),
    stations = as.integer(uniqueN(x$station_number)),
    shape = shape0,

    linear_change_pct_decade = 100 * (exp(beta) - 1),
    linear_ci_low = 100 * (exp(beta_low) - 1),
    linear_ci_high = 100 * (exp(beta_high) - 1),
    linear_interval_excludes_zero = beta_low > 0 || beta_high < 0,
    linear_LR_p = as.numeric(pchisq(lr1, df = 1L, lower.tail = FALSE)),
    linear_delta_AIC = as.numeric(m1$AIC - m0$AIC),

    middle_vs_early_pct = 100 * (exp(middle) - 1),
    recent_vs_early_pct = 100 * (exp(recent) - 1),
    period_LR_p = as.numeric(pchisq(lr2, df = 2L, lower.tail = FALSE)),
    period_delta_AIC = as.numeric(m2$AIC - m0$AIC),

    m0_valid = m0$fit_valid,
    m1_valid = m1$fit_valid,
    m2_valid = m2$fit_valid,
    all_models_valid = m0$fit_valid & m1$fit_valid & m2$fit_valid
  )
}

# ------------------------------------------------------------
# Fit 9 region-threshold combinations
# ------------------------------------------------------------

results_list <- vector("list", nrow(design))
models_list <- vector("list", nrow(design))

for (i in seq_len(nrow(design))) {

  d <- design[i]
  x <- pot[state == d$region & threshold_id == d$threshold_id]

  message(
    "\n========================================",
    "\nRegional threshold sensitivity ", i, " of ", nrow(design),
    ": ", d$region, " / ", d$threshold_id,
    "\n========================================"
  )

  if (nrow(x) < 500L) {
    stop(
      "Too few events for ", d$region, " / ", d$threshold_id,
      ": ", nrow(x),
      call. = FALSE
    )
  }

  message("  Fitting M0")
  m0 <- fit_gpd_m0(x$normalized_excess)

  message("  Fitting M1")
  m1 <- fit_gpd_m1(
    excess = x$normalized_excess,
    year_decade = x$year_decade
  )

  message("  Fitting M2")
  m2 <- fit_gpd_m2(
    excess = x$normalized_excess,
    climate_period = x$climate_period
  )

  results_list[[i]] <- extract_result(
    region = d$region,
    threshold_id = d$threshold_id,
    x = x,
    m0 = m0,
    m1 = m1,
    m2 = m2
  )

  model_key <- paste(d$region, d$threshold_id, sep = "_")
  models_list[[i]] <- list(M0 = m0, M1 = m1, M2 = m2)
  names(models_list)[i] <- model_key
}

results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

results[, region_order := match(region, regions)]
results[, threshold_order := match(threshold_id, thresholds)]
setorder(results, region_order, threshold_order)
results[, c("region_order", "threshold_order") := NULL]

# ------------------------------------------------------------
# Robustness flags by region
# ------------------------------------------------------------

results[, linear_direction := fcase(
  linear_ci_low > 0, "increase",
  linear_ci_high < 0, "decrease",
  default = "no clear trend"
)]

results[, period_model_preferred :=
  period_delta_AIC < -2 & period_LR_p < 0.05
]

results[, middle_lower := middle_vs_early_pct < 0]
results[, recent_above_early := recent_vs_early_pct > 0]

region_robustness <- results[, .(
  specifications = .N,
  all_models_valid = all(all_models_valid),
  linear_increases = sum(linear_direction == "increase"),
  linear_decreases = sum(linear_direction == "decrease"),
  linear_unclear = sum(linear_direction == "no clear trend"),
  period_models_preferred = sum(period_model_preferred),
  middle_lower_count = sum(middle_lower),
  recent_above_early_count = sum(recent_above_early),
  min_linear_change_pct_decade = min(linear_change_pct_decade),
  max_linear_change_pct_decade = max(linear_change_pct_decade),
  min_shape = min(shape),
  max_shape = max(shape)
), by = region]

# ------------------------------------------------------------
# Compact publication-style output
# ------------------------------------------------------------

compact <- results[, .(
  region,
  threshold_id,
  events,
  stations,
  shape,
  linear_change_pct_decade,
  linear_ci_low,
  linear_ci_high,
  linear_LR_p,
  linear_delta_AIC,
  middle_vs_early_pct,
  recent_vs_early_pct,
  period_LR_p,
  period_delta_AIC,
  linear_direction,
  period_model_preferred,
  all_models_valid
)]

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  list(
    design = design,
    sample_audit = sample_audit,
    results = results,
    models = models_list,
    region_robustness = region_robustness
  ),
  file.path(
    model_output_dir,
    "regional_gpd_threshold_sensitivity_models.rds"
  ),
  compress = "xz"
)

fwrite(
  sample_audit,
  file.path(
    table_output_dir,
    "regional_gpd_threshold_sensitivity_sample_audit.csv"
  )
)

fwrite(
  results,
  file.path(
    table_output_dir,
    "regional_gpd_threshold_sensitivity_results.csv"
  )
)

fwrite(
  compact,
  file.path(
    table_output_dir,
    "regional_gpd_threshold_sensitivity_compact.csv"
  )
)

fwrite(
  region_robustness,
  file.path(
    table_output_dir,
    "regional_gpd_threshold_sensitivity_summary.csv"
  )
)

message(
  "\n========================================",
  "\nREGIONAL GPD THRESHOLD SENSITIVITY COMPLETE",
  "\n========================================"
)

print(compact)
print(region_robustness)
