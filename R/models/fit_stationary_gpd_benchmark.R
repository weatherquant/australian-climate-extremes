# ============================================================
# fit_stationary_gpd_benchmark.R
#
# Fast stationary GPD benchmark for Australian extreme rainfall.
#
# Primary specification:
#   Threshold:         station-specific P97.5
#   Declustering:      3-day runs
#   Quality rule:      PRCP_QUAL 0-4
#   Analysis period:   1950-2025
#
# Sensitivities:
#   Thresholds:        P95, P97.5, P99
#   Run lengths:       2, 3, 5 days
#
# Computational design:
#   - no spatial random field
#   - no station-by-station GPD optimisation
#   - one national pooled fit per threshold/run combination
#   - optional pooled state fits
#
# Excesses are divided by the station-specific threshold before
# fitting. This avoids pooling raw millimetre excesses from very
# different Australian rainfall climates.
#
# The fitted GPD parameters are:
#   scale_ratio = GPD scale / station threshold
#   shape       = common pooled tail-shape parameter
#
# Station-specific millimetre scale is recovered as:
#   scale_mm = scale_ratio * threshold_mm
# ============================================================


# ============================================================
# Configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L
analysis_years <- analysis_end_year - analysis_start_year + 1L

primary_threshold_id <- "p975"
primary_run_length <- 3L

accepted_rainfall_quality <- c(
  0L,
  1L,
  2L,
  3L,
  4L
)

return_periods_years <- c(
  2,
  5,
  10,
  20,
  50,
  100
)

pot_file <- file.path(
  "data-derived",
  "rainfall",
  "pot",
  "rainfall_pot_events_all_thresholds_runs_1950_2025.parquet"
)

model_output_dir <- file.path(
  "data-derived",
  "rainfall",
  "models",
  "gpd"
)

table_output_dir <- file.path(
  "tables",
  "rainfall_models",
  "gpd"
)


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "data.table",
  "arrow"
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

library(data.table)
library(arrow)


dir.create(
  model_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# GPD likelihood
# ============================================================

gpd_negative_log_likelihood <- function(
    parameters,
    excess) {
  
  log_scale <- parameters[1L]
  shape <- parameters[2L]
  
  scale <- exp(log_scale)
  
  if (
    !is.finite(scale) ||
    scale <= 0
  ) {
    return(1e30)
  }
  
  excess <- as.numeric(excess)
  
  if (
    length(excess) == 0L ||
    any(!is.finite(excess)) ||
    any(excess <= 0)
  ) {
    return(1e30)
  }
  
  # Exponential limit as shape approaches zero.
  if (abs(shape) < 1e-7) {
    return(
      length(excess) * log(scale) +
        sum(excess) / scale
    )
  }
  
  support <- 1 + shape * excess / scale
  
  if (
    any(!is.finite(support)) ||
    any(support <= 0)
  ) {
    return(1e30)
  }
  
  length(excess) * log(scale) +
    (
      1 / shape + 1
    ) *
    sum(
      log(support)
    )
}


# ============================================================
# Fit one pooled GPD
# ============================================================

fit_stationary_gpd <- function(
    excess,
    confidence_level = 0.95) {
  
  excess <- as.numeric(
    excess[
      is.finite(excess) &
        excess > 0
    ]
  )
  
  n <- length(excess)
  
  if (n < 50L) {
    stop(
      "At least 50 positive excesses are required.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Starting values
  # ----------------------------------------------------------
  
  starting_scale <- max(
    mean(excess),
    1e-4
  )
  
  
  # ----------------------------------------------------------
  # Unconstrained parameterisation
  #
  # theta[1] = log(scale)
  #
  # shape is transformed to lie within approximately
  # (-0.75, 0.75):
  #
  # shape = 0.75 * tanh(theta[2])
  #
  # This avoids the L-BFGS-B behaviour observed previously,
  # where shape remained exactly at its starting value.
  # ----------------------------------------------------------
  
  objective <- function(theta) {
    
    log_scale <- theta[1]
    
    shape <- 0.75 *
      tanh(
        theta[2]
      )
    
    gpd_negative_log_likelihood(
      parameters = c(
        log_scale,
        shape
      ),
      excess = excess
    )
  }
  
  
  # ----------------------------------------------------------
  # Fit from multiple starting shapes
  #
  # Deliberately inexpensive: only three two-parameter fits.
  # ----------------------------------------------------------
  
  starting_shapes <- c(
    -0.10,
    0.00,
    0.10
  )
  
  
  candidate_fits <- lapply(
    starting_shapes,
    function(start_shape) {
      
      raw_shape_start <-
        atanh(
          start_shape /
            0.75
        )
      
      tryCatch(
        
        stats::optim(
          par = c(
            log(starting_scale),
            raw_shape_start
          ),
          fn = objective,
          method = "BFGS",
          hessian = TRUE,
          control = list(
            maxit = 2000L,
            reltol = 1e-10
          )
        ),
        
        error = function(e) {
          NULL
        }
      )
    }
  )
  
  
  candidate_fits <- candidate_fits[
    !vapply(
      candidate_fits,
      is.null,
      logical(1)
    )
  ]
  
  
  if (length(candidate_fits) == 0L) {
    stop(
      "All GPD optimisation attempts failed.",
      call. = FALSE
    )
  }
  
  
  objective_values <- vapply(
    candidate_fits,
    function(x) {
      as.numeric(
        x$value
      )
    },
    numeric(1)
  )
  
  
  fit <- candidate_fits[[
    which.min(
      objective_values
    )
  ]]
  
  
  # ----------------------------------------------------------
  # Back-transform estimates
  # ----------------------------------------------------------
  
  log_scale_hat <-
    fit$par[1]
  
  raw_shape_hat <-
    fit$par[2]
  
  scale_hat <-
    exp(
      log_scale_hat
    )
  
  shape_hat <-
    0.75 *
    tanh(
      raw_shape_hat
    )
  
  
  # ----------------------------------------------------------
  # Numerical gradient at optimum
  # ----------------------------------------------------------
  
  gradient <- tryCatch(
    
    stats::optimHess(
      par = fit$par,
      fn = objective
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Covariance on transformed optimisation scale
  # ----------------------------------------------------------
  
  covariance_raw <- tryCatch(
    
    solve(
      fit$hessian
    ),
    
    error = function(e) {
      matrix(
        NA_real_,
        2L,
        2L
      )
    }
  )
  
  
  # ----------------------------------------------------------
  # Delta-method transformation
  # ----------------------------------------------------------
  
  derivative_scale <-
    scale_hat
  
  derivative_shape <-
    0.75 *
    (
      1 -
        tanh(
          raw_shape_hat
        )^2
    )
  
  
  se_log_scale <-
    sqrt(
      max(
        covariance_raw[1, 1],
        0
      )
    )
  
  se_raw_shape <-
    sqrt(
      max(
        covariance_raw[2, 2],
        0
      )
    )
  
  
  se_scale <-
    derivative_scale *
    se_log_scale
  
  se_shape <-
    derivative_shape *
    se_raw_shape
  
  
  # ----------------------------------------------------------
  # Confidence intervals
  # ----------------------------------------------------------
  
  z <- stats::qnorm(
    1 -
      (
        1 -
          confidence_level
      ) /
      2
  )
  
  
  scale_conf_low <-
    exp(
      log_scale_hat -
        z *
        se_log_scale
    )
  
  
  scale_conf_high <-
    exp(
      log_scale_hat +
        z *
        se_log_scale
    )
  
  
  raw_shape_conf_low <-
    raw_shape_hat -
    z *
    se_raw_shape
  
  
  raw_shape_conf_high <-
    raw_shape_hat +
    z *
    se_raw_shape
  
  
  shape_conf_low <-
    0.75 *
    tanh(
      raw_shape_conf_low
    )
  
  
  shape_conf_high <-
    0.75 *
    tanh(
      raw_shape_conf_high
    )
  
  
  # ----------------------------------------------------------
  # Hessian diagnostic
  # ----------------------------------------------------------
  
  hessian_eigenvalues <- tryCatch(
    
    eigen(
      fit$hessian,
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    
    error = function(e) {
      NA_real_
    }
  )
  
  
  hessian_positive_definite <-
    all(
      is.finite(
        hessian_eigenvalues
      )
    ) &&
    all(
      hessian_eigenvalues > 0
    )
  
  
  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  
  list(
    
    n =
      as.integer(n),
    
    scale =
      as.numeric(
        scale_hat
      ),
    
    shape =
      as.numeric(
        shape_hat
      ),
    
    scale_standard_error =
      as.numeric(
        se_scale
      ),
    
    shape_standard_error =
      as.numeric(
        se_shape
      ),
    
    scale_conf_low =
      as.numeric(
        scale_conf_low
      ),
    
    scale_conf_high =
      as.numeric(
        scale_conf_high
      ),
    
    shape_conf_low =
      as.numeric(
        shape_conf_low
      ),
    
    shape_conf_high =
      as.numeric(
        shape_conf_high
      ),
    
    negative_log_likelihood =
      as.numeric(
        fit$value
      ),
    
    AIC =
      as.numeric(
        4 +
          2 *
          fit$value
      ),
    
    convergence =
      as.integer(
        fit$convergence
      ),
    
    convergence_message =
      as.character(
        fit$message
      ),
    
    maximum_gradient =
      NA_real_,
    
    hessian_positive_definite =
      hessian_positive_definite,
    
    starting_shapes_tested =
      paste(
        starting_shapes,
        collapse = ","
      )
  )
}


# ============================================================
# GPD quantile of threshold excess
# ============================================================

gpd_excess_quantile <- function(
    probability,
    scale,
    shape) {
  
  probability <- as.numeric(probability)
  
  if (
    any(probability <= 0) ||
    any(probability >= 1)
  ) {
    stop(
      "probability must be strictly between zero and one.",
      call. = FALSE
    )
  }
  
  if (abs(shape) < 1e-7) {
    return(
      -scale *
        log(
          1 -
            probability
        )
    )
  }
  
  scale / shape *
    (
      (
        1 -
          probability
      )^(-shape) -
        1
    )
}


# ============================================================
# Return level from Poisson rate plus GPD magnitude
# ============================================================

pot_return_level <- function(
    threshold_mm,
    scale_ratio,
    shape,
    event_rate_per_year,
    return_period_years) {
  
  threshold_mm <- as.numeric(
    threshold_mm
  )
  
  scale_mm <-
    scale_ratio *
    threshold_mm
  
  expected_exceedances <-
    event_rate_per_year *
    return_period_years
  
  result <- rep(
    NA_real_,
    length(
      expected_exceedances
    )
  )
  
  valid <-
    is.finite(expected_exceedances) &
    expected_exceedances > 1 &
    is.finite(threshold_mm) &
    threshold_mm > 0 &
    is.finite(scale_mm) &
    scale_mm > 0
  
  if (abs(shape) < 1e-7) {
    
    result[valid] <-
      threshold_mm[valid] +
      scale_mm[valid] *
      log(
        expected_exceedances[valid]
      )
    
  } else {
    
    result[valid] <-
      threshold_mm[valid] +
      scale_mm[valid] /
      shape *
      (
        expected_exceedances[valid]^shape -
          1
      )
  }
  
  result
}


# ============================================================
# Model diagnostics
# ============================================================

calculate_gpd_diagnostics <- function(
    excess,
    scale,
    shape,
    probabilities = c(
      0.50,
      0.75,
      0.90,
      0.95,
      0.99
    )) {
  
  excess <- sort(
    as.numeric(
      excess[
        is.finite(excess) &
          excess > 0
      ]
    )
  )
  
  empirical <- as.numeric(
    stats::quantile(
      excess,
      probs = probabilities,
      names = FALSE,
      type = 8
    )
  )
  
  fitted <- gpd_excess_quantile(
    probability = probabilities,
    scale = scale,
    shape = shape
  )
  
  data.table(
    probability = probabilities,
    empirical_excess = empirical,
    fitted_excess = fitted,
    difference = fitted - empirical,
    relative_difference =
      (
        fitted -
          empirical
      ) /
      empirical
  )
}


# ============================================================
# Load and prepare POT events
# ============================================================

if (!file.exists(pot_file)) {
  stop(
    "POT event dataset does not exist: ",
    pot_file,
    call. = FALSE
  )
}

pot_events <- as.data.table(
  arrow::read_parquet(
    pot_file
  )
)

required_columns <- c(
  "station_number",
  "station_name",
  "state",
  "latitude",
  "longitude",
  "threshold_id",
  "threshold_mm",
  "run_length_days",
  "event_peak_date",
  "event_peak_mm",
  "excess_mm",
  "rainfall_quality_peak"
)

missing_columns <- setdiff(
  required_columns,
  names(pot_events)
)

if (length(missing_columns) > 0L) {
  stop(
    "POT events are missing columns: ",
    paste(
      missing_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

pot_events <- pot_events[
  rainfall_quality_peak %in%
    accepted_rainfall_quality &
    is.finite(excess_mm) &
    excess_mm > 0 &
    is.finite(threshold_mm) &
    threshold_mm > 0
]

pot_events[
  ,
  normalized_excess :=
    excess_mm /
    threshold_mm
]


# ============================================================
# National threshold and run-length sensitivity fits
# ============================================================

fit_groups <- unique(
  pot_events[
    ,
    .(
      threshold_id,
      run_length_days
    )
  ]
)

data.table::setorder(
  fit_groups,
  threshold_id,
  run_length_days
)

national_fit_list <- vector(
  mode = "list",
  length = nrow(fit_groups)
)

national_model_objects <- vector(
  mode = "list",
  length = nrow(fit_groups)
)

for (
  i in seq_len(
    nrow(fit_groups)
  )
) {
  
  threshold_value <-
    fit_groups$threshold_id[i]
  
  run_value <-
    fit_groups$run_length_days[i]
  
  message(
    "Fitting national stationary GPD: ",
    threshold_value,
    ", run length ",
    run_value
  )
  
  x <- pot_events[
    threshold_id ==
      threshold_value &
      run_length_days ==
      run_value
  ]
  
  fit <- fit_stationary_gpd(
    x$normalized_excess
  )
  
  model_name <- paste0(
    threshold_value,
    "_run",
    run_value
  )
  
  national_model_objects[[i]] <- fit
  names(national_model_objects)[i] <-
    model_name
  
  national_fit_list[[i]] <- data.table(
    
    model_name =
      model_name,
    
    threshold_id =
      threshold_value,
    
    run_length_days =
      as.integer(
        run_value
      ),
    
    events =
      as.integer(
        fit$n
      ),
    
    stations =
      as.integer(
        uniqueN(
          x$station_number
        )
      ),
    
    scale_ratio =
      fit$scale,
    
    shape =
      fit$shape,
    
    scale_standard_error =
      fit$scale_standard_error,
    
    shape_standard_error =
      fit$shape_standard_error,
    
    scale_conf_low =
      fit$scale_conf_low,
    
    scale_conf_high =
      fit$scale_conf_high,
    
    shape_conf_low =
      fit$shape_conf_low,
    
    shape_conf_high =
      fit$shape_conf_high,
    
    negative_log_likelihood =
      fit$negative_log_likelihood,
    
    AIC =
      fit$AIC,
    
    convergence =
      fit$convergence,
    
    maximum_gradient =
      fit$maximum_gradient,
    
    hessian_positive_definite =
      fit$hessian_positive_definite
  )
}

national_gpd_results <- rbindlist(
  national_fit_list,
  use.names = TRUE,
  fill = TRUE
)


# ============================================================
# Primary national benchmark
# ============================================================

primary_events <- pot_events[
  threshold_id ==
    primary_threshold_id &
    run_length_days ==
    primary_run_length
]

primary_fit_name <- paste0(
  primary_threshold_id,
  "_run",
  primary_run_length
)

primary_fit <-
  national_model_objects[[
    primary_fit_name
    ]]

if (is.null(primary_fit)) {
  stop(
    "Primary GPD fit was not found: ",
    primary_fit_name,
    call. = FALSE
  )
}

primary_diagnostics <- calculate_gpd_diagnostics(
  excess =
    primary_events$normalized_excess,
  
  scale =
    primary_fit$scale,
  
  shape =
    primary_fit$shape
)


# ============================================================
# Station event rates and return levels
# ============================================================

station_event_rates <- primary_events[
  ,
  .(
    events =
      as.integer(.N),
    
    event_rate_per_year =
      as.numeric(
        .N /
          analysis_years
      ),
    
    threshold_mm =
      as.numeric(
        data.table::first(
          threshold_mm
        )
      )
  ),
  by = .(
    station_number,
    station_name,
    state,
    latitude,
    longitude
  )
]

station_return_levels <- rbindlist(
  lapply(
    return_periods_years,
    function(return_period) {
      
      station_event_rates[
        ,
        .(
          station_number,
          station_name,
          state,
          latitude,
          longitude,
          events,
          event_rate_per_year,
          threshold_mm,
          
          return_period_years =
            as.numeric(
              return_period
            ),
          
          return_level_mm =
            pot_return_level(
              threshold_mm =
                threshold_mm,
              
              scale_ratio =
                primary_fit$scale,
              
              shape =
                primary_fit$shape,
              
              event_rate_per_year =
                event_rate_per_year,
              
              return_period_years =
                return_period
            )
        )
      ]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)


# ============================================================
# National return-level summaries
# ============================================================

return_level_summary <- station_return_levels[
  is.finite(return_level_mm),
  .(
    stations =
      as.integer(.N),
    
    median_return_level_mm =
      as.numeric(
        median(
          return_level_mm
        )
      ),
    
    q10_return_level_mm =
      as.numeric(
        quantile(
          return_level_mm,
          0.10,
          names = FALSE
        )
      ),
    
    q90_return_level_mm =
      as.numeric(
        quantile(
          return_level_mm,
          0.90,
          names = FALSE
        )
      ),
    
    minimum_return_level_mm =
      as.numeric(
        min(
          return_level_mm
        )
      ),
    
    maximum_return_level_mm =
      as.numeric(
        max(
          return_level_mm
        )
      )
  ),
  by = return_period_years
]


# ============================================================
# Fast pooled state fits
#
# These are descriptive benchmark fits, not final hierarchical
# regional models.
# ============================================================

state_groups <- primary_events[
  ,
  .N,
  by = state
][
  N >= 1000L,
  state
]

state_fit_list <- lapply(
  state_groups,
  function(state_value) {
    
    message(
      "Fitting stationary state benchmark: ",
      state_value
    )
    
    x <- primary_events[
      state ==
        state_value
    ]
    
    fit <- fit_stationary_gpd(
      x$normalized_excess
    )
    
    data.table(
      
      state =
        state_value,
      
      events =
        as.integer(
          fit$n
        ),
      
      stations =
        as.integer(
          uniqueN(
            x$station_number
          )
        ),
      
      scale_ratio =
        fit$scale,
      
      shape =
        fit$shape,
      
      scale_conf_low =
        fit$scale_conf_low,
      
      scale_conf_high =
        fit$scale_conf_high,
      
      shape_conf_low =
        fit$shape_conf_low,
      
      shape_conf_high =
        fit$shape_conf_high,
      
      convergence =
        fit$convergence,
      
      maximum_gradient =
        fit$maximum_gradient,
      
      hessian_positive_definite =
        fit$hessian_positive_definite
    )
  }
)

state_gpd_results <- rbindlist(
  state_fit_list,
  use.names = TRUE,
  fill = TRUE
)


# ============================================================
# Validation
# ============================================================

primary_fit_valid <-
  primary_fit$convergence == 0L &&
  isTRUE(
    primary_fit$hessian_positive_definite
  ) &&
  is.finite(
    primary_fit$scale
  ) &&
  primary_fit$scale > 0 &&
  is.finite(
    primary_fit$shape
  )

if (!primary_fit_valid) {
  warning(
    "Primary stationary GPD did not pass all validation checks.",
    call. = FALSE
  )
}

benchmark_audit <- data.table(
  
  analysis_start_year =
    analysis_start_year,
  
  analysis_end_year =
    analysis_end_year,
  
  analysis_years =
    analysis_years,
  
  primary_threshold_id =
    primary_threshold_id,
  
  primary_run_length =
    primary_run_length,
  
  accepted_quality_codes =
    paste(
      accepted_rainfall_quality,
      collapse = ","
    ),
  
  events =
    nrow(
      primary_events
    ),
  
  stations =
    uniqueN(
      primary_events$station_number
    ),
  
  scale_ratio =
    primary_fit$scale,
  
  shape =
    primary_fit$shape,
  
  convergence =
    primary_fit$convergence,
  
  maximum_gradient =
    primary_fit$maximum_gradient,
  
  hessian_positive_definite =
    primary_fit$hessian_positive_definite,
  
  fit_valid =
    primary_fit_valid
)


# ============================================================
# Save outputs
# ============================================================

saveRDS(
  list(
    national_models =
      national_model_objects,
    
    primary_model_name =
      primary_fit_name,
    
    primary_model =
      primary_fit,
    
    configuration =
      list(
        analysis_start_year =
          analysis_start_year,
        
        analysis_end_year =
          analysis_end_year,
        
        primary_threshold_id =
          primary_threshold_id,
        
        primary_run_length =
          primary_run_length,
        
        accepted_rainfall_quality =
          accepted_rainfall_quality
      )
  ),
  file.path(
    model_output_dir,
    "stationary_gpd_benchmark_models.rds"
  ),
  compress = "xz"
)

fwrite(
  benchmark_audit,
  file.path(
    table_output_dir,
    "stationary_gpd_benchmark_audit.csv"
  )
)

fwrite(
  national_gpd_results,
  file.path(
    table_output_dir,
    "stationary_gpd_threshold_run_sensitivity.csv"
  )
)

fwrite(
  primary_diagnostics,
  file.path(
    table_output_dir,
    "stationary_gpd_primary_quantile_diagnostics.csv"
  )
)

fwrite(
  station_event_rates,
  file.path(
    table_output_dir,
    "stationary_gpd_station_event_rates.csv"
  )
)

fwrite(
  station_return_levels,
  file.path(
    table_output_dir,
    "stationary_gpd_station_return_levels.csv"
  )
)

fwrite(
  return_level_summary,
  file.path(
    table_output_dir,
    "stationary_gpd_return_level_summary.csv"
  )
)

fwrite(
  state_gpd_results,
  file.path(
    table_output_dir,
    "stationary_gpd_state_benchmarks.csv"
  )
)


# ============================================================
# Console summary
# ============================================================

message(
  "\n========================================",
  "\nSTATIONARY GPD BENCHMARK COMPLETE",
  "\n========================================"
)

print(
  benchmark_audit
)

print(
  national_gpd_results[
    order(
      threshold_id,
      run_length_days
    )
  ]
)

print(
  primary_diagnostics
)

print(
  return_level_summary
)

print(
  state_gpd_results[
    order(state)
  ]
)
