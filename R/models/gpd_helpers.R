# ============================================================
# Configuration
# ============================================================

analysis_start_year <- 1950L
analysis_end_year <- 2025L
analysis_years <- analysis_end_year - analysis_start_year + 1L

centre_year <- mean(
  c(
    analysis_start_year,
    analysis_end_year
  )
)

primary_threshold_id <- "p975"
primary_run_length <- 3L

accepted_rainfall_quality <- c(
  0L,
  1L,
  2L,
  3L,
  4L
)

representative_years <- c(
  1960L,
  1990L,
  2020L
)

return_periods_years <- c(
  10,
  20,
  50,
  100
)

pot_file <- file.path(
  "data-derived",
  "rainfall",
  "pot",
  "rainfall_pot_events_p975_run3_1950_2025.parquet"
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
# Climate-period definition
# ============================================================

assign_climate_period <- function(year) {
  
  data.table::fcase(
    
    year >= 1950L & year <= 1975L,
    "1950-1975",
    
    year >= 1976L & year <= 2000L,
    "1976-2000",
    
    year >= 2001L & year <= 2025L,
    "2001-2025",
    
    default = NA_character_
  )
}


# ============================================================
# Parameter transformations
# ============================================================

shape_from_raw <- function(raw_shape) {
  0.75 * tanh(raw_shape)
}

raw_from_shape <- function(shape) {
  
  if (
    !is.finite(shape) ||
    abs(shape) >= 0.75
  ) {
    stop(
      "Starting shape must be finite and strictly between -0.75 and 0.75.",
      call. = FALSE
    )
  }
  
  atanh(shape / 0.75)
}

shape_derivative <- function(raw_shape) {
  0.75 * (1 - tanh(raw_shape)^2)
}


# ============================================================
# GPD negative log likelihood with observation-specific scale
# ============================================================

gpd_nll_variable_scale <- function(
    excess,
    log_scale,
    shape) {
  
  excess <- as.numeric(excess)
  log_scale <- as.numeric(log_scale)
  shape <- as.numeric(shape)
  
  if (
    length(excess) == 0L ||
    length(log_scale) != length(excess) ||
    any(!is.finite(excess)) ||
    any(excess <= 0) ||
    any(!is.finite(log_scale)) ||
    !is.finite(shape)
  ) {
    return(1e30)
  }
  
  scale <- exp(log_scale)
  
  if (
    any(!is.finite(scale)) ||
    any(scale <= 0)
  ) {
    return(1e30)
  }
  
  # Exponential limit.
  if (abs(shape) < 1e-7) {
    return(
      sum(
        log_scale +
          excess / scale
      )
    )
  }
  
  support <- 1 + shape * excess / scale
  
  if (
    any(!is.finite(support)) ||
    any(support <= 0)
  ) {
    return(1e30)
  }
  
  sum(
    log_scale +
      (1 / shape + 1) * log(support)
  )
}


# ============================================================
# Numerical gradient
# ============================================================

numerical_gradient <- function(
    par,
    fn,
    rel_step = 1e-5) {
  
  par <- as.numeric(par)
  gradient <- numeric(length(par))
  
  for (i in seq_along(par)) {
    
    step <- rel_step * max(1, abs(par[i]))
    
    plus <- par
    minus <- par
    
    plus[i] <- plus[i] + step
    minus[i] <- minus[i] - step
    
    gradient[i] <-
      (
        fn(plus) -
          fn(minus)
      ) /
      (2 * step)
  }
  
  gradient
}


# ============================================================
# Generic lightweight GPD optimiser
# ============================================================

fit_gpd_model <- function(
    objective,
    starting_vectors,
    parameter_names,
    model_name,
    n_observations,
    confidence_level = 0.95) {
  
  candidate_fits <- lapply(
    starting_vectors,
    function(start) {
      
      tryCatch(
        
        stats::optim(
          par = start,
          fn = objective,
          method = "BFGS",
          hessian = TRUE,
          control = list(
            maxit = 3000L,
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
      "All optimisation attempts failed for model: ",
      model_name,
      call. = FALSE
    )
  }
  
  objective_values <- vapply(
    candidate_fits,
    function(x) {
      as.numeric(x$value)
    },
    numeric(1)
  )
  
  fit <- candidate_fits[[
    which.min(objective_values)
  ]]
  
  names(fit$par) <- parameter_names
  
  covariance_raw <- tryCatch(
    solve(fit$hessian),
    error = function(e) {
      matrix(
        NA_real_,
        nrow = length(fit$par),
        ncol = length(fit$par)
      )
    }
  )
  
  rownames(covariance_raw) <- parameter_names
  colnames(covariance_raw) <- parameter_names
  
  standard_errors_raw <- sqrt(
    pmax(
      diag(covariance_raw),
      0
    )
  )
  
  names(standard_errors_raw) <- parameter_names
  
  gradient <- tryCatch(
    numerical_gradient(
      par = fit$par,
      fn = objective
    ),
    error = function(e) {
      rep(NA_real_, length(fit$par))
    }
  )
  
  names(gradient) <- parameter_names
  
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
    all(is.finite(hessian_eigenvalues)) &&
    all(hessian_eigenvalues > 0)
  
  maximum_gradient <- if (any(is.finite(gradient))) {
    max(abs(gradient), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  parameter_count <- length(fit$par)
  
  list(
    model_name = model_name,
    n = as.integer(n_observations),
    parameters_raw = fit$par,
    covariance_raw = covariance_raw,
    standard_errors_raw = standard_errors_raw,
    gradient = gradient,
    maximum_gradient = as.numeric(maximum_gradient),
    negative_log_likelihood = as.numeric(fit$value),
    log_likelihood = as.numeric(-fit$value),
    parameter_count = as.integer(parameter_count),
    AIC = as.numeric(2 * parameter_count + 2 * fit$value),
    convergence = as.integer(fit$convergence),
    convergence_message = if (length(fit$message) == 0L) {
      ""
    } else {
      as.character(fit$message)
    },
    hessian_positive_definite = hessian_positive_definite,
    fit_valid =
      fit$convergence == 0L &&
      hessian_positive_definite &&
      all(
        is.finite(
          fit$par
        )
      ),
    starting_values_tested = length(starting_vectors),
    confidence_level = confidence_level
  )
}


# ============================================================
# Fit M0: stationary GPD
# ============================================================

fit_gpd_m0 <- function(
    excess,
    starting_shapes = c(-0.10, 0.00, 0.10)) {
  
  excess <- as.numeric(
    excess[
      is.finite(excess) &
        excess > 0
    ]
  )
  
  starting_scale <- max(mean(excess), 1e-4)
  
  objective <- function(theta) {
    
    log_scale <- theta[1L]
    shape <- shape_from_raw(theta[2L])
    
    gpd_nll_variable_scale(
      excess = excess,
      log_scale = rep(log_scale, length(excess)),
      shape = shape
    )
  }
  
  starting_vectors <- lapply(
    starting_shapes,
    function(start_shape) {
      c(
        log_scale_intercept = log(starting_scale),
        raw_shape = raw_from_shape(start_shape)
      )
    }
  )
  
  fit_gpd_model(
    objective = objective,
    starting_vectors = starting_vectors,
    parameter_names = c(
      "log_scale_intercept",
      "raw_shape"
    ),
    model_name = "M0_stationary",
    n_observations = length(excess)
  )
}


# ============================================================
# Fit M1: linear trend in GPD scale
# ============================================================

fit_gpd_m1 <- function(
    excess,
    year_decade,
    starting_shapes = c(-0.10, 0.00, 0.10),
    starting_time_slopes = c(-0.02, 0.00, 0.02)) {
  
  valid <-
    is.finite(excess) &
    excess > 0 &
    is.finite(year_decade)
  
  excess <- as.numeric(excess[valid])
  year_decade <- as.numeric(year_decade[valid])
  
  starting_scale <- max(mean(excess), 1e-4)
  
  objective <- function(theta) {
    
    log_scale <-
      theta[1L] +
      theta[2L] * year_decade
    
    shape <- shape_from_raw(theta[3L])
    
    gpd_nll_variable_scale(
      excess = excess,
      log_scale = log_scale,
      shape = shape
    )
  }
  
  starting_grid <- data.table::CJ(
    start_shape = starting_shapes,
    start_time = starting_time_slopes
  )
  
  starting_vectors <- lapply(
    seq_len(nrow(starting_grid)),
    function(i) {
      c(
        log_scale_intercept = log(starting_scale),
        scale_time_per_decade = starting_grid$start_time[i],
        raw_shape = raw_from_shape(starting_grid$start_shape[i])
      )
    }
  )
  
  fit_gpd_model(
    objective = objective,
    starting_vectors = starting_vectors,
    parameter_names = c(
      "log_scale_intercept",
      "scale_time_per_decade",
      "raw_shape"
    ),
    model_name = "M1_linear_scale_trend",
    n_observations = length(excess)
  )
}


# ============================================================
# Fit M2: three climate-period effects in GPD scale
# ============================================================

fit_gpd_m2 <- function(
    excess,
    climate_period,
    starting_shapes = c(-0.10, 0.00, 0.10)) {
  
  period_levels <- c(
    "1950-1975",
    "1976-2000",
    "2001-2025"
  )
  
  climate_period <- factor(
    climate_period,
    levels = period_levels
  )
  
  valid <-
    is.finite(excess) &
    excess > 0 &
    !is.na(climate_period)
  
  excess <- as.numeric(excess[valid])
  climate_period <- droplevels(climate_period[valid])
  
  period_middle <- as.numeric(
    climate_period == "1976-2000"
  )
  
  period_recent <- as.numeric(
    climate_period == "2001-2025"
  )
  
  starting_scale <- max(mean(excess), 1e-4)
  
  objective <- function(theta) {
    
    log_scale <-
      theta[1L] +
      theta[2L] * period_middle +
      theta[3L] * period_recent
    
    shape <- shape_from_raw(theta[4L])
    
    gpd_nll_variable_scale(
      excess = excess,
      log_scale = log_scale,
      shape = shape
    )
  }
  
  starting_vectors <- lapply(
    starting_shapes,
    function(start_shape) {
      c(
        log_scale_intercept = log(starting_scale),
        period_1976_2000 = 0,
        period_2001_2025 = 0,
        raw_shape = raw_from_shape(start_shape)
      )
    }
  )
  
  fit_gpd_model(
    objective = objective,
    starting_vectors = starting_vectors,
    parameter_names = c(
      "log_scale_intercept",
      "period_1976_2000",
      "period_2001_2025",
      "raw_shape"
    ),
    model_name = "M2_three_period_scale",
    n_observations = length(excess)
  )
}


# ============================================================
# Extract interpretable model parameters
# ============================================================

extract_m0_parameters <- function(fit) {
  
  par <- fit$parameters_raw
  covariance <- fit$covariance_raw
  
  log_scale <- par["log_scale_intercept"]
  raw_shape <- par["raw_shape"]
  
  scale <- exp(log_scale)
  shape <- shape_from_raw(raw_shape)
  
  se_log_scale <- sqrt(max(covariance["log_scale_intercept", "log_scale_intercept"], 0))
  se_raw_shape <- sqrt(max(covariance["raw_shape", "raw_shape"], 0))
  
  z <- qnorm(0.975)
  
  data.table(
    model = fit$model_name,
    term = c("scale_ratio", "shape"),
    estimate = c(scale, shape),
    std_error = c(
      scale * se_log_scale,
      shape_derivative(raw_shape) * se_raw_shape
    ),
    conf_low = c(
      exp(log_scale - z * se_log_scale),
      shape_from_raw(raw_shape - z * se_raw_shape)
    ),
    conf_high = c(
      exp(log_scale + z * se_log_scale),
      shape_from_raw(raw_shape + z * se_raw_shape)
    )
  )
}


extract_m1_parameters <- function(fit) {
  
  par <- fit$parameters_raw
  covariance <- fit$covariance_raw
  z <- qnorm(0.975)
  
  log_scale <- par["log_scale_intercept"]
  beta_time <- par["scale_time_per_decade"]
  raw_shape <- par["raw_shape"]
  
  se_log_scale <- sqrt(max(covariance["log_scale_intercept", "log_scale_intercept"], 0))
  se_time <- sqrt(max(covariance["scale_time_per_decade", "scale_time_per_decade"], 0))
  se_raw_shape <- sqrt(max(covariance["raw_shape", "raw_shape"], 0))
  
  data.table(
    model = fit$model_name,
    term = c(
      "scale_ratio_at_centre_year",
      "scale_time_per_decade",
      "shape"
    ),
    estimate = c(
      exp(log_scale),
      beta_time,
      shape_from_raw(raw_shape)
    ),
    std_error = c(
      exp(log_scale) * se_log_scale,
      se_time,
      shape_derivative(raw_shape) * se_raw_shape
    ),
    conf_low = c(
      exp(log_scale - z * se_log_scale),
      beta_time - z * se_time,
      shape_from_raw(raw_shape - z * se_raw_shape)
    ),
    conf_high = c(
      exp(log_scale + z * se_log_scale),
      beta_time + z * se_time,
      shape_from_raw(raw_shape + z * se_raw_shape)
    )
  )
}


extract_m2_parameters <- function(fit) {
  
  par <- fit$parameters_raw
  covariance <- fit$covariance_raw
  z <- qnorm(0.975)
  
  log_scale_early <- par["log_scale_intercept"]
  log_scale_middle <- log_scale_early + par["period_1976_2000"]
  log_scale_recent <- log_scale_early + par["period_2001_2025"]
  
  contrast_variance <- function(weights) {
    as.numeric(
      t(weights) %*%
        covariance %*%
        weights
    )
  }
  
  weights_early <- c(1, 0, 0, 0)
  weights_middle <- c(1, 1, 0, 0)
  weights_recent <- c(1, 0, 1, 0)
  
  se_early <- sqrt(max(contrast_variance(weights_early), 0))
  se_middle <- sqrt(max(contrast_variance(weights_middle), 0))
  se_recent <- sqrt(max(contrast_variance(weights_recent), 0))
  
  raw_shape <- par["raw_shape"]
  se_raw_shape <- sqrt(max(covariance["raw_shape", "raw_shape"], 0))
  
  data.table(
    model = fit$model_name,
    term = c(
      "scale_ratio_1950_1975",
      "scale_ratio_1976_2000",
      "scale_ratio_2001_2025",
      "shape"
    ),
    estimate = c(
      exp(log_scale_early),
      exp(log_scale_middle),
      exp(log_scale_recent),
      shape_from_raw(raw_shape)
    ),
    std_error = c(
      exp(log_scale_early) * se_early,
      exp(log_scale_middle) * se_middle,
      exp(log_scale_recent) * se_recent,
      shape_derivative(raw_shape) * se_raw_shape
    ),
    conf_low = c(
      exp(log_scale_early - z * se_early),
      exp(log_scale_middle - z * se_middle),
      exp(log_scale_recent - z * se_recent),
      shape_from_raw(raw_shape - z * se_raw_shape)
    ),
    conf_high = c(
      exp(log_scale_early + z * se_early),
      exp(log_scale_middle + z * se_middle),
      exp(log_scale_recent + z * se_recent),
      shape_from_raw(raw_shape + z * se_raw_shape)
    )
  )
}


# ============================================================
# GPD return-level function
# ============================================================

pot_return_level <- function(
    threshold_mm,
    scale_ratio,
    shape,
    event_rate_per_year,
    return_period_years) {
  
  scale_mm <- scale_ratio * threshold_mm
  expected_exceedances <- event_rate_per_year * return_period_years
  
  result <- rep(NA_real_, length(expected_exceedances))
  
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
      scale_mm[valid] * log(expected_exceedances[valid])
    
  } else {
    
    result[valid] <-
      threshold_mm[valid] +
      scale_mm[valid] / shape *
      (
        expected_exceedances[valid]^shape -
          1
      )
  }
  
  result
}
