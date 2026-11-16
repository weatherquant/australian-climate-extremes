# ============================================================
# 10_rainfall_spatial_duration_robustness.R
# Quasi-Poisson robustness check for temporal changes in
# additional duration of coherent P97.5 / 200-km rainfall episodes.
# ============================================================

library(data.table)
library(arrow)

episode_file <- file.path(
  "data-derived", "rainfall", "spatial_extremes",
  "rainfall_spatial_extreme_episodes_p975_200km_1950_2025.parquet"
)

table_output_dir <- file.path("tables", "rainfall_spatial_extremes")
model_output_file <- file.path(
  "data-derived", "rainfall", "spatial_extremes",
  "rainfall_spatial_duration_robustness_models.rds"
)

dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(episode_file)) {
  stop("Missing episode file: ", episode_file, call. = FALSE)
}

episodes <- as.data.table(arrow::read_parquet(episode_file))

centre_year <- mean(c(1950, 2022))

episodes[, duration_extra_days := as.integer(duration_days - 1L)]
episodes[, year_decade := (start_year - centre_year) / 10]
episodes[, primary_analysis_period := start_year >= 1950L & start_year <= 2022L]
episodes[, sensitivity_analysis_period := start_year >= 1950L & start_year <= 2023L]

episodes[, spatial_climate_period := data.table::fcase(
  start_year >= 1950L & start_year <= 1975L, "1950-1975",
  start_year >= 1976L & start_year <= 2000L, "1976-2000",
  start_year >= 2001L & start_year <= 2022L, "2001-2022",
  default = NA_character_
)]

episodes[, spatial_climate_period := factor(
  spatial_climate_period,
  levels = c("1950-1975", "1976-2000", "2001-2022")
)]

primary <- episodes[
  primary_analysis_period == TRUE &
    all_days_primary_coverage == TRUE &
    is.finite(duration_extra_days) &
    duration_extra_days >= 0L
]

sensitivity <- episodes[
  sensitivity_analysis_period == TRUE &
    all_days_primary_coverage == TRUE &
    is.finite(duration_extra_days) &
    duration_extra_days >= 0L
]

# Primary quasi-Poisson models.
qp_m0 <- glm(
  duration_extra_days ~ 1,
  family = quasipoisson(link = "log"),
  data = primary
)

qp_m1 <- glm(
  duration_extra_days ~ year_decade,
  family = quasipoisson(link = "log"),
  data = primary
)

qp_m2 <- glm(
  duration_extra_days ~ spatial_climate_period,
  family = quasipoisson(link = "log"),
  data = primary
)

# Ordinary Poisson M1 retained only to document dispersion and compare SE.
pois_m1 <- glm(
  duration_extra_days ~ year_decade,
  family = poisson(link = "log"),
  data = primary
)

pearson_dispersion <- sum(
  residuals(pois_m1, type = "pearson")^2
) / df.residual(pois_m1)

extract_linear <- function(fit, period) {
  beta <- coef(fit)["year_decade"]
  se <- sqrt(vcov(fit)["year_decade", "year_decade"])
  z <- qnorm(0.975)

  data.table(
    model = "quasi_poisson_linear",
    analysis_period = period,
    log_rate_change_per_decade = as.numeric(beta),
    standard_error = as.numeric(se),
    rate_ratio_per_decade = as.numeric(exp(beta)),
    percent_change_per_decade = as.numeric(100 * (exp(beta) - 1)),
    conf_low_percent = as.numeric(100 * (exp(beta - z * se) - 1)),
    conf_high_percent = as.numeric(100 * (exp(beta + z * se) - 1)),
    p_value = as.numeric(summary(fit)$coefficients["year_decade", "Pr(>|t|)"]),
    dispersion = as.numeric(summary(fit)$dispersion)
  )
}

qp_linear_effect <- extract_linear(qp_m1, "1950-2022")

period_terms <- c(
  "spatial_climate_period1976-2000",
  "spatial_climate_period2001-2022"
)
period_beta <- coef(qp_m2)[period_terms]
period_se <- sqrt(diag(vcov(qp_m2)))[period_terms]
z <- qnorm(0.975)

qp_period_effects <- data.table(
  comparison = c(
    "1976-2000 vs 1950-1975",
    "2001-2022 vs 1950-1975"
  ),
  log_rate_difference = as.numeric(period_beta),
  standard_error = as.numeric(period_se),
  rate_ratio = as.numeric(exp(period_beta)),
  percent_difference = as.numeric(100 * (exp(period_beta) - 1)),
  conf_low_percent = as.numeric(100 * (exp(period_beta - z * period_se) - 1)),
  conf_high_percent = as.numeric(100 * (exp(period_beta + z * period_se) - 1)),
  p_value = as.numeric(
    summary(qp_m2)$coefficients[period_terms, "Pr(>|t|)"]
  ),
  dispersion = as.numeric(summary(qp_m2)$dispersion)
)

# Quasi-likelihood has no valid likelihood/AIC. Compare nested temporal
# structures with analysis-of-deviance F tests instead.
qp_linear_test <- anova(qp_m0, qp_m1, test = "F")
qp_period_test <- anova(qp_m0, qp_m2, test = "F")

qp_model_tests <- data.table(
  comparison = c("M0 vs M1 linear", "M0 vs M2 three-period"),
  df = c(
    qp_linear_test$Df[2],
    qp_period_test$Df[2]
  ),
  deviance_difference = c(
    qp_linear_test$Deviance[2],
    qp_period_test$Deviance[2]
  ),
  F_statistic = c(
    qp_linear_test$F[2],
    qp_period_test$F[2]
  ),
  p_value = c(
    qp_linear_test$`Pr(>F)`[2],
    qp_period_test$`Pr(>F)`[2]
  )
)

# 1950-2023 sensitivity.
qp_m1_2023 <- glm(
  duration_extra_days ~ year_decade,
  family = quasipoisson(link = "log"),
  data = sensitivity
)

qp_linear_sensitivity_2023 <- extract_linear(
  qp_m1_2023,
  "1950-2023_sensitivity"
)

# Direct comparison of Poisson and quasi-Poisson uncertainty.
poisson_se <- sqrt(vcov(pois_m1)["year_decade", "year_decade"])
quasi_se <- sqrt(vcov(qp_m1)["year_decade", "year_decade"])

uncertainty_comparison <- data.table(
  metric = c(
    "Pearson dispersion",
    "Poisson SE",
    "Quasi-Poisson SE",
    "SE inflation ratio"
  ),
  value = c(
    pearson_dispersion,
    poisson_se,
    quasi_se,
    quasi_se / poisson_se
  )
)

fwrite(
  qp_linear_effect,
  file.path(table_output_dir, "rainfall_spatial_duration_quasipoisson_linear.csv")
)

fwrite(
  qp_period_effects,
  file.path(table_output_dir, "rainfall_spatial_duration_quasipoisson_periods.csv")
)

fwrite(
  qp_model_tests,
  file.path(table_output_dir, "rainfall_spatial_duration_quasipoisson_tests.csv")
)

fwrite(
  qp_linear_sensitivity_2023,
  file.path(table_output_dir, "rainfall_spatial_duration_quasipoisson_sensitivity_2023.csv")
)

fwrite(
  uncertainty_comparison,
  file.path(table_output_dir, "rainfall_spatial_duration_uncertainty_comparison.csv")
)

saveRDS(
  list(
    quasi_M0 = qp_m0,
    quasi_M1 = qp_m1,
    quasi_M2 = qp_m2,
    quasi_M1_2023 = qp_m1_2023,
    poisson_M1 = pois_m1,
    linear_effect = qp_linear_effect,
    period_effects = qp_period_effects,
    model_tests = qp_model_tests,
    sensitivity_2023 = qp_linear_sensitivity_2023,
    uncertainty_comparison = uncertainty_comparison
  ),
  model_output_file,
  compress = "xz"
)

message("\n========================================")
message("QUASI-POISSON DURATION ROBUSTNESS COMPLETE")
message("========================================\n")

message("Primary linear effect:")
print(qp_linear_effect)

message("\nPeriod effects:")
print(qp_period_effects)

message("\nQuasi-Poisson temporal-structure tests:")
print(qp_model_tests)

message("\n1950-2023 sensitivity:")
print(qp_linear_sensitivity_2023)

message("\nPoisson vs quasi-Poisson uncertainty:")
print(uncertainty_comparison)
