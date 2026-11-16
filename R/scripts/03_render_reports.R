### RUN REPORTS

###################################
#### RUN REPORT - RAW VS HOMOGENISED
quarto::quarto_render(
  "analysis/03_raw_vs_homogenised.qmd"
)

# ------------------------------------------------------------
# Rainfall observing network
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/04_rainfall_network.qmd"
)

# ------------------------------------------------------------
# Extreme rainfall
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/05_extreme_rainfall.qmd"
)

# ------------------------------------------------------------
# Rainfall POT diagnostics
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/06_rainfall_pot_diagnostics.qmd"
)

# ------------------------------------------------------------
# Rainfall GDP diagnostics
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/07_rainfall_gpd_results.qmd"
)

# ------------------------------------------------------------
# Rainfall GDP dependence sensitivity
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/08_rainfall_gpd_dependence_sensitivity.qmd"
)

# ------------------------------------------------------------
# Extreme rainfall concurrence
# ------------------------------------------------------------
quarto::quarto_render(
  "analysis/09_rainfall_spatial_concurrence.qmd"
)

# High-level story
# The evidence does not support one simplistic national narrative.
# Instead:
# - rainfall extremes exhibit strong regional heterogeneity;
# - conditional extreme magnitude has no robust monotonic national increase;
# - NSW shows a robust decline in conditional tail magnitude;
# - Victoria has a positive directional magnitude signal, but dependence-aware uncertainty weakens the inference;
# - Queensland shows more multi-decadal than monotonic behaviour;
# - coherent spatial episodes are not detectably increasing in footprint;
# - episode frequency is approximately stationary to modestly declining;
# - the clearest spatial-event change is shorter persistence;
# - dependence and network exposure materially affect uncertainty and therefore need explicit treatment.
