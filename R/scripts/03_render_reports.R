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
