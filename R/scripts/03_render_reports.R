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
