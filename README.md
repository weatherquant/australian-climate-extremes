# Australian Climate Extremes

Reproducible R-first research project for analysing changes in Australian climate extremes, including explicit comparison of raw and homogenised temperature observations, marginal non-stationarity, persistence, compound events and spatial dependence.

## Principles

1. Raw source data are immutable.
2. Raw and homogenised temperature are parallel analysis streams, never overwritten.
3. Large/restricted datasets are stored outside Git and located through `config.local.yml`.
4. Derived analytical datasets use open column schemas and preferably Parquet.
5. The first outputs are data inventory and raw-vs-homogenised audits, not climate trend claims.
6. Every analytical result should be reproducible through the `targets` pipeline.

## Quick start

1. Install R and open this directory as an R project/workspace.
2. Copy `config.example.yml` to `config.local.yml` and edit local data paths.
3. In R, install `renv`, then run `renv::restore()` after a lockfile exists, or `renv::init()` for the initial environment.
4. Add source data to the configured external folders. Do not commit restricted/raw data.
5. Run `targets::tar_make()`.
6. Render/examine the reports in `analysis/`.

## Initial research sequence

- Data availability and station-history audit
- Raw versus homogenised temperature audit
- Temperature marginal/extreme-value analysis
- Persistence and spatial coherence
- Compound/dependence modelling
- Soil-moisture and teleconnection covariates

See `metadata/variable_dictionary.csv` and `metadata/dataset_registry.csv` for the canonical Tier-1 schemas.
