library(targets)

# Keep this list lean until source formats are confirmed.
tar_option_set(
  packages = c("data.table", "dplyr", "yaml", "arrow", "ggplot2")
)

invisible(lapply(list.files("R", pattern = "[.]R$", recursive = TRUE, full.names = TRUE), source))

read_project_config <- function() {
  path <- if (file.exists("config.local.yml")) "config.local.yml" else "config.yml"
  yaml::read_yaml(path)
}

list(
  tar_target(config, read_project_config()),

  tar_target(
    dataset_registry,
    data.table::fread("metadata/dataset_registry.csv")
  ),

  tar_target(
    variable_dictionary,
    data.table::fread("metadata/variable_dictionary.csv")
  ),

  # These targets intentionally stop short of source ingestion until actual
  # Bureau file layouts are supplied and source-specific parsers are tested.
  tar_target(
    project_schema_check,
    {
      stopifnot(nrow(dataset_registry) > 0L)
      stopifnot(nrow(variable_dictionary) > 0L)
      TRUE
    }
  )
)
