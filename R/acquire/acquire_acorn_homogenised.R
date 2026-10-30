# ============================================================
# acquire_acorn_homogenised.R
#
# Acquisition of published homogenised ACORN-SAT daily
# Tmax and Tmin products.
#
# Initial analysis uses ACORN-SAT v2.3.0 to match the
# v2.3 raw-data and supporting-information bundle.
# ============================================================


ACORN_HOMOGENISED_VERSION <- "2.3.0"

ACORN_HOMOGENISED_FILES <- c(tmax = "acorn_sat_v2.3.0_daily_tmax.tar.gz", tmin = "acorn_sat_v2.3.0_daily_tmin.tar.gz")


acquire_acorn_homogenised <- function(cfg) {
  destination_dir <- file.path(cfg$paths$raw_root, cfg$paths$acorn$homogenised)
  
  dir.create(destination_dir,
             recursive = TRUE,
             showWarnings = FALSE)
  
  results <- lapply(names(ACORN_HOMOGENISED_FILES), function(variable) {
    filename <-
      ACORN_HOMOGENISED_FILES[[variable]]
    
    url <- paste0(ACORN_BASE_URL, filename)
    
    destination <- file.path(destination_dir, filename)
    
    message("Acquiring homogenised ACORN ", variable, ": ", filename)
    
    result <- download_file_checked(url = url,
                                    dest = destination,
                                    overwrite = FALSE)
    
    result$variable <- variable
    result$source_version <-
      ACORN_HOMOGENISED_VERSION
    
    result
  })
  
  data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
}


inventory_acorn_tar <- function(tar_file) {
  if (!file.exists(tar_file)) {
    stop("ACORN archive does not exist: ", tar_file, call. = FALSE)
  }
  
  members <- utils::untar(tar_file, list = TRUE)
  
  data.table::data.table(member = members)[, filename := basename(member)][, extension :=
                                                                             tolower(tools::file_ext(filename))]
}