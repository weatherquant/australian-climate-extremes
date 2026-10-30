# ============================================================
# acquire_acorn.R
#
# Acquisition and inventory of ACORN-SAT source data
# ============================================================


ACORN_BASE_URL <-
  paste0("ftp://ftp.bom.gov.au/",
         "anon/home/ncc/www/change/",
         "ACORN_SAT_daily/")


ACORN_RAW_BUNDLE <-
  "v2.3-raw-data-info-bundle.zip"


ACORN_VERSION <- "2.3"



acquire_acorn_raw <- function(cfg) {
  raw_root <-
    cfg$paths$raw_root
  
  
  destination_dir <-
    file.path(raw_root, cfg$paths$acorn$raw)
  
  
  dir.create(destination_dir,
             recursive = TRUE,
             showWarnings = FALSE)
  
  
  destination <-
    file.path(destination_dir, ACORN_RAW_BUNDLE)
  
  
  url <-
    paste0(ACORN_BASE_URL, ACORN_RAW_BUNDLE)
  
  
  result <-
    download_file_checked(url = url,
                          dest = destination,
                          overwrite = FALSE)
  
  
  result
}



inventory_acorn_zip <- function(zip_file) {
  if (!file.exists(zip_file)) {
    stop("ACORN archive does not exist: ", zip_file, call. = FALSE)
    
  }
  
  
  inventory <-
    unzip(zip_file, list = TRUE)
  
  
  inventory |>
    dplyr::transmute(
      member = Name,
      
      bytes = Length,
      
      modified = Date,
      
      extension =
        tolower(tools::file_ext(Name))
      
    )
}



extract_acorn_archive <- function(zip_file, destination, overwrite = FALSE) {
  if (!file.exists(zip_file)) {
    stop("ACORN archive does not exist: ", zip_file, call. = FALSE)
    
  }
  
  
  if (dir.exists(destination) &&
      length(list.files(destination)) > 0 &&
      !overwrite) {
    message("ACORN extraction directory already populated: ",
            destination)
    
    return(normalizePath(destination, winslash = "/"))
  }
  
  
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  
  
  unzip(zip_file, exdir = destination)
  
  
  message("ACORN archive extracted to: ", destination)
  
  
  normalizePath(destination, winslash = "/")
}