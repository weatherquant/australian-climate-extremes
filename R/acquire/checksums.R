# ============================================================
# checksums.R
#
# File integrity utilities
# ============================================================


file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot calculate checksum. File does not exist: ",
         path,
         call. = FALSE)
  }
  
  unname(tools::md5sum(path))
}


file_inventory <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }
  
  info <- file.info(path)
  
  tibble::tibble(
    filename = basename(path),
    
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    
    bytes = info$size,
    
    modified = info$mtime,
    
    md5 = file_md5(path)
    
  )
}