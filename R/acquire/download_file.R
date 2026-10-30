# ============================================================
# download_file.R
#
# Safe and reproducible file downloading
# ============================================================


download_file_checked <- function(url,
                                  dest,
                                  overwrite = FALSE,
                                  quiet = FALSE) {
  dir.create(dirname(dest),
             recursive = TRUE,
             showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Existing file
  # ----------------------------------------------------------
  
  if (file.exists(dest) && !overwrite) {
    message("File already exists: ", dest)
    
    info <- file.info(dest)
    
    return(
      tibble::tibble(
        url = url,
        
        destination = normalizePath(dest, winslash = "/", mustWork = TRUE),
        
        status = "existing",
        
        bytes = info$size,
        
        md5 = unname(tools::md5sum(dest)),
        
        retrieved_at = as.POSIXct(NA)
        
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Temporary download
  # ----------------------------------------------------------
  
  tmp <- paste0(dest, ".part")
  
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  
  
  message("Downloading: ", url)
  
  
  curl::curl_download(
    url = url,
    
    destfile = tmp,
    
    quiet = quiet,
    
    mode = "wb"
    
  )
  
  
  # ----------------------------------------------------------
  # Validate download
  # ----------------------------------------------------------
  
  if (!file.exists(tmp)) {
    stop("Download failed: temporary file was not created.", call. = FALSE)
    
  }
  
  
  size <- file.info(tmp)$size
  
  
  if (is.na(size) || size <= 0) {
    unlink(tmp)
    
    stop("Download failed: downloaded file is empty.", call. = FALSE)
    
  }
  
  
  # ----------------------------------------------------------
  # Atomic-ish rename
  # ----------------------------------------------------------
  
  success <- file.rename(tmp, dest)
  
  
  if (!success) {
    stop(
      "Download completed but temporary file could not ",
      "be renamed to the destination.",
      call. = FALSE
    )
    
  }
  
  
  info <- file.info(dest)
  
  
  message("Download complete: ",
          basename(dest),
          " (",
          round(info$size / 1024^2, 2),
          " MB)")
  
  
  tibble::tibble(
    url = url,
    
    destination = normalizePath(dest, winslash = "/", mustWork = TRUE),
    
    status = "downloaded",
    
    bytes = info$size,
    
    md5 = unname(tools::md5sum(dest)),
    
    retrieved_at = Sys.time()
    
  )
}