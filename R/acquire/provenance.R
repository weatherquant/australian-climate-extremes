# ============================================================
# provenance.R
#
# Dataset provenance utilities
# ============================================================


git_commit_hash <- function() {
  result <- tryCatch(
    system2(
      "git",
      args = c("rev-parse", "--short", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    ),
    
    error = function(e)
      NA_character_
    
  )
  
  if (length(result) == 0) {
    return(NA_character_)
  }
  
  result[[1]]
}



create_provenance_record <- function(dataset_id,
                                     source_url,
                                     source_version,
                                     path,
                                     status) {
  if (!file.exists(path)) {
    stop("Cannot create provenance record. File does not exist: ",
         path,
         call. = FALSE)
    
  }
  
  info <- file.info(path)
  
  tibble::tibble(
    dataset_id = dataset_id,
    
    retrieved_at = Sys.time(),
    
    source_url = source_url,
    
    source_version = source_version,
    
    filename = basename(path),
    
    local_path = normalizePath(path, winslash = "/", mustWork = TRUE),
    
    bytes = info$size,
    
    md5 = unname(tools::md5sum(path)),
    
    git_commit = git_commit_hash(),
    
    status = status
    
  )
}