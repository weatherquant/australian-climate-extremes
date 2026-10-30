# ============================================================
# read_acorn_primary_sites.R
#
# Parse ACORN-SAT primarysites.txt.
# ============================================================


find_acorn_primary_sites_member <- function(zip_file) {
  inventory <- unzip(zip_file, list = TRUE)
  
  member <- inventory$Name[grepl("primarysites[.]txt$", inventory$Name, ignore.case = TRUE)]
  
  if (length(member) != 1L) {
    stop("Expected exactly one primarysites.txt file; found ",
         length(member),
         ".",
         call. = FALSE)
  }
  
  member[[1]]
}


read_acorn_primary_sites <- function(zip_file, source_version = "2.3") {
  member <- find_acorn_primary_sites_member(zip_file)
  
  connection <- unz(zip_file, member, open = "rt")
  
  on.exit(close(connection), add = TRUE)
  
  lines <- readLines(connection, warn = FALSE)
  
  sites <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    col.names = c(
      "acorn_station_number",
      "site_1_number",
      "site_1_start",
      "site_1_end",
      "site_2_number",
      "site_2_start",
      "site_2_end",
      "site_3_number",
      "site_3_start",
      "site_3_end"
    ),
    colClasses = "character",
    showProgress = FALSE
  )
  
  date_columns <- c(
    "site_1_start",
    "site_1_end",
    "site_2_start",
    "site_2_end",
    "site_3_start",
    "site_3_end"
  )
  
  for (column in date_columns) {
    value <- sites[[column]]
    
    value[value == "99999999"] <- NA_character_
    
    sites[, (column) :=
            as.Date(value, format = "%Y%m%d")]
  }
  
  station_columns <- c("site_1_number", "site_2_number", "site_3_number")
  
  for (column in station_columns) {
    value <- sites[[column]]
    
    value[value == "999999"] <- NA_character_
    
    sites[, (column) := value]
  }
  
  sites[, `:=`(
    source_dataset = "acorn_primary_sites",
    source_version = source_version,
    source_member = member
  )]
  
  sites
}


acorn_primary_sites_long <- function(primary_sites) {
  required <- c(
    "acorn_station_number",
    "site_1_number",
    "site_1_start",
    "site_1_end",
    "site_2_number",
    "site_2_start",
    "site_2_end",
    "site_3_number",
    "site_3_start",
    "site_3_end"
  )
  
  missing_columns <- setdiff(required, names(primary_sites))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Primary-site table is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  site_1 <- primary_sites[, .(
    acorn_station_number,
    primary_site_order = 1L,
    source_station_number = site_1_number,
    start_date = site_1_start,
    end_date = site_1_end
  )]
  
  site_2 <- primary_sites[, .(
    acorn_station_number,
    primary_site_order = 2L,
    source_station_number = site_2_number,
    start_date = site_2_start,
    end_date = site_2_end
  )]
  
  site_3 <- primary_sites[, .(
    acorn_station_number,
    primary_site_order = 3L,
    source_station_number = site_3_number,
    start_date = site_3_start,
    end_date = site_3_end
  )]
  
  result <- data.table::rbindlist(list(site_1, site_2, site_3))
  
  result <- result[!is.na(source_station_number)]
  
  data.table::setorder(result, acorn_station_number, start_date)
  
  result
}