# ============================================================
# match_acorn_raw_homogenised.R
#
# Match the unadjusted ACORN primary-site composite against
# the published homogenised ACORN-SAT daily series.
#
# Adjustment convention:
#
#   adjustment = homogenised - raw
#
# ============================================================


match_acorn_raw_homogenised <- function(raw_composite, homogenised) {
  required_raw <- c(
    "acorn_station_number",
    "date",
    "source_station_number",
    "primary_site_order",
    "tmax",
    "tmin"
  )
  
  required_hom <- c("acorn_station_number",
                    "station_name",
                    "date",
                    "tmax_hom",
                    "tmin_hom")
  
  missing_raw <- setdiff(required_raw, names(raw_composite))
  
  missing_hom <- setdiff(required_hom, names(homogenised))
  
  if (length(missing_raw) > 0L) {
    stop("Raw composite is missing columns: ",
         paste(missing_raw, collapse = ", "),
         call. = FALSE)
  }
  
  if (length(missing_hom) > 0L) {
    stop("Homogenised data are missing columns: ",
         paste(missing_hom, collapse = ", "),
         call. = FALSE)
  }
  
  
  # ----------------------------------------------------------
  # Keep only fields needed for comparison
  # ----------------------------------------------------------
  
  raw <- raw_composite[, .(
    acorn_station_number,
    date,
    source_station_number,
    primary_site_order,
    
    tmax_raw = tmax,
    tmin_raw = tmin
  )]
  
  
  hom <- homogenised[, .(acorn_station_number, station_name, date, tmax_hom, tmin_hom)]
  
  
  # ----------------------------------------------------------
  # Protect against duplicate join keys
  # ----------------------------------------------------------
  
  raw_duplicates <- raw[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  hom_duplicates <- hom[, .N, by = .(acorn_station_number, date)][N > 1L]
  
  if (nrow(raw_duplicates) > 0L) {
    stop("Raw composite contains duplicate station-date keys.",
         call. = FALSE)
  }
  
  if (nrow(hom_duplicates) > 0L) {
    stop("Homogenised data contain duplicate station-date keys.",
         call. = FALSE)
  }
  
  
  # ----------------------------------------------------------
  # Full join
  #
  # Deliberately retain dates that occur in either source.
  # ----------------------------------------------------------
  
  result <- merge(raw,
                  hom,
                  by = c("acorn_station_number", "date"),
                  all = TRUE)
  
  
  # ----------------------------------------------------------
  # Availability and adjustment variables
  # ----------------------------------------------------------
  
  result[, `:=`(
    tmax_adjustment =
      tmax_hom - tmax_raw,
    
    tmin_adjustment =
      tmin_hom - tmin_raw,
    
    raw_tmax_available =
      !is.na(tmax_raw),
    
    hom_tmax_available =
      !is.na(tmax_hom),
    
    raw_tmin_available =
      !is.na(tmin_raw),
    
    hom_tmin_available =
      !is.na(tmin_hom),
    
    matched_tmax =
      !is.na(tmax_raw) &
      !is.na(tmax_hom),
    
    matched_tmin =
      !is.na(tmin_raw) &
      !is.na(tmin_hom)
  )]
  
  
  # ----------------------------------------------------------
  # Categorise matching state
  # ----------------------------------------------------------
  
  result[, tmax_match_status :=
           data.table::fcase(
             matched_tmax,
             "matched",
             
             raw_tmax_available &
               !hom_tmax_available,
             "raw_only",!raw_tmax_available &
               hom_tmax_available,
             "homogenised_only",
             
             default = "missing_both"
           )]
  
  
  result[, tmin_match_status :=
           data.table::fcase(
             matched_tmin,
             "matched",
             
             raw_tmin_available &
               !hom_tmin_available,
             "raw_only",!raw_tmin_available &
               hom_tmin_available,
             "homogenised_only",
             
             default = "missing_both"
           )]
  
  
  data.table::setorder(result, acorn_station_number, date)
  
  result
}


# ============================================================
# Station-level adjustment summary
# ============================================================

summarise_acorn_adjustments <- function(comparison) {
  comparison[, .(
    first_date = min(date, na.rm = TRUE),
    
    last_date = max(date, na.rm = TRUE),
    
    # ----------------------------
    # Tmax coverage
    # ----------------------------
    
    n_raw_tmax =
      sum(!is.na(tmax_raw)),
    
    n_hom_tmax =
      sum(!is.na(tmax_hom)),
    
    n_matched_tmax =
      sum(matched_tmax),
    
    # ----------------------------
    # Tmin coverage
    # ----------------------------
    
    n_raw_tmin =
      sum(!is.na(tmin_raw)),
    
    n_hom_tmin =
      sum(!is.na(tmin_hom)),
    
    n_matched_tmin =
      sum(matched_tmin),
    
    # ----------------------------
    # Tmax adjustments
    # ----------------------------
    
    mean_tmax_adjustment =
      mean(tmax_adjustment, na.rm = TRUE),
    
    median_tmax_adjustment =
      stats::median(tmax_adjustment, na.rm = TRUE),
    
    sd_tmax_adjustment =
      stats::sd(tmax_adjustment, na.rm = TRUE),
    
    max_abs_tmax_adjustment =
      max(abs(tmax_adjustment), na.rm = TRUE),
    
    proportion_tmax_changed =
      mean(abs(tmax_adjustment) > 1e-10, na.rm = TRUE),
    
    # ----------------------------
    # Tmin adjustments
    # ----------------------------
    
    mean_tmin_adjustment =
      mean(tmin_adjustment, na.rm = TRUE),
    
    median_tmin_adjustment =
      stats::median(tmin_adjustment, na.rm = TRUE),
    
    sd_tmin_adjustment =
      stats::sd(tmin_adjustment, na.rm = TRUE),
    
    max_abs_tmin_adjustment =
      max(abs(tmin_adjustment), na.rm = TRUE),
    
    proportion_tmin_changed =
      mean(abs(tmin_adjustment) > 1e-10, na.rm = TRUE)
  ), by = .(acorn_station_number, station_name)]
}


# ============================================================
# Network-wide matching audit
# ============================================================

audit_acorn_matching <- function(comparison) {
  comparison[, .(
    rows = .N,
    
    stations =
      data.table::uniqueN(acorn_station_number),
    
    matched_tmax =
      sum(tmax_match_status == "matched"),
    
    raw_only_tmax =
      sum(tmax_match_status == "raw_only"),
    
    homogenised_only_tmax =
      sum(tmax_match_status ==
            "homogenised_only"),
    
    matched_tmin =
      sum(tmin_match_status == "matched"),
    
    raw_only_tmin =
      sum(tmin_match_status == "raw_only"),
    
    homogenised_only_tmin =
      sum(tmin_match_status ==
            "homogenised_only")
  )]
}