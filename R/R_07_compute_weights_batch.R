
#' Batch compute NAYD weights for a set of counties and years
#'
#' @param nass_df data.frame with columns county, year, area_planted_ha, area_harvested_ha.
#' @param years integer vector of years.
#' @param counties_sf sf with county polygons and a COUNTY column.
#' @param crop_name character.
#' @param crop_code integer CDL code.
#' @param sensor character, e.g. "Landsat".
#' @param season_start,season_end character "MM-DD".
#' @param CRS_TARGET character CRS string for the template.
#' @param lower_thresh Lower bound for CDL/NASS area ratio in mask QC.
#' @param upper_thresh Upper bound for CDL/NASS area ratio in mask QC.
#' @param area_buffer_frac Fractional buffer for NDVI eligibility area (e.g. 0.2).
#' @export
compute_weights_batch <- function(
    nass_df,
    years,
    counties_sf,
    crop_name,
    crop_code,
    sensor,
    season_start,
    season_end,
    CRS_TARGET = "EPSG:3857",
    lower_thresh    = 0.5,
    upper_thresh    = 3.5,
    area_buffer_frac = 0.2
) {
  years <- as.integer(years)
  
  for (yr in years) {
    cdl_fp_year <- cdl_for_processing(yr)
    if (is.na(cdl_fp_year)) {
      message("Year ", yr, ": no CDL mosaic, skipping.")
      next
    }
    cdl_TX <- terra::rast(cdl_fp_year)
    
    # Use all counties in counties_sf, or a subset if you pass it in filtered
    cnties <- unique(counties_sf$COUNTY)
    
    for (c in cnties) {
      cname     <- norm_county(c)
      county_sf <- counties_sf[counties_sf$COUNTY == c, ]
      
      out_fp <- .cache_path("weight", cname, yr, norm_county(crop_name), sensor)
      if (file.exists(out_fp)) {
        message("✓ skip existing weight: ", basename(out_fp))
        next
      }
      
      row_i <- nass_df[nass_df$county == cname & nass_df$year == yr, , drop = FALSE]
      if (nrow(row_i) == 0L) {
        message("– skip ", cname, " ", yr, ": no NASS row")
        next
      }
      row_i <- row_i[which.max(row_i$area_harvested_ha), , drop = FALSE]
      
      NASS_planted_ha  <- row_i$area_planted_ha
      NASS_harvest_ha  <- row_i$area_harvested_ha
      if (!is.finite(NASS_harvest_ha) || NASS_harvest_ha <= 0) {
        message("– skip ", cname, " ", yr, ": no NASS harvested area")
        next
      }
      
      wt <- tryCatch(
        get_weights_for_county_year(
          county_sf            = county_sf,
          county_name          = cname,
          sensor               = sensor,
          year                 = yr,
          season_start         = season_start,
          season_end           = season_end,
          CRS_TARGET           = CRS_TARGET,
          crop_code            = crop_code,
          crop_name            = crop_name,
          a_ndvi               = 0.7,
          a_et                 = 0.3,
          NASS_planted_area_ha = NASS_planted_ha,   # fixed name
          NASS_harvest_area_ha = NASS_harvest_ha,
          lower_thresh         = lower_thresh,
          upper_thresh         = upper_thresh,
          area_buffer_frac     = area_buffer_frac,
          cdl_dir              = cdl_dir,           # must be defined in compute_weights_batch()
          cache_dir            = cache_dir          # must be defined in compute_weights_batch()
        ),
        error = function(e) {
          message("– ERROR ", cname, " ", yr, ": ", conditionMessage(e))
          NULL
        }
      )
      
      if (is.null(wt) || !("weight" %in% names(wt))) {
        message("– skip ", cname, " ", yr, ": weights NULL")
        next
      }
      
      wsum <- as.numeric(terra::global(wt$weight, "sum", na.rm = TRUE))
      message("✓ ", cname, " ", yr,
              " (weights sum = ", sprintf("%.6f", wsum), 
              ", NASS_harvest_ha = ", round(NASS_harvest_ha, 1), ")")
    }
  }
  
  invisible(TRUE)
}
