#' Run NAYD disaggregation and clustering for multiple counties and years
#'
#' This is a high-level batch wrapper around `process_county_year()` that:
#' \enumerate{
#'   \item reads cached weight rasters for each county–year,
#'   \item segments harvest units and disaggregates NASS yields,
#'   \item clusters units into production zones,
#'   \item returns a stacked table of cluster-level yields.
#' }
#'
#' It assumes that weight rasters have already been created by
#' `get_weights_for_county_year()` and written via
#' `.cache_path("weight", ...)`.
#'
#' @param years Integer vector of years to process.
#' @param counties Character vector of county names (as they appear in
#'   `counties_sf$COUNTY`).
#' @param nass_table A data frame with at least the columns
#'   `county`, `year`, `area_harvested_ha`, and `yield_kg_ha`. County names
#'   should already be normalised in the same way used when creating the
#'   weight rasters.
#' @param counties_sf An `sf` object with county polygons and a `COUNTY`
#'   column matching the (unnormalised) county names used in `counties`.
#' @param crop_name Character crop label used when constructing cache paths
#'   (e.g., `"winter wheat"`, `"cotton"`).
#' @param sensor Character sensor tag used when constructing cache paths
#'   (e.g., `"LANDSAT"`).
#' @param betas Numeric vector of beta exponents to evaluate.
#' @param cluster_dist_m Maximum complete-linkage distance (in metres) used
#'   when clustering units into production zones.
#' @param max_field_ha_thresh,tile_target_ha,min_field_ha,max_field_ha_final
#'   Segmentation and tiling parameters passed on to `segment_fields_diag()`.
#' @param min_frac,max_frac Optional relative yield bounds passed through to
#'   `disagg_units_yield_safe()`. Use `NULL` to rely purely on crop-specific
#'   caps.
#' @param y_min_crop,y_max_crop,min_abs Optional crop-specific caps passed to
#'   `disagg_units_yield_safe()`. If `NULL`, defaults are chosen internally
#'   based on the crop.
#'
#' @return A tibble with one row per cluster–county–year–beta combination,
#'   containing cluster centroid coordinates, area, production, yield, and the
#'   corresponding NASS reference area and yield.
#'
#' @export

run_disagg_clustering_batch <- function(
    years,
    counties,
    nass_table,
    counties_sf,
    crop_name,
    sensor           = "LANDSAT",
    betas            = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8),
    cluster_dist_m   = 10000,
    max_field_ha_thresh = 300,
    tile_target_ha      = 200,
    min_field_ha        = 5,
    max_field_ha_final  = 500,
    min_frac            = NULL,
    max_frac            = NULL,
    y_min_crop          = NULL,
    y_max_crop          = NULL,
    min_abs             = NULL
) {
  years    <- sort(unique(as.integer(years)))
  counties <- sort(unique(as.character(counties)))
  
  # filter counties_sf up front
  counties_sf <- counties_sf[counties_sf$COUNTY %in% counties, ]
  
  # filter NASS table to relevant subset
  nass_sub <- nass_table[
    nass_table$county %in% norm_county(counties) &
      nass_table$year   %in% years,
    ,
    drop = FALSE
  ]
  
  if (nrow(nass_sub) == 0L) {
    stop("nass_table has no rows matching the requested counties/years.")
  }
  
  res_list <- list()
  idx <- 1L
  
  for (yr in years) {
    for (c in counties) {
      
      county_i <- c
      year_i   <- yr
      
      # 1. NASS row for this county–year ------------------------------
      nass_rows <- nass_sub[
        nass_sub$county == norm_county(county_i) &
          nass_sub$year   == year_i,
        ,
        drop = FALSE
      ]
      
      if (nrow(nass_rows) == 0L) {
        message("[SKIP] No NASS row for ", county_i, " ", year_i)
        next
      }
      
      # if multiple entries exist (e.g. different irrig statuses),
      # keep row with max harvested area
      nass_rows <- nass_rows[which.max(nass_rows$area_harvested_ha), , drop = FALSE]
      
      nass_area_ha_i     <- nass_rows$area_harvested_ha
      nass_yield_kg_ha_i <- nass_rows$yield_kg_ha
      
      if (!is.finite(nass_area_ha_i) || nass_area_ha_i <= 0) {
        message("[SKIP] ", county_i, " ", year_i, ": invalid NASS harvested area")
        next
      }
      
      # 2. County geometry (sf) ---------------------------------------
      county_geom_i <- counties_sf[counties_sf$COUNTY == county_i, ]
      if (nrow(county_geom_i) == 0L) {
        message("[SKIP] No county geometry for ", county_i)
        next
      }
      
      # 3. Weight raster path (cached by get_weights_for_county_year) --
      w_path_i <- .cache_path(
        kind      = "weight",
        county    = norm_county(county_i),
        year      = year_i,
        crop_name = crop_name,
        sensor    = sensor
      )
      
      if (!file.exists(w_path_i)) {
        message("[SKIP] No weight raster for ", county_i, " ", year_i,
                " at ", w_path_i)
        next
      }
      
      w_norm_orig_i <- terra::rast(w_path_i)
      
      # 4. Loop over beta values ---------------------------------------
      for (b in betas) {
        
        message(sprintf(
          "[INFO] %s %d – beta = %.2f",
          county_i, year_i, b
        ))
        
        out <- process_county_year(
          w_norm_orig        = w_norm_orig_i,
          county_geom        = county_geom_i,
          max_field_ha_thresh = max_field_ha_thresh,
          tile_target_ha      = tile_target_ha,
          min_field_ha        = min_field_ha,
          max_field_ha_final  = max_field_ha_final,
          cluster_dist_m      = cluster_dist_m,
          nass_area_ha        = nass_area_ha_i,
          nass_yield_kg_ha    = nass_yield_kg_ha_i,
          county_name         = county_i,
          crop_name           = crop_name,
          year_val            = year_i,
          beta                = b,
          min_frac            = min_frac,
          max_frac            = max_frac,
          y_min_crop          = y_min_crop,
          y_max_crop          = y_max_crop,
          min_abs             = min_abs
        )
        
        cl_tab <- out$clusters$clusters_table |>
          dplyr::mutate(
            beta   = b
          )
        
        res_list[[idx]] <- cl_tab
        idx <- idx + 1L
      }
    }
  }
  
  if (length(res_list) == 0L) {
    warning("No county–year–beta combinations produced output.")
    return(tibble::tibble())
  }
  
  dplyr::bind_rows(res_list)
}
