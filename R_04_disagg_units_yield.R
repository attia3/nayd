
# helper: default caps by crop ---

.get_crop_caps <- function(crop){
  crop_low <- tolower(crop)
  if (crop_low %in% c("cotton", "lint", "upland cotton")){
    list(y_min = 100, y_max = 4000, min_abs = 80)
  }else if(crop_low %in% c("winter wheat","wheat")){
    list(y_min = 500, y_max = 9000, min_abs = 400)
  }else if(crop_low %in% c("corn","maize")){
    list(y_min = 700, y_max = 16000, min_abs = 600)
  }else {
    # generic, very wide caps
    list(y_min = 100,   y_max = 10000, min_abs = 50)
  }
}

disagg_units_yield_safe <- function(units_sf,
                                    nass_area_ha,
                                    nass_yield_kg_ha,
                                    county,
                                    year,
                                    beta = NULL,
                                    crop = crop_name,
                                    to_crs = "EPSG:4326",
                                    # relative bounds (optional, can be relaxed)
                                    min_frac = NULL,
                                    max_frac = NULL,
                                    # OPTIONAL crop-specific caps; if NULL, use defaults
                                    y_min_crop = NULL,   # e.g. cotton floor
                                    y_max_crop = NULL,  # e.g. cotton ceiling
                                    min_abs  = NULL       # legacy absolute floor
) {
  stopifnot(inherits(units_sf, "sf"))
  needed_cols <- c("unit_id", "unit_area_ha", "unit_w_sum")
  if (!all(needed_cols %in% names(units_sf))) {
    stop("units_sf must have columns: unit_id, unit_area_ha, unit_w_sum")
  }
  
  ## --- 0) resolve crop-specific defaults ------------------------
  caps <- .get_crop_caps(crop)
  
  if (is.null(y_min_crop)) y_min_crop <- caps$y_min
  if (is.null(y_max_crop)) y_max_crop <- caps$y_max
  if (is.null(min_abs))    min_abs    <- caps$min_abs  
  
  ## 1) scale area to NASS (unchanged)
  area_total_raw <- sum(units_sf$unit_area_ha, na.rm = TRUE)
  if (!is.finite(area_total_raw) || area_total_raw <= 0)
    stop("Bad area total in units_sf")
  
  scaleA <- nass_area_ha / area_total_raw
  units_sf$area_ha_scaled <- units_sf$unit_area_ha * scaleA
  
  ## 2) weight with beta (unchanged)
  units_sf$w_beta <- (units_sf$unit_w_sum)^beta
  units_sf$w_beta[!is.finite(units_sf$w_beta) | units_sf$w_beta < 0] <- 0
  
  wb_mean <- sum(units_sf$w_beta * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  yield_raw <- nass_yield_kg_ha * (units_sf$w_beta / wb_mean)
  
  ## --- 3) Safeguards: crop caps (+ optional relative bounds) ---
  # relative bounds
  y_floor_rel <- nass_yield_kg_ha * min_frac
  y_ceiling_rel <- nass_yield_kg_ha * max_frac
  
  # combine: final floor = max(relative floor, legacy min_abs, crop minimum)
  y_floor <- max(y_min_crop, min_abs, y_floor_rel)
  
  # final ceiling = min(relative ceiling, crop maximum)
  # (so crop-specific cap cannot be broken by a high NASS)
  y_ceiling <- min(y_ceiling_rel, y_max_crop)
  
  # if for some bad county y_floor > y_ceiling, we shrink floor
  if (y_floor > y_ceiling) {
    warning("For ", county, " ", year,
            " floor > ceiling; relaxing floor to ceiling.")
    y_floor <- y_ceiling
  }
  
  yield_clamped <- pmin(pmax(yield_raw, y_floor), y_ceiling)
  
  ## 4) Recenter to NASS mean *if possible*
  aw_mean_clamped <- sum(yield_clamped * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  adj_factor <- if (is.finite(aw_mean_clamped) && aw_mean_clamped > 0) {
    nass_yield_kg_ha / aw_mean_clamped
  } else 1
  
  yield_adj <- yield_clamped * adj_factor
  
  # tiny final floor, same as before, but won’t exceed y_ceiling
  yield_final <- pmin(pmax(yield_adj, y_floor), y_ceiling)
  
  ## 5) production, centroids, metadata (unchanged)
  prod_final <- yield_final * units_sf$area_ha_scaled
  units_sf$yield_kg_ha <- yield_final
  units_sf$prod_kg     <- prod_final
  
  if (!is.null(to_crs)) {
    units_ll <- sf::st_transform(units_sf, to_crs)
  } else {
    units_ll <- units_sf
  }
  cents <- sf::st_centroid(units_ll)
  xy    <- sf::st_coordinates(cents)
  units_sf$lon <- xy[,1]
  units_sf$lat <- xy[,2]
  
  units_sf$county <- county
  units_sf$year   <- as.integer(year)
  
  ## 6) QA
  total_area_scaled <- sum(units_sf$area_ha_scaled, na.rm = TRUE)
  total_prod_kg     <- sum(units_sf$prod_kg,        na.rm = TRUE)
  mean_yield_check  <- total_prod_kg / total_area_scaled
  
  # how far are we from NASS after enforcing caps?
  nass_dev_pct <- 100 * (mean_yield_check - nass_yield_kg_ha) / nass_yield_kg_ha
  
  check <- data.frame(
    county = county,
    year   = year,
    nass_area_ha               = nass_area_ha,
    units_total_area_scaled_ha = total_area_scaled,
    nass_yield_kg_ha           = nass_yield_kg_ha,
    area_weighted_mean_yield   = mean_yield_check,
    expected_prod_kg           = nass_area_ha * nass_yield_kg_ha,
    units_total_prod_kg        = total_prod_kg,
    y_floor_used_kg_ha         = y_floor,
    y_ceiling_used_kg_ha       = y_ceiling,
    adj_factor_applied         = adj_factor,
    nass_deviation_pct         = nass_dev_pct
  )
  
  list(
    units_with_yield = units_sf %>%
      dplyr::select(
        unit_id,
        county,
        year,
        lon, lat,
        area_ha_scaled,
        yield_kg_ha,
        prod_kg,
        unit_area_ha,
        unit_w_sum,
        geometry
      ),
    check = check
  )
}

