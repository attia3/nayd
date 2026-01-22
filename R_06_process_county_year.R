# 06_process_county_year.R
#
# Purpose:
#   High-level wrapper: for a given county-year,
#   1) take weight raster,
#   2) segment into harvest units,
#   3) disaggregate NASS yield,
#   4) cluster into production zones,
#   5) return cluster table and QA info.
#
# Main API:
#   process_county_year(
#     w_norm_orig,
#     county_geom,
#     max_field_ha_thresh,
#     tile_target_ha,
#     min_field_ha,
#     max_field_ha_final,
#     cluster_dist_m,
#     nass_area_ha,
#     nass_yield_kg_ha,
#     county_name,
#     year_val,
#     beta,
#     min_y_frac,
#     max_y_frac,
#     min_abs_dynamic
#   )

process_county_year <- function(
    w_norm_orig,        # SpatRaster for this county+year (EPSG:5070)
    county_geom,
    max_field_ha_thresh,
    tile_target_ha,
    min_field_ha,
    max_field_ha_final,
    cluster_dist_m, # sf polygon of county in EPSG:5070
    nass_area_ha,
    nass_yield_kg_ha,
    county_name,
    crop_name,
    year_val,
    beta = NULL,
    min_frac = NULL,
    max_frac = NULL,
    # absolute crop-specific caps
    y_min_crop = NULL,   # e.g. cotton floor
    y_max_crop = NULL,  # e.g. cotton ceiling
    min_abs  = NULL       # legacy absolute floor
) {
  # crop/mask raster to county boundary just to be safe
  
  w_norm_orig <- terra::project(w_norm_orig,"EPSG:5070",method = "bilinear")
  county_geom <- sf::st_transform(county_geom, crs = terra::crs(w_norm_orig))
  
  w_norm_crop <- terra::crop(w_norm_orig, terra::vect(county_geom))
  w_norm_crop <- terra::mask(w_norm_crop, terra::vect(county_geom))
  
  # step 2: segment fields and tile giants
  seg <- segment_fields_diag(
    w_smooth          = w_norm_crop,
    max_field_ha_thresh = max_field_ha_thresh,
    tile_target_ha      = tile_target_ha,
    min_field_ha        = min_field_ha,
    max_field_ha_final  = max_field_ha_final
  )
  
  all_units_sf <- seg$all_units_sf
  
  # step 3: pick harvested units to match NASS area,
  # build w_cap, and push renormalized weights back to polygons
  harv <- build_harvest_units(
    all_units_sf = all_units_sf,
    w_smooth     = w_norm_crop,
    target_ha    = nass_area_ha    # harvested ha target
  )
  
  harvest_units_sf <- prune_harvest_units_by_weight(harv$harvest_units_sf,
                                                    nass_harvest_area_ha = nass_area_ha,  
                                                    min_units = 2)
  # w_cap <- harv$w_cap  # available if you want to save the raster
  
  # step 4: disaggregate NASS yield into each harvested unit (safe version)
  min_abs_dynamic <- max(80, 0.2 * nass_yield_kg_ha)
  dis <- disagg_units_yield_safe(
    units_sf         = harvest_units_sf,
    nass_area_ha     = nass_area_ha,
    nass_yield_kg_ha = nass_yield_kg_ha,
    county           = county_name,
    year             = year_val,
    beta             = beta,
    crop             = crop_name,
    to_crs           = "EPSG:4326",
    min_frac         = min_frac,
    max_frac         = max_frac,
    # absolute crop-specific caps
    y_min_crop = y_min_crop,   # e.g. cotton floor
    y_max_crop = y_max_crop,  # e.g. cotton ceiling
    min_abs          = min_abs
  )
  
  units_with_yield <- dis$units_with_yield
  qa_check         <- dis$check
  
  # step 5: cluster units into ~5 km production zones
  cl <- cluster_units_5km(units_with_yield, cluster_dist_m = cluster_dist_m)
  
  eps <- 1e-6
  floors_before <- sum(round(cl$cluster_sf$cluster_yield_kg_ha,6) <= round(min_abs_dynamic,6)+eps)
  
  merged <- merge_excess_floor_clusters(
    clusters_sf     = cl$cluster_sf,
    min_abs         = min_abs_dynamic,
    keep_floor      = 1,           # keep only ONE floor cluster
    prefer_nonfloor = TRUE,
    county          = county_name,
    year            = year_val,
    add_lonlat      = TRUE
  )
  
  merged$cluster_table <- merged$cluster_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(county, year, lon, lat, cluster_area_ha, cluster_prod_kg, cluster_yield_kg_ha, n_units)
  
  clusters_sf_final  <- merged$cluster_sf
  clusters_table_out <- merged$cluster_table |> 
    dplyr::mutate(
      nass_area_ha = nass_area_ha,
      nass_yield_kg_ha = nass_yield_kg_ha
    )
  
  floors_after <- sum(round(clusters_sf_final$cluster_yield_kg_ha,6) <= round(min_abs_dynamic,6)+eps)
  
  # --- Print status line ---
  n_units    <- nrow(units_with_yield)
  n_clusters <- nrow(clusters_sf_final)
  cat(sprintf(
    "[%s %d]  NASS: area=%.0f ha | yield=%.1f kg/ha (%.1f kg/ha) | clusters=%d (floors: %d→%d) || calc_check: area=%.0f ha, mean=%.1f kg/ha\n",
    county_name, as.integer(year_val),
    nass_area_ha,
    nass_yield_kg_ha,
    nass_yield_kg_ha,
    n_clusters, floors_before, floors_after,
    sum(clusters_table_out$cluster_area_ha),
    mean(clusters_table_out$cluster_yield_kg_ha) 
  ))
  
  # --- Return a bundle you can append/save ---
  list(
    prep =list(
      w_smooth          = w_norm_crop,
      w_cap             = harv$w_cap %||% NULL  # if you saved it
    ),
    seg = list(
      all_units_sf      = seg$all_units_sf,
      smallmed_poly_sf  = seg$smallmed_poly_sf %||% NULL,
      pseudo_sf2        = seg$pseudo_sf2 %||% NULL
    ),
    harv = list(
      harvest_units_sf  = harvest_units_sf,
      sel_units_full    = harv$sel_units_full %||% NULL
    ),
    dis = list(
      units_with_yield  = units_with_yield,
      check             = qa_check
    ),
    clusters = list(
      clusters_sf = clusters_sf_final, # merged clusters (with lon/lat)
      clusters_table = clusters_table_out
    ),
    info = list(
      county   = county_name,
      year     = as.integer(year_val))
  )
  
}
