# 03_segment_fields_pseudo_units.R
#
# Purpose:
#   Convert a continuous weight raster (w_smooth) into a set of “harvest units”
#   by:
#     1) segmentation into contiguous patches,
#     2) tiling very large patches with a fishnet,
#     3) enforcing min/max area per unit,
#     4) selecting and renormalizing “harvest units” up to a target area.
#
# Main APIs:
#   segment_fields_diag()
#   build_harvest_units()
#   prune_harvest_units_by_weight()
#
# Assumptions:
#   - w_smooth is a single-layer SpatRaster in an equal-area CRS (e.g. EPSG:5070)
#   - w_smooth has NA outside the eligible crop mask
#   - the layer name is “w_norm” (for the zonal summarization + extract)

# --------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------

#' Create a fishnet grid over a polygon
#'
#' @param poly_sf sf polygon(s)
#' @param tile_size_m numeric tile side length in meters
#' @return sf with square tiles covering the bounding box of poly_sf
make_fishnet <- function(poly_sf, tile_size_m) {
  bb <- sf::st_bbox(poly_sf)
  
  xs <- seq(bb["xmin"], bb["xmax"], by = tile_size_m)
  ys <- seq(bb["ymin"], bb["ymax"], by = tile_size_m)
  
  cells <- vector("list", (length(xs) - 1) * (length(ys) - 1))
  k <- 1L
  for (i in seq_len(length(xs) - 1L)) {
    for (j in seq_len(length(ys) - 1L)) {
      x1 <- xs[i];   x2 <- xs[i + 1L]
      y1 <- ys[j];   y2 <- ys[j + 1L]
      cells[[k]] <- sf::st_polygon(list(rbind(
        c(x1, y1),
        c(x2, y1),
        c(x2, y2),
        c(x1, y2),
        c(x1, y1)
      )))
      k <- k + 1L
    }
  }
  
  grid <- sf::st_sfc(cells, crs = sf::st_crs(poly_sf))
  sf::st_sf(tile_id = seq_along(grid), geometry = grid)
}

#' Summarize weight and area inside a tile polygon
#'
#' @param tile_sf 1-row sf polygon
#' @param w_smooth SpatRaster of weights (same CRS as tile_sf)
#' @param cell_ha SpatRaster of cell areas in ha
#' @return list(area_ha, w_sum) or NULL if tile area <= 0
summarize_tile <- function(tile_sf, w_smooth, cell_ha) {
  v <- terra::vect(tile_sf)
  
  wv <- terra::values(terra::mask(w_smooth, v))
  av <- terra::values(terra::mask(cell_ha,  v))
  
  tile_area_ha <- sum(av, na.rm = TRUE)
  tile_w_sum   <- sum(wv, na.rm = TRUE)
  
  if (!is.finite(tile_area_ha) || tile_area_ha <= 0) return(NULL)
  
  list(area_ha = tile_area_ha,
       w_sum   = tile_w_sum)
}

# --------------------------------------------------------------------
# 1) Segment support into patches + pseudo-units
# --------------------------------------------------------------------

#' Segment fields into harvest units (patches + pseudo-fields)
#'
#' @param w_smooth SpatRaster (equal-area CRS, single band "w_norm")
#' @param max_field_ha_thresh numeric, max area for “small/medium” fields
#' @param tile_target_ha numeric, target tile area (ha) for giant fields
#' @param min_field_ha numeric, minimum allowed unit area (ha)
#' @param max_field_ha_final numeric, maximum allowed unit area (ha)
#'
#' @return list with patch_id raster, cell_ha raster, smallmed_poly_sf,
#'   pseudo_sf2, and all_units_sf (pre-harvest selection).
#' @export
segment_fields_diag <- function(w_smooth,
                                max_field_ha_thresh,
                                tile_target_ha,
                                min_field_ha,
                                max_field_ha_final) {
  
  # 0. Base area of support -------------------------------------------
  support_mask <- !is.na(w_smooth)
  support_mask <- terra::ifel(support_mask, 1, NA)
  
  cell_ha_raw <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha_sup <- terra::mask(cell_ha_raw, support_mask)
  
  support_area_ha <- terra::global(cell_ha_sup, "sum", na.rm = TRUE)[1, 1]
  cat(">> support (from mask) area (ha):", round(support_area_ha, 1), "\n")
  
  # 1. Segment into patches -------------------------------------------
  patch_id <- terra::patches(support_mask, directions = 4)
  
  cell_ha <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha <- terra::mask(cell_ha, patch_id)
  area_patches_before <- terra::global(cell_ha, "sum", na.rm = TRUE)[1, 1]
  cat(">> area in patches before sieve (ha):", round(area_patches_before, 1), "\n")
  
  # remove microscopic specks (< 50 cells)
  patch_id <- terra::sieve(patch_id, 50)
  
  cell_ha <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha <- terra::mask(cell_ha, patch_id)
  area_patches_after <- terra::global(cell_ha, "sum", na.rm = TRUE)[1, 1]
  cat(">> area in patches after sieve (ha): ", round(area_patches_after, 1), "\n")
  
  # 2. Zonal stats per patch ------------------------------------------
  df_area <- as.data.frame(
    terra::zonal(cell_ha, patch_id, fun = "sum", na.rm = TRUE)
  ) |>
    dplyr::rename(patch = patches, area_ha = area)
  
  df_wsum <- as.data.frame(
    terra::zonal(w_smooth, patch_id, fun = "sum", na.rm = TRUE)
  ) |>
    # assumes the band in w_smooth is called "w_norm"
    dplyr::rename(patch = patches, w_sum = w_norm)
  
  fields <- df_area |>
    dplyr::inner_join(df_wsum, by = "patch") |>
    dplyr::filter(
      is.finite(.data$area_ha), .data$area_ha > 0,
      is.finite(.data$w_sum)
    ) |>
    dplyr::mutate(
      w_density = .data$w_sum / .data$area_ha
    )
  
  cat(">> total area in fields (ha):",
      round(sum(fields$area_ha, na.rm = TRUE), 1), "\n")
  
  # 3. Split into small/medium vs giant -------------------------------
  fields_smallmed <- fields |>
    dplyr::filter(.data$area_ha <= max_field_ha_thresh)
  fields_giant <- fields |>
    dplyr::filter(.data$area_ha > max_field_ha_thresh)
  
  cat(">> small/med fields area (ha):",
      round(sum(fields_smallmed$area_ha, na.rm = TRUE), 1), "\n")
  cat(">> giant fields area (ha):   ",
      round(sum(fields_giant$area_ha,    na.rm = TRUE), 1), "\n")
  
  # keep only small/med patches in this mask
  keep_ids <- fields_smallmed$patch
  
  smallmed_mask <- patch_id
  smallmed_mask[!(patch_id[] %in% keep_ids)] <- NA
  
  smallmed_poly <- terra::as.polygons(smallmed_mask, dissolve = TRUE, minArea = 0)
  smallmed_poly <- smallmed_poly[!is.na(smallmed_poly$patches), ]
  smallmed_poly_sf <- sf::st_as_sf(smallmed_poly)
  
  smallmed_poly_sf <- smallmed_poly_sf |>
    dplyr::left_join(
      fields_smallmed |> dplyr::select(.data$patch, .data$area_ha, .data$w_sum),
      by = c("patches" = "patch")
    ) |>
    dplyr::rename(
      unit_id      = patches,
      unit_area_ha = area_ha,
      unit_w_sum   = w_sum
    )
  
  smallmed_poly_sf <- sf::st_set_crs(
    smallmed_poly_sf,
    sf::st_crs(terra::crs(w_smooth))
  )
  
  cat(">> area in smallmed_poly_sf (ha):",
      round(sum(smallmed_poly_sf$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  # 4. Tilings for giant fields ---------------------------------------
  tile_side_m <- sqrt(tile_target_ha * 10000)  # ha -> m² -> side length
  pseudo_list <- list()
  
  for (p in fields_giant$patch) {
    this_patch_mask <- terra::ifel(patch_id == p, 1, NA)
    this_poly       <- terra::as.polygons(this_patch_mask, dissolve = TRUE, minArea = 0)
    
    id_col <- names(this_poly)[1]
    this_poly <- this_poly[!is.na(this_poly[[id_col]]), ]
    if (nrow(this_poly) == 0) {
      cat("patch", p, " -> no polygon rows after mask\n")
      next
    }
    
    this_poly_sf <- sf::st_as_sf(this_poly)
    this_poly_sf <- sf::st_set_crs(this_poly_sf, sf::st_crs(terra::crs(patch_id)))
    this_poly_sf <- sf::st_make_valid(this_poly_sf)
    
    grid_sf <- make_fishnet(this_poly_sf, tile_side_m)
    grid_sf <- sf::st_make_valid(grid_sf)
    
    tiles_clipped <- suppressWarnings(sf::st_intersection(grid_sf, this_poly_sf))
    if (nrow(tiles_clipped) == 0) {
      cat("patch", p, " -> tiles_clipped empty after intersection\n")
      next
    }
    
    for (i in seq_len(nrow(tiles_clipped))) {
      tile_i <- tiles_clipped[i, , drop = FALSE]
      s      <- summarize_tile(tile_i, w_smooth, cell_ha)
      if (is.null(s)) next
      
      pseudo_list[[length(pseudo_list) + 1L]] <- sf::st_sf(
        patch_parent = p,
        unit_id      = paste0(p, "_", i),
        unit_area_ha = s$area_ha,
        unit_w_sum   = s$w_sum,
        geometry     = sf::st_geometry(tile_i),
        crs          = sf::st_crs(this_poly_sf)
      )
    }
    
    cat("patch", p, ": tiles added, pseudo_list length =",
        length(pseudo_list), "\n")
  }
  
  if (length(pseudo_list) > 0) {
    pseudo_sf2 <- do.call(rbind, pseudo_list)
  } else {
    pseudo_sf2 <- sf::st_sf(
      patch_parent = character(0),
      unit_id      = character(0),
      unit_area_ha = numeric(0),
      unit_w_sum   = numeric(0),
      geometry     = sf::st_sfc(crs = sf::st_crs(terra::crs(patch_id)))
    )
  }
  
  pseudo_sf2 <- sf::st_set_crs(pseudo_sf2, sf::st_crs(smallmed_poly_sf))
  
  cat(">> area in pseudo_sf2 (ha):",
      round(sum(pseudo_sf2$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  smallmed_poly_sf$unit_id <- as.character(smallmed_poly_sf$unit_id)
  pseudo_sf2$unit_id       <- as.character(pseudo_sf2$unit_id)
  
  # 5. Combine and filter units ---------------------------------------
  all_units0 <- dplyr::bind_rows(
    smallmed_poly_sf[, c("unit_id","unit_area_ha","unit_w_sum","geometry")],
    pseudo_sf2[,        c("unit_id","unit_area_ha","unit_w_sum","geometry")]
  ) |>
    dplyr::mutate(
      unit_w_density = .data$unit_w_sum / .data$unit_area_ha
    )
  
  cat(">> area in all_units0 (before size filter, ha):",
      round(sum(all_units0$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  all_units_sf <- all_units0 |>
    dplyr::filter(
      .data$unit_area_ha >= min_field_ha,
      .data$unit_area_ha <= max_field_ha_final
    ) |>
    dplyr::arrange(dplyr::desc(.data$unit_w_sum)) |>
    dplyr::mutate(
      cum_area = cumsum(.data$unit_area_ha)
    )
  
  cat(">> area in all_units_sf (after size filter, ha):",
      round(sum(all_units_sf$unit_area_ha, na.rm = TRUE), 1), "\n")
  cat(">> max cum_area (ha):",
      round(max(all_units_sf$cum_area, na.rm = TRUE), 1), "\n")
  
  list(
    patch_id         = patch_id,
    cell_ha          = cell_ha,
    smallmed_poly_sf = smallmed_poly_sf,
    pseudo_sf2       = pseudo_sf2,
    all_units_sf     = all_units_sf
  )
}

# --------------------------------------------------------------------
# 2) Build harvest subset + renormalize weights
# --------------------------------------------------------------------

#' Build capped / renormalized harvest units
#'
#' @param all_units_sf An sf object with columns
#'   `unit_id`, `unit_area_ha`, `unit_w_sum`, and `geometry`.
#' @param w_smooth SpatRaster of weights (same CRS as units)
#' @param target_ha numeric, NASS harvested hectares target for this county-year
#'
#' @return list(harvest_units_sf, sel_units_full, w_cap)
#' @export
build_harvest_units <- function(all_units_sf, w_smooth, target_ha) {
  stopifnot(inherits(all_units_sf, "sf"))
  req <- c("unit_id","unit_area_ha","unit_w_sum","geometry")
  if (!all(req %in% names(all_units_sf))) {
    stop("all_units_sf must contain: ", paste(req, collapse = ", "))
  }
  if (!inherits(w_smooth, "SpatRaster")) stop("w_smooth must be a SpatRaster")
  if (!is.numeric(target_ha) || target_ha <= 0) stop("target_ha must be > 0")
  
  # rank by unit_w_sum (cotton strength) and compute cumulative area
  harvest_cands <- all_units_sf |>
    dplyr::arrange(dplyr::desc(unit_w_sum)) |>
    dplyr::mutate(cum_area = cumsum(unit_area_ha))
  
  A_total  <- sum(harvest_cands$unit_area_ha, na.rm = TRUE)
  A_target <- min(target_ha, A_total)
  
  k <- which(harvest_cands$cum_area >= A_target)[1]
  
  sel_units_full <- if (is.na(k)) {
    harvest_cands$unit_id
  } else if (k > 1L) {
    harvest_cands$unit_id[seq_len(k - 1L)]
  } else {
    character(0)
  }
  
  if (length(sel_units_full) == 0 && nrow(harvest_cands) > 0) {
    sel_units_full <- harvest_cands$unit_id[1]
  }
  
  harvest_units_sf <- all_units_sf |>
    dplyr::filter(rlang::.data$unit_id %in% sel_units_full)
  
  if (nrow(harvest_units_sf) == 0) {
    stop("No units selected for harvest subset; check inputs/target_ha.")
  }
  
  # rasterize selected polygons -> mask -> normalize to sum ~ 1
  sel_vect <- terra::vect(harvest_units_sf)
  sel_mask <- w_smooth * NA
  sel_mask <- terra::rasterize(sel_vect, sel_mask, field = 1, touches = TRUE)
  
  w_cap <- terra::mask(w_smooth, sel_mask)
  s     <- as.numeric(terra::global(w_cap, "sum", na.rm = TRUE))
  if (!is.na(s) && s > 0) w_cap <- w_cap / s
  
  # push w_cap back to polygons: sum of pixel weights per unit
  cap_extract <- terra::extract(w_cap, sel_vect, fun = "sum", na.rm = TRUE)
  # assumes the band is named "w_norm"
  harvest_units_sf$unit_w_sum_cap <- cap_extract[, "w_norm"]
  
  total_cap <- sum(harvest_units_sf$unit_w_sum_cap, na.rm = TRUE)
  if (!is.finite(total_cap) || total_cap <= 0) {
    stop("Capped weight per-unit sums to 0; inspect w_cap/selection.")
  }
  harvest_units_sf$unit_w_sum_cap <- harvest_units_sf$unit_w_sum_cap / total_cap
  
  harvest_units_sf <- harvest_units_sf |>
    dplyr::mutate(unit_w_sum = .data$unit_w_sum_cap) |>
    dplyr::select(-.data$unit_w_sum_cap)
  
  list(
    harvest_units_sf = harvest_units_sf,
    sel_units_full   = sel_units_full,
    w_cap            = w_cap
  )
}

# --------------------------------------------------------------------
# 3) Optional pruning by cumulative weight
# --------------------------------------------------------------------

#' Prune harvest units by cumulative weight
#'
#' @param harvest_units_sf sf with column `unit_w_sum` (summing to ~1)
#' @param keep_cumw numeric, target cumulative weight to keep (e.g. 0.95)
#' @param min_units integer, minimum number of units to retain
#'
#' @return sf with a subset of units whose unit_w_sum has been renormalized.
#' @export
prune_harvest_units_by_weight <- function(harvest_units_sf,
                                          keep_cumw = 0.95,
                                          min_units = 2) {
  stopifnot("unit_w_sum" %in% names(harvest_units_sf))
  
  df <- harvest_units_sf |>
    dplyr::arrange(dplyr::desc(rlang::.data$unit_w_sum)) |>
    dplyr::mutate(cumw = cumsum(rlang::.data$unit_w_sum))
  
  k <- which(df$cumw >= keep_cumw)[1]
  if (is.na(k)) k <- nrow(df)
  k <- max(k, min_units)
  
  kept <- df[seq_len(k), ]
  kept$unit_w_sum <- kept$unit_w_sum / sum(kept$unit_w_sum, na.rm = TRUE)
  kept
}
