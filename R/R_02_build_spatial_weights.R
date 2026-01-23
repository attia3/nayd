
# 02_get_weights_county_year.R
#
# Build CDL-based crop mask and NDVI/ET-based weights for a single
# county-year and crop, aligned to a template grid.
#
# Main API:
#   get_weights_for_county_year(
#     county_sf, county_name,
#     sensor, year,
#     season_start, season_end,
#     CRS_TARGET,
#     crop_code, crop_name,
#     a_ndvi, a_et,
#     NASS_planted_area_ha, NASS_harvest_area_ha,
#     lower_thresh, upper_thresh,
#     area_buffer_frac,
#     cdl_dir, cache_dir
#   )
#
# `cdl_dir` should contain CDL mosaics named like "CDL_2016_48.tif" (for Texas = 48).
# `cache_dir` is where NDVI/ET/weight rasters will be written.
#
# Assumes:
#   - rgee is installed and ee_Initialize() has been called
#   - helpers `ndvi_season_image()` and `eemetric_et_season_image()` exist in the package.

# ----------------------------------------------------------------------
# small utilities

# ----------------------------------------------------------------------
# CDL-based crop mask builder (single county)
# ----------------------------------------------------------------------

make_crop_mask_year <- function(template,
                                county_v,
                                year,
                                NASS_planted_area_ha = NULL,
                                lower_thresh = NULL,
                                upper_thresh = NULL,
                                crop_code = NULL,
                                cdl_dir) {
  stopifnot(inherits(template, "SpatRaster"),
            inherits(county_v,  "SpatVector"))
  
  # path helper: use 2008 CDL for any year < 2008 (for pre-2008 years)
  cdl_path_for_year <- function(y) {
    y_eff <- if (y < 2008L) 2008L else as.integer(y)
    file.path(cdl_dir, sprintf("CDL_%d_48.tif", y_eff))
  }
  
  .to_bin <- function(r) terra::ifel(r == crop_code, 1, NA)
  
  .proj_crop_mask <- function(r) {
    terra::project(r, template, method = "near")
  }
  
  .area_ha <- function(binmask) {
    a <- terra::cellSize(binmask, unit = "m")
    a <- terra::mask(a, binmask)
    as.numeric(terra::global(a, "sum", na.rm = TRUE)) / 1e4
  }
  
  .load_cdl_bin <- function(y) {
    fp <- cdl_path_for_year(as.integer(y))
    if (!file.exists(fp)) return(NULL)
    r  <- try(terra::rast(fp), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    .to_bin(.proj_crop_mask(r))
  }
  
  .union_bin <- function(years_vec) {
    years_vec <- unique(as.integer(years_vec))
    rs <- lapply(years_vec, .load_cdl_bin)
    rs <- Filter(Negate(is.null), rs)
    if (length(rs) == 0) return(NULL)
    if (length(rs) == 1) return(rs[[1]])
    
    u <- do.call(terra::mosaic, c(rs, list(fun = "max")))
    terra::ifel(u == 1, 1, NA)
  }
  
  .intersect_bin <- function(years_vec) {
    years_vec <- unique(as.integer(years_vec))
    rs <- lapply(years_vec, .load_cdl_bin)
    rs <- Filter(Negate(is.null), rs)
    if (length(rs) == 0) return(NULL)
    if (length(rs) == 1) return(rs[[1]])
    
    b <- rs[[1]]
    for (k in 2:length(rs)) {
      b <- terra::ifel(!is.na(b) & !is.na(rs[[k]]), 1, NA)
    }
    b
  }
  
  .win3 <- function(year) {
    if (year < 2008L) 2008L:2010L else (year - 1L):(year + 1L)
  }
  
  used_strategy <- NA_character_
  
  mask_year <- if (year >= 2008L) .load_cdl_bin(year) else NULL
  
  if (!is.null(mask_year)) {
    if (!is.null(NASS_planted_area_ha) && is.finite(NASS_planted_area_ha)) {
      area_y <- .area_ha(mask_year)
      ratio  <- area_y / NASS_planted_area_ha
      
      # Case 1: CDL << NASS → 3-year union
      if (!is.na(area_y) && ratio < lower_thresh) {
        yrs_win <- .win3(year)
        mask3   <- .union_bin(yrs_win)
        if (!is.null(mask3)) {
          used_strategy <- sprintf("union_%s", paste(yrs_win, collapse = "_"))
          return(list(mask = mask3,
                      strategy = used_strategy,
                      area_ha = .area_ha(mask3)))
        }
      }
      # Case 2: CDL >> NASS → 3-year intersection
      if (!is.na(ratio) && ratio >= upper_thresh) {
        yrs_win    <- .win3(year)
        mask3_int  <- .intersect_bin(yrs_win)
        if (!is.null(mask3_int)) {
          used_strategy <- sprintf("intersect_%s", paste(yrs_win, collapse = "_"))
          return(list(
            mask     = mask3_int,
            strategy = used_strategy,
            area_ha  = .area_ha(mask3_int)
          ))
        }
      }
    }
    
    used_strategy <- sprintf("exact_%d", year)
    return(list(
      mask     = mask_year,
      strategy = used_strategy,
      area_ha  = .area_ha(mask_year)
    ))
  }
  
  # Fallback: 3-year union if exact-year fails
  yrs_win <- .win3(year)
  mask3   <- .union_bin(yrs_win)
  if (!is.null(mask3)) {
    used_strategy <- sprintf("union_%s", paste(yrs_win, collapse = "_"))
    return(list(
      mask     = mask3,
      strategy = used_strategy,
      area_ha  = .area_ha(mask3)
    ))
  }
  
  list(mask = NULL, strategy = "failed_to_build_mask", area_ha = NA_real_)
}

# ----------------------------------------------------------------------
# NDVI eligibility builder
# ----------------------------------------------------------------------

.ndvi_eligibility <- function(ndvi_r,
                              crop_mask,
                              NASS_harvest_area_ha,
                              area_buffer_frac = 0,
                              candidate_pcts = seq(0.90, 0.10, by = -0.02)) {
  stopifnot(inherits(ndvi_r, "SpatRaster"),
            inherits(crop_mask, "SpatRaster"))
  
  ndv <- terra::mask(ndvi_r, crop_mask)
  vals <- terra::values(ndv, mat = FALSE)
  if (all(is.na(vals))) {
    return(list(mask = NULL, threshold = NA_real_, area_ha = NA_real_))
  }
  
  a_m2 <- terra::cellSize(ndv, unit = "m")
  a_m2 <- terra::mask(a_m2, ndv)
  a_ha <- as.numeric(terra::values(a_m2, mat = FALSE)) / 1e4
  nd   <- as.numeric(vals)
  
  if (is.null(NASS_harvest_area_ha) || !is.finite(NASS_harvest_area_ha) || NASS_harvest_area_ha <= 0) {
    thr  <- stats::quantile(nd, probs = 0.50, na.rm = TRUE)
    keep <- nd >= thr
    return(list(
      mask      = terra::ifel(ndv >= thr, 1, NA),
      threshold = as.numeric(thr),
      area_ha   = sum(a_ha[keep], na.rm = TRUE)
    ))
  }
  
  target <- NASS_harvest_area_ha * (1 + area_buffer_frac)
  
  for (p in candidate_pcts) {
    thr  <- stats::quantile(nd, probs = p, na.rm = TRUE)
    keep <- nd >= thr
    A    <- sum(a_ha[keep], na.rm = TRUE)
    if (is.finite(A) && A >= target) {
      return(list(
        mask      = terra::ifel(ndv >= thr, 1, NA),
        threshold = as.numeric(thr),
        area_ha   = A
      ))
    }
  }
  
  thr  <- stats::quantile(nd, probs = utils::tail(candidate_pcts, 1), na.rm = TRUE)
  keep <- nd >= thr
  A    <- sum(a_ha[keep], na.rm = TRUE)
  list(
    mask      = terra::ifel(ndv >= thr, 1, NA),
    threshold = as.numeric(thr),
    area_ha   = A
  )
}

# ----------------------------------------------------------------------
# Main API: one county-year weights
# ----------------------------------------------------------------------

#' Build NDVI/ET-based weights for a single county-year
#'
#' @param county_sf sf polygon for the county (single feature)
#' @param county_name character name (used in logging and filenames)
#' @param sensor character, e.g. "Landsat", "Sentinel-2"
#' @param year integer year
#' @param season_start character, e.g. "04-01"
#' @param season_end   character, e.g. "09-30"
#' @param CRS_TARGET target CRS for processing (e.g. "EPSG:5070")
#' @param crop_code integer CDL code for crop
#' @param crop_name character crop label, used in cache dir
#' @param a_ndvi numeric weight for NDVI (0–1)
#' @param a_et   numeric weight for ET (0–1)
#' @param NASS_planted_area_ha numeric planted area (ha)
#' @param NASS_harvest_area_ha numeric harvested area (ha)
#' @param lower_thresh,upper_thresh thresholds for CDL/NASS area ratio
#' @param area_buffer_frac NDVI eligibility buffer (e.g. 0.05)
#' @param cdl_dir directory that holds CDL mosaics "CDL_<year>_48.tif"
#' @param cache_dir directory to write NDVI/ET/weight rasters
#'
#' @return list with weight raster, ndvi, et, masks, diagnostics.
#' @export
get_weights_for_county_year <- function(county_sf, county_name,
                                        sensor,
                                        year,
                                        season_start, season_end,
                                        CRS_TARGET,
                                        crop_code, crop_name,
                                        a_ndvi, a_et,
                                        NASS_planted_area_ha = NULL,
                                        NASS_harvest_area_ha = NULL,
                                        lower_thresh = 0.6,
                                        upper_thresh = 2.0,
                                        area_buffer_frac = 0.05,
                                        cdl_dir,
                                        cache_dir) {
  stopifnot(is.character(county_name), length(county_name) == 1)
  
  # ---- runtime dependency check for rgee ----
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop(
      "NAYD::get_weights_for_county_year() requires the 'rgee' package.\n",
      "Please install it with install.packages('rgee') and configure Earth Engine.",
      call. = FALSE
    )
  }
  
  # local alias for GEE handle
  ee <- rgee::ee
  
  mode <- rs_mode_for_year(year)
  
  aoi      <- rgee::sf_as_ee(sf::st_transform(county_sf, 4326))
  aoi_geom <- aoi$geometry()
  d1       <- ee$Date(sprintf("%d-%s", year, season_start))
  d2       <- ee$Date(sprintf("%d-%s", year, season_end))
  
  # --- CDL statewide -> county ---
  cdl_fp <- file.path(cdl_dir, sprintf("CDL_%d_48.tif", if (year < 2008L) 2008L else as.integer(year)))
  if (!file.exists(cdl_fp)) {
    message("Year ", year, ": no CDL mosaic at ", cdl_fp)
    return(NULL)
  }
  cdl_TX  <- terra::rast(cdl_fp)
  county_v <- county_sf |>
    sf::st_transform(terra::crs(cdl_TX)) |>
    terra::vect()
  
  cdl_cnty <- terra::mask(terra::crop(cdl_TX, county_v), county_v)
  cdl_cnty <- terra::project(cdl_cnty, CRS_TARGET)
  
  if (!terra::hasValues(cdl_cnty)) return(NULL)
  template <- cdl_cnty
  
  # --- Crop mask via CDL + NASS consistency ---
  res <- make_crop_mask_year(
    template            = template,
    county_v            = county_v,
    year                = year,
    NASS_planted_area_ha = NASS_planted_area_ha,
    lower_thresh        = lower_thresh,
    upper_thresh        = upper_thresh,
    crop_code           = crop_code,
    cdl_dir             = cdl_dir
  )
  
  crop_mask   <- res$mask
  mask_ha     <- res$area_ha
  nass_ha     <- NASS_harvest_area_ha
  nass_pl_ha  <- NASS_planted_area_ha
  
  pct_cov <- if (!is.null(nass_ha) && is.finite(nass_ha) && nass_ha > 0 &&
                 !is.null(mask_ha) && is.finite(mask_ha)) {
    100 * mask_ha / nass_ha
  } else {
    NA_real_
  }
  
  cat(
    sprintf(
      "✓ Built via %s for county %s, year %d — NASS planted area = %s ha, crop mask = %s ha%s\n",
      res$strategy,
      county_name,
      year,
      fmt_ha(nass_pl_ha),
      fmt_ha(mask_ha),
      if (is.finite(pct_cov)) sprintf(" (%.1f%% of NASS)", pct_cov) else ""
    )
  )
  
  if (is.null(crop_mask)) {
    cat("No mask:", res$strategy, "\n")
    return(NULL)
  }
  
  n_pix <- sum(terra::values(crop_mask, mat = FALSE) == 1, na.rm = TRUE)
  if (is.na(n_pix) || n_pix == 0) return(NULL)
  
  has_et <- (mode != "NDVI_ONLY_2008mask")
  
  # ------------------------------------------------------------------
  # ET (OpenET eeMETRIC)
  # ------------------------------------------------------------------
  et_cot     <- NULL
  et_screen  <- NULL
  scale_et   <- NA_integer_
  
  if (has_et) {
    et_obj <- eemetric_et_season_image(
      start = d1,
      end   = d2,
      geom  = aoi_geom
    )
    et_img   <- et_obj$et
    scale_et <- et_obj$scale_m %||% 30L
    
    et_r <- rgee::ee_as_raster(
      et_img,
      region = aoi_geom,
      via    = "drive",
      dsn    = tempfile(fileext = ".tif"),
      scale  = scale_et
    ) |>
      terra::rast()
    
    et_r <- terra::project(et_r, template, method = "bilinear")
    
    q_et5 <- stats::quantile(
      terra::values(et_r, mat = FALSE),
      probs = 0.05,
      na.rm = TRUE
    )
    et_screen <- terra::ifel(et_r > q_et5, et_r, NA)
  }
  
  # ------------------------------------------------------------------
  # NDVI
  # ------------------------------------------------------------------
  ndvi_cot      <- NULL
  ndvi_screen   <- NULL
  eligible_mask <- crop_mask
  
  ndvi_fp <- nayd_cache_path(cache_dir, "ndvi", county_name, year, crop_name, sensor)
  
  ndvi_season <- ndvi_season_image(
    sensor = sensor,
    start  = d1,
    end    = d2,
    geom   = aoi_geom,
    yr     = year
  )
  
  if (is.null(ndvi_season) || is.null(ndvi_season$ndvi)) {
    warning("NDVI composite missing; falling back to ET-only inside crop mask.")
    eligible_mask <- crop_mask
  } else {
    ndvi_img   <- ndvi_season$ndvi
    scale_ndvi <- ndvi_season$scale_m
    
    ndvi_r <- rgee::ee_as_raster(
      ndvi_img,
      region = aoi_geom,
      via    = "drive",
      dsn    = tempfile(fileext = ".tif"),
      scale  = scale_ndvi
    ) |>
      terra::rast()
    
    ndvi_r <- terra::project(ndvi_r, template, method = "bilinear")
    
    if (all(is.na(terra::values(ndvi_r, mat = FALSE)))) {
      warning("NDVI raster is empty after projection; using crop mask only.")
      eligible_mask <- crop_mask
    } else {
      q_nd5 <- stats::quantile(
        terra::values(ndvi_r, mat = FALSE),
        probs = 0.05,
        na.rm = TRUE
      )
      ndvi_screen <- terra::ifel(ndvi_r > q_nd5, ndvi_r, NA)
      
      elig <- .ndvi_eligibility(
        ndvi_r               = ndvi_r,
        crop_mask            = crop_mask,
        NASS_harvest_area_ha = NASS_harvest_area_ha,
        area_buffer_frac     = area_buffer_frac,
        candidate_pcts       = seq(0.90, 0.10, by = -0.01)
      )
      
      ndvi_elig_mask <- elig$mask
      
      if (!is.null(ndvi_elig_mask)) {
        eligible_mask <- terra::ifel(!is.na(ndvi_elig_mask) & !is.na(crop_mask), 1, NA)
        cat(
          sprintf(
            "  NDVI eligibility: thr=%.3f, area=%.0f ha (NASS=%s ha)\n",
            elig$threshold %||% NaN,
            elig$area_ha   %||% NaN,
            fmt_ha(NASS_harvest_area_ha)
          )
        )
      } else {
        eligible_mask <- crop_mask
        cat("NDVI eligibility failed → using crop mask only.\n")
      }
      
      n_elg <- sum(terra::values(eligible_mask, mat = FALSE) == 1, na.rm = TRUE)
      if (is.na(n_elg) || n_elg == 0) return(NULL)
      
      ndvi_cot <- terra::mask(ndvi_screen, eligible_mask)
      terra::writeRaster(
        ndvi_cot, ndvi_fp, overwrite = TRUE,
        wopt = list(datatype = "FLT4S",
                    gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
      )
    }
  }
  
  # ------------------------------------------------------------------
  # Build weights inside eligibility
  # ------------------------------------------------------------------
  if (has_et) {
    et_cot <- terra::mask(et_screen, eligible_mask)
    et_fp  <- nayd_cache_path(cache_dir, "et", county_name, year, crop_name, sensor)
    terra::writeRaster(
      et_cot, et_fp, overwrite = TRUE,
      wopt = list(datatype = "FLT4S",
                  gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
    )
    if (all(is.na(terra::values(et_cot)))) return(NULL)
  }
  
  w_ndvi <- norm_sensor(ndvi_cot, new_min = 0.01, new_max = 1)
  names(w_ndvi) <- "w_ndvi"
  
  if (has_et) {
    w_et   <- norm_sensor(et_cot, new_min = 0.01, new_max = 1)
    names(w_et) <- "w_et"
    w_norm <- a_et * w_et + a_ndvi * w_ndvi
  } else {
    w_norm <- w_ndvi
  }
  
  w_pos <- terra::ifel(!is.na(w_norm), w_norm, NA)
  s     <- as.numeric(terra::global(w_pos, "sum", na.rm = TRUE))
  
  if (!is.na(s) && s > 0) {
    w_norm_final <- w_pos / s
  } else {
    w_norm_final <- w_pos
  }
  
  names(w_norm_final) <- "w_norm"
  
  weight_fp <- nayd_cache_path(cache_dir, "weight", county_name, year, crop_name, sensor)
  terra::writeRaster(
    w_norm_final, weight_fp, overwrite = TRUE,
    wopt = list(datatype = "FLT4S",
                gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
  )
  
  cat(
    sprintf(
      "  ✓ %s %d: crop pix=%d, sum(w_norm_final)=%.6f\n",
      county_name, year, n_pix,
      as.numeric(terra::global(w_norm_final, "sum", na.rm = TRUE))
    )
  )
  
  cdl_mask_area_ha  <- res$area_ha
  elig_mask_area_ha <- .area_ha_bin(eligible_mask)
  
  ratio_cdl_to_nass  <- if (is.finite(NASS_harvest_area_ha) && NASS_harvest_area_ha > 0) {
    cdl_mask_area_ha / NASS_harvest_area_ha
  } else {
    NA_real_
  }
  ratio_elig_to_nass <- if (is.finite(NASS_harvest_area_ha) && NASS_harvest_area_ha > 0) {
    elig_mask_area_ha / NASS_harvest_area_ha
  } else {
    NA_real_
  }
  
  cat(
    sprintf(
      "  mask areas (ha): CDL=%s (=%s of NASS) | Eligible=%s (=%s of NASS)\n",
      fmt_ha(cdl_mask_area_ha),  fmt_pct(ratio_cdl_to_nass),
      fmt_ha(elig_mask_area_ha), fmt_pct(ratio_elig_to_nass)
    )
  )
  
  ndvi_sensor <- if (!is.null(ndvi_season) && !is.null(ndvi_season$ndvi)) sensor else NA_character_
  
  diag <- list(
    county             = county_name,
    year               = year,
    mode               = mode,
    cdl_strategy       = res$strategy,
    NASS_harvest_area_ha       = NASS_harvest_area_ha,
    CDL_mask_area_ha   = cdl_mask_area_ha,
    Eligible_area_ha   = elig_mask_area_ha,
    CDL_over_NASS      = ratio_cdl_to_nass,
    Eligible_over_NASS = ratio_elig_to_nass,
    NDVI_thre          = if (exists("elig")) elig$threshold else NA_real_,
    ndvi_sensor        = ndvi_sensor,
    ndvi_scale_m       = if (exists("scale_ndvi")) as.integer(scale_ndvi) else NA_integer_,
    et_scale_m         = if (exists("scale_et"))   as.integer(scale_et)   else NA_integer_,
    a_ndvi             = a_ndvi,
    a_et               = a_et,
    q05_ndvi           = if (exists("q_nd5")) as.numeric(q_nd5) else NA_real_,
    q05_et             = if (exists("q_et5")) as.numeric(q_et5) else NA_real_
  )
  
  list(
    weight      = w_norm_final,
    ndvi        = if (!is.null(ndvi_cot)) ndvi_cot else NULL,
    et          = if (!is.null(et_cot))   et_cot   else NULL,
    county_v    = county_v,
    template    = template,
    crop_mask   = crop_mask,
    eligibility = eligible_mask,
    diagnostics = diag
  )
}

