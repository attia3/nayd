
# 02_get_weights_county_year.R
#
# Purpose:
#   Build CDL-based crop mask and NDVI/ET-based weights for a single
#   county-year and crop, aligned to a template grid.
#
# Main API:
#   get_weights_for_county_year(county_sf, county_name, sensor, year,
#                               season_start, season_end,
#                               cdl_rast, CRS_TARGET,
#                               crop_code, crop_name,
#                               a_ndvi, a_et,
#                               NASS_planed_area_ha, NASS_harvest_area_ha,
#                               lower_thresh, upper_thresh,
#                               area_buffer_frac, pre2008_mask_strategy)
#
# Side effects:
#   - writes NDVI, ET, and weight rasters to cache_dir via .cache_path()

#---------ET/NDVI download z-score calculation and weight-----------------

# helper function 

norm_sensor <- function(x_cot, new_min = 0.01, new_max = 1) {
  
  # get min/max over cotton pixels
  r <- terra::global(x_cot, c("min", "max"), na.rm = TRUE)
  x_min <- as.numeric(r[1, "min"])
  x_max <- as.numeric(r[1, "max"])
  
  # no valid range → return NA (or mid value if you prefer)
  if (!is.finite(x_min) || !is.finite(x_max) || x_max <= x_min) {
    # Option 1: all NA
    x_norm <- x_cot * NA
    
    # Option 2 (alternative): everything gets mid weight
    # mid <- (new_min + new_max) / 2
    # x_norm <- classify(x_cot, cbind(-Inf, Inf, mid))
    
    names(x_norm) <- "weight"
    return(x_norm)
  }
  
  # standard min–max scaling to [0,1]
  x_s <- (x_cot - x_min) / (x_max - x_min)
  
  # stretch to [new_min, new_max]
  x_norm <- x_s * (new_max - new_min) + new_min
  names(x_norm) <- "weight"
  
  x_norm
}

norm_county <- function(x) {
  x |>
    toupper() |>
    str_replace("\\s+COUNTY$", "") |>
    str_trim()
}

moa_CDL_TX_dir <- "C:/Users/ahmed.attia/OneDrive - Texas A&M AgriLife/ML_study/SpatialData/CDL_rasters/TX_state"
# exact-year or 2008-as-reference for pre-2008
cdl_path_for_year <- function(year_eff) {
  file.path(moa_CDL_TX_dir, sprintf("CDL_%d_48.tif", year_eff))
}

# No per-county writes, just return the path to the right statewide mosaic
cdl_for_processing <- function(year) {
  # Use 2008 for any year < 2008; use 2023 for any year > 2023
  if (year < 2008L) {
    yr_eff <- 2008L
  } else {
    yr_eff <- as.integer(year)
  }
  fp <- cdl_path_for_year(yr_eff)
  message(sprintf("[cdl_for_processing] year %d -> %s", yr_eff, fp))
  if (!file.exists(fp)) return(NA_character_)
  fp
}

library("rgee","reticulate","geojsonio")  # if not already

# Create & configure a new conda env named "rgee"
#rgee::ee_install(py_env = "rgee", conda = "auto")   # downloads miniconda if needed

# where to write cached rasters
cache_dir <- file.path(
  "C:/Users/ahmed.attia/OneDrive - Texas A&M AgriLife/ML_study/SpatialData/newFiles","cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

.sanitize_county <- function(x) {
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

.cache_path <- function(kind, county, year, crop_name, sensor = NULL) {
  stopifnot(kind %in% c("ndvi", "et", "weight"))
  
  sub <- file.path(cache_dir, kind, .sanitize_county(crop_name))  # crop subdir
  dir.create(sub, recursive = TRUE, showWarnings = FALSE)
  
  sensor_tag <- if (!is.null(sensor) && nzchar(sensor)) paste0("_", sensor) else ""
  
  file.path(
    sub,
    sprintf(
      "%s_%s_%d%s.tif",
      kind,
      .sanitize_county(county),
      as.integer(year),
      sensor_tag
    )
  )
}

use_gee <- TRUE
#rgee::ee_Initialize(drive = TRUE)

rs_mode_for_year <- function(y) {
  if (y >= 2008) "NDVI_ET"
  else if (y > 1999) "NDVI_ET_2008mask"
  else "NDVI_ONLY_2008mask"   # for 1984–1999
}

# Build a crop mask for a county, aligned to a template (CDL @ 30 m)
# strategy: "y2008", "union2008_2012", "dilate", or "auto"
make_crop_mask_year <- function(template, county_v, year,
                                  NASS_planed_area_ha = NULL,
                                  lower_thresh = NULL,     # 60% threshold
                                  upper_thresh  = NULL,
                                  crop_code = NULL) {      # CDL code for crop
  stopifnot(inherits(template, "SpatRaster"),
            inherits(county_v,  "SpatVector"))
  
  # ---- helpers ----
  .to_bin <- function(r) terra::ifel(r == crop_code, 1, NA)
  .proj_crop_mask <- function(r) {
    r |>
      terra::project(template, method = "near") 
  }
  .area_ha <- function(binmask) {
    # sum cell areas (in m^2) where mask==1 (non-NA), convert to ha
    a <- terra::cellSize(binmask, unit = "m")     # area for every cell
    a <- terra::mask(a, binmask)                  # keep area only where mask has data (crop)
    as.numeric(terra::global(a, "sum", na.rm = TRUE)) / 1e4  # m² -> ha
  }
  .load_cdl_bin <- function(y) {
    fp <- cdl_path_for_year(as.integer(y))
    if (is.null(fp) || !file.exists(fp)) return(NULL)
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
    
    # pixel-wise OR (max); keep 1/NA output
    u <- do.call(terra::mosaic, c(rs, list(fun = "max")))
    terra::ifel(u == 1, 1, NA)
  }
  # 6) intersection of multiple years (pixel-wise AND)
  .intersect_bin <- function(years_vec) {
    years_vec <- unique(as.integer(years_vec))
    rs <- lapply(years_vec, .load_cdl_bin)
    rs <- Filter(Negate(is.null), rs)
    if (length(rs) == 0) return(NULL)
    if (length(rs) == 1) return(rs[[1]])
    
    # start with first mask, keep only pixels present in all
    b <- rs[[1]]
    for (k in 2:length(rs)) {
      b <- terra::ifel(!is.na(b) & !is.na(rs[[k]]), 1, NA)
    }
    b
  }
  # ---- logic ----
   used_strategy <- NA_character_
  .win3 <- function(year) {
    if (year < 2008L) 2008L:2010L else (year-1L):(year+1L)
  }
  
  # Try exact-year CDL if year >= 2008 and file exists
  mask_year <- if (year >= 2008L) .load_cdl_bin(year) else NULL
  
  if (!is.null(mask_year)) {
    if (!is.null(NASS_planed_area_ha) && is.finite(NASS_planed_area_ha)) {
      area_y <- .area_ha(mask_year)
      ratio <- area_y/NASS_planed_area_ha
      # Case 1: CDL area << NASS planted → likely lower-detection → use 3-year union
      if (!is.na(area_y) && ratio < lower_thresh) {  
        # too small vs NASS → 3-year window per rule
        yrs_win <- .win3(year)
        mask3   <- .union_bin(yrs_win)
        if (!is.null(mask3)) {
          used_strategy <- sprintf("union_%s", paste(yrs_win, collapse = "_"))
          return(list(mask = mask3, strategy = used_strategy,
                      area_ha = .area_ha(mask3)))
        }
        # fallback continues below...
      }
      # Case 2: CDL area >> NASS planted → likely over-detection → use 3-year intersection
      if (!is.na(ratio) && ratio >= upper_thresh) {
        yrs_win <- .win3(year)
        # Option A: intersection of 3 years
        mask3_int <- .intersect_bin(yrs_win)
        if (!is.null(mask3_int)) {
          return(list(
            mask     = mask3_int,
            used_strategy = sprintf("intersect_%s", paste(yrs_win, collapse = "_")),
            area_ha  = .area_ha(mask3_int)))}
        
     }
    }
    
    used_strategy <- sprintf("exact_%d", year)
    return(list(
      mask     = mask_year,
      strategy = used_strategy,
      area_ha  = .area_ha(mask_year)
    ))
  }
  
  # No exact-year CDL (or year < 2008): build 3-year union around that year.
  yrs_win <- .win3(year)
  mask3 <- .union_bin(yrs_win)
  if (!is.null(mask3)) {
    used_strategy <- sprintf("union_%s", paste(yrs_win, collapse = "_"))
    return(list(mask = mask3, strategy = used_strategy,
                area_ha = .area_ha(mask3)))
    }
  
  # Last resort
  list(mask = NULL, strategy = "failed_to_build_mask", area_ha = NA_real_)
  
}
  
fmt_ha <- function(x) {
  if (is.null(x) || is.na(x) || !is.finite(x)) return("NA")
  format(round(as.numeric(x), 1), big.mark = ",")
}

# Build NDVI eligibility threshold so eligible area ≳ NASS area (with a small buffer)
# ndvi_r: SpatRaster (band "NDVI" in [0,1]) aligned to template
# crop_mask: 1/NA SpatRaster aligned to template
.ndvi_eligibility <- function(ndvi_r, crop_mask, NASS_harvest_area_ha,
                              area_buffer_frac = NULL, # allow ~x% slack
                              candidate_pcts = seq(0.90, 0.10, by = -0.02)) {
  stopifnot(inherits(ndvi_r, "SpatRaster"), inherits(crop_mask, "SpatRaster"))
  # restrict to crop pixels
  ndv <- terra::mask(ndvi_r, crop_mask)
  vals <- terra::values(ndv, mat = FALSE)
  if (all(is.na(vals))) return(list(mask = NULL, threshold = NA_real_, area_ha = NA_real_))
  
  # areas per pixel (ha), then mask to NDVI (non-NA) so counting area is exact
  a_m2 <- terra::cellSize(ndv, unit = "m"); a_m2 <- terra::mask(a_m2, ndv)
  a_ha <- as.numeric(terra::values(a_m2, mat = FALSE)) / 1e4
  nd  <- as.numeric(vals)
  
  # if no NASS given, use a middle quantile (e.g., 0.5) for a sensible eligibility
  if (is.null(NASS_harvest_area_ha) || !is.finite(NASS_harvest_area_ha) || NASS_harvest_area_ha <= 0) {
    thr <- stats::quantile(nd, probs = 0.50, na.rm = TRUE)
    keep <- nd >= thr
    return(list(mask = terra::ifel(ndv >= thr, 1, NA),
                threshold = as.numeric(thr),
                area_ha = sum(a_ha[keep], na.rm = TRUE)))
  }
  
  target <- NASS_harvest_area_ha * (1 + area_buffer_frac)
  
  # try high → low thresholds until we cover target area
  for (p in candidate_pcts) {
    thr <- stats::quantile(nd, probs = p, na.rm = TRUE)
    keep <- nd >= thr
    A   <- sum(a_ha[keep], na.rm = TRUE)
    if (is.finite(A) && A >= target) {
      return(list(mask = terra::ifel(ndv >= thr, 1, NA),
                  threshold = as.numeric(thr),
                  area_ha = A))
    }
  }
  
  # if none covered, just take the lowest candidate as a fallback
  thr <- stats::quantile(nd, probs = utils::tail(candidate_pcts, 1), na.rm = TRUE)
  keep <- nd >= thr
  A   <- sum(a_ha[keep], na.rm = TRUE)
  list(mask = terra::ifel(ndv >= thr, 1, NA),
       threshold = as.numeric(thr),
       area_ha = A)
}


get_weights_for_county_year <- function(county_sf, county_name, 
                                        sensor = NULL,year,
                                        season_start, season_end,
                                        cdl_rast, CRS_TARGET,
                                        crop_code = NULL,crop_name = NULL, # CDL code for crop
                                        a_ndvi = NULL, a_et = NULL, 
                                        NASS_planed_area_ha = NULL,
                                        NASS_harvest_area_ha = NULL, 
                                        lower_thresh = NULL,
                                        upper_thresh = NULL,                  # minimum fraction of NASS harvested area 
                                        area_buffer_frac = NULL,              # that the single-yr CDL must reach e.g. 80%
                                        pre2008_mask_strategy = "auto") # *** "auto","y2008","union2008_2012","dilate"
{
  stopifnot(is.character(county_name), length(county_name) == 1)
  if (!use_gee) stop("This helper uses rgee. Set use_gee = TRUE and call ee_Initialize().")
  
  mode <- rs_mode_for_year(year)
  
  aoi <- sf_as_ee(sf::st_transform(county_sf, 4326)) 
  aoi_geom <- aoi$geometry()
  d1  <- ee$Date(sprintf("%d-%s", year, season_start))
  d2  <- ee$Date(sprintf("%d-%s", year, season_end))
  
  fp <- cdl_for_processing(year)
  if (is.na(fp)) { cat("Year", year, ": no CDL mosaic\n"); return(NULL) }
  
  # CDL (statewide) -> county
  cdl_TX <- terra::rast(fp)
  county_v <- county_sf |> sf::st_transform(terra::crs(cdl_TX)) |> terra::vect()
  cdl_cnty <- terra::mask(terra::crop(cdl_TX, county_v), county_v)
  cdl_cnty <- terra::project(cdl_cnty, CRS_TARGET)
  
  # TEMPLATE = CDL grid (30 m)
  template <- cdl_cnty
  if (!terra::hasValues(cdl_cnty)) return(NULL)
  
  res <- make_crop_mask_year(
    template     = template,      # your 30 m template
    county_v     = county_v,    # this county geometry
    year         = year,           # target year
    NASS_planed_area_ha = NASS_planed_area_ha, # county-year harvested area
    lower_thresh = lower_thresh,
    upper_thresh = upper_thresh,
    crop_code = crop_code 
  )
  
  mask_ha <- res$area_ha
  nass_ha <- NASS_harvest_area_ha
  nass_planted_ha <- NASS_planed_area_ha
  crop_mask <- res$mask
  
  pct_cov <- if (!is.null(nass_ha) && is.finite(nass_ha) && nass_ha > 0 &&
                 !is.null(mask_ha) && is.finite(mask_ha)) 100 * mask_ha / nass_ha else NA_real_
  
  cat(sprintf("✓ Built via %s for county %s, year %d — NASS planted area = %s ha, crop mask = %s ha%s\n",
              res$strategy,
              county_name,
              year,
              fmt_ha(nass_planted_ha),
              fmt_ha(mask_ha),
              if (is.finite(pct_cov)) sprintf(" (%.1f%% of NASS)", pct_cov) else ""))
  
  
  if (is.null(res$mask)) {
    cat("No mask:", res$strategy, "\n")
  } else {
    cat("Built via", res$strategy, "area =", round(res$area_ha, 1), "ha\n")
  }
  
  # quick pixel guard
  n_pix <- sum(terra::values(crop_mask, mat = FALSE) == 1, na.rm = TRUE)
  if (is.na(n_pix) || n_pix == 0) return(NULL)
  
  has_et <- (mode != "NDVI_ONLY_2008mask")
  if (has_et) {
  # --- OpenET eeMETRIC ET (30 m) instead of GRIDMET ---
  et_obj <- eemetric_et_season_image(
    start = ee$Date(sprintf("%d-%s", year, season_start)),
    end   = ee$Date(sprintf("%d-%s", year, season_end)),
    geom  = aoi_geom
  )
  et_img   <- et_obj$et
  scale_et <- et_obj$scale_m %||% 30L   # tiny helper: %||% = default if NULL
  
  et_fp <- .cache_path("et", county_name, year, crop_name,sensor)
  
  et_r <- ee_as_raster(
    et_img,
    region = aoi_geom,
    via    = "drive",
    dsn    = tempfile(fileext = ".tif"),
    scale  = scale_et
  ) |>
    terra::rast()
  
  # project to CDL grid
  et_r <- terra::project(et_r, template, method = "bilinear")
  
  # light lower-tail screen, as before
  q_et5 <- stats::quantile(
    terra::values(et_r, mat = FALSE),
    probs = 0.05,
    na.rm = TRUE
  )
  et_screen <- terra::ifel(et_r > q_et5, et_r, NA)
  }
  # NDVI 
  ndvi_cot    <- NULL
  ndvi_screen <- NULL          # initialize so we can test !is.null later
  eligible_mask <- crop_mask # default fallback

  ndvi_fp <- .cache_path("ndvi", county_name, year,crop_name,sensor)
  ndvi_season <- ndvi_season_image(sensor=sensor,start=d1,end=d2,geom= aoi_geom, yr = year)
  if (is.null(ndvi_season) || is.null(ndvi_season$ndvi)) {
      warning("NDVI composite missing; falling back to ET-only inside crop mask.")
      eligible_mask <- crop_mask
    } else {
      # Use the image and its native scale
      ndvi_img   <- ndvi_season$ndvi
      scale_ndvi <- ndvi_season$scale_m 
      ndvi_r <- ee_as_raster(
        ndvi_img, region = aoi_geom, via = "drive",
        dsn = tempfile(fileext = ".tif"), scale = scale_ndvi) |> terra::rast()
      ndvi_r <- terra::project(ndvi_r, template, method = "bilinear")
    
    # If NDVI is all NA, fall back
    if (all(is.na(terra::values(ndvi_r, mat = FALSE)))) {
      warning("NDVI raster is empty after projection; using crop mask only.")
      eligible_mask <- crop_mask
    } else {
      
      q_nd5  <- stats::quantile(terra::values(ndvi_r, mat = FALSE), probs = 0.05, na.rm = TRUE)
      ndvi_screen <- terra::ifel(ndvi_r > q_nd5, ndvi_r, NA)
    
    # ----- 3) NDVI eligibility inside crop (modes with NDVI) -----

      elig <- .ndvi_eligibility(
        ndvi_r      = ndvi_r,
        crop_mask = crop_mask,
        NASS_harvest_area_ha= NASS_harvest_area_ha,
        area_buffer_frac = area_buffer_frac,
        candidate_pcts   = seq(0.90, 0.10, by = -0.01)
      )
    
      ndvi_elig_mask <- elig$mask
      
      if (!is.null(ndvi_elig_mask)) {
        # eligibility = crop ∩ NDVI-eligible
        eligible_mask <- terra::ifel(!is.na(ndvi_elig_mask) & !is.na(crop_mask), 1, NA)
        cat(sprintf("  NDVI eligibility: thr=%.3f, area=%.0f ha (NASS=%s ha)\n",
                    elig$threshold %||% NaN,
                    elig$area_ha %||% NaN,
                    fmt_ha(NASS_harvest_area_ha)))
      } else {
        # fallback to crop only if eligibility failed
        eligible_mask <- crop_mask
        cat("NDVI eligibility failed → using crop mask only.\n")
      }

    # quick guard on eligibility
    n_elg <- sum(terra::values(eligible_mask, mat = FALSE) == 1, na.rm = TRUE)
    if (is.na(n_elg) || n_elg == 0) return(NULL)
    
    # *** now mask NDVI using the (already built) crop_mask
    ndvi_cot <- terra::mask(ndvi_screen, eligible_mask)
    terra::writeRaster(ndvi_cot, ndvi_fp, overwrite = TRUE,
                       wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW","PREDICTOR=2")))
    }
   }
  # ----- 4) Build weights INSIDE eligibility -----

  # ET on crop pixels (use NDVI screen when available)
  et_cot <- NULL
  et_cot <- if (has_et) {
    terra::mask(et_screen, ndvi_cot)}
  
  if(has_et){
    terra::writeRaster(et_cot, et_fp, overwrite = TRUE,
                       wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW","PREDICTOR=2")))
    if (all(is.na(terra::values(et_cot)))) return(NULL)
  }

# --- re-scaling and mix
  # ndvi_cot might be NULL if NDVI failed

  w <- NULL  # final weight raster
  
  if (has_et) {
    w_et   <- norm_sensor(et_cot, new_min = 0.01, new_max = 1)
    names(w_et) <- "w_et"
  }
  
  w_ndvi <- norm_sensor(ndvi_cot, new_min = 0.01, new_max = 1)
  names(w_ndvi) <- "w_ndvi"
  
  w_norm <- w_ndvi * NA
  
  if(has_et){
    
    w_norm <- (a_et * w_et + a_ndvi * w_ndvi)    # assuming a_et + a_ndvi = 1
    
  }else{
    
    w_norm <- w_ndvi
    
  }
  
  # renormalize so sum of w_norm over cotton = 1
  w_pos <- ifel(!is.na(w_norm), w_norm, NA)
  s <- as.numeric(global(w_pos, "sum", na.rm = TRUE))
  
  if (!is.na(s) && s > 0) {
    w_norm_final <- w_pos / s
  } else {
    w_norm_final <- w_pos  # fallback
  }
  
  names(w_norm_final) <- "w_norm"
  
  weight_fp <- .cache_path("weight", county_name, year,crop_name,sensor)
  terra::writeRaster(
    w_norm_final, weight_fp, overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW","PREDICTOR=2"))
  )
  
  # quick debug print
  cat(sprintf("  ✓ %s %d: crop pix=%d, sum(w_norm_final)=%.6f\n",
              county_name, year, n_pix, as.numeric(terra::global(w_norm_final, "sum", na.rm = TRUE))))
  
#  --- some diagnosis---
  
  # --- helpers (top of file once) ---
  .area_ha_bin <- function(binmask) {
    if (is.null(binmask)) return(NA_real_)
    a <- terra::cellSize(binmask, unit = "m")
    a <- terra::mask(a, binmask)
    as.numeric(terra::global(a, "sum", na.rm = TRUE)) / 1e4  # m² -> ha
  }
  fmt_pct <- function(x) ifelse(is.finite(x), sprintf("%.1f%%", 100*x), "NA")
  
  # --- inside get_weights_for_county_year(), after you have: 
  # crop_mask, eligible_mask, NASS_harvest_area_ha, ndvi_cot (maybe),
  # ndvi_img/scale_ndvi (maybe), et_cot, w_norm, and strings like `mode`
  # and the CDL mask build 'strategy' from make_crop_mask_year() call (res$strategy) ---
  
  cdl_mask_area_ha  <- res$area_ha                     # from make_crop_mask_year()
  elig_mask_area_ha <- .area_ha_bin(eligible_mask)
  
  # Ratios vs NASS (harvested)
  ratio_cdl_to_nass  <- if (is.finite(NASS_harvest_area_ha) && NASS_harvest_area_ha > 0) cdl_mask_area_ha / NASS_harvest_area_ha else NA_real_
  ratio_elig_to_nass <- if (is.finite(NASS_harvest_area_ha) && NASS_harvest_area_ha > 0) elig_mask_area_ha / NASS_harvest_area_ha else NA_real_
  
  # Simple console line
  cat(sprintf(
    "  mask areas (ha): CDL=%s (=%s of NASS) | Eligible=%s (=%s of NASS)\n",
    fmt_ha(cdl_mask_area_ha),  fmt_pct(ratio_cdl_to_nass),
    fmt_ha(elig_mask_area_ha), fmt_pct(ratio_elig_to_nass)
  ))
  
  # Optional: quick sensor meta for later comparisons
  ndvi_sensor <- if (exists("ndvi_img")) {
    # your caller sets `sensor` already ("MODIS","LANDSAT","SENTINEL2")
    sensor
  } else {
    NA_character_
  }
  
  # Optional: stash screening thresholds used
  diag <- list(
    county             = county_name,
    year               = year,
    mode               = mode,                 # "NDVI_ET", "NDVI_ET_2008mask", "ET_ONLY_2008mask"
    cdl_strategy       = res$strategy,         # e.g., "exact_2015" or "union_2014_2015_2016"
    NASS_harvest_area_ha       = NASS_harvest_area_ha,
    CDL_mask_area_ha   = cdl_mask_area_ha,
    Eligible_area_ha   = elig_mask_area_ha,
    CDL_over_NASS      = ratio_cdl_to_nass,
    Eligible_over_NASS = ratio_elig_to_nass,
    NDVI_thre = elig$threshold,
    ndvi_sensor        = ndvi_sensor,          # "MODIS"/"LANDSAT"/"SENTINEL2" or NA
    ndvi_scale_m       = if (exists("scale_ndvi")) as.integer(scale_ndvi) else NA_integer_,
    et_scale_m         = if (exists("scale_et")) as.integer(scale_et) else NA_integer_,                 # as used
    a_ndvi             = a_ndvi,
    a_et               = a_et,
    q05_ndvi           = if (exists("q_nd5"))  as.numeric(q_nd5)  else NA_real_,
    q05_et             = if (exists("q_et5"))  as.numeric(q_et5)  else NA_real_
  )
  
  # (Optionally) append one-line CSV per run for later analysis
  # readr::write_csv(as.data.frame(diag), "weights_diag_log.csv", append = file.exists("weights_diag_log.csv"))
  
  # --- final return updated ---
  return(list(
    weight        = w_norm_final,
    ndvi          = if (!is.null(ndvi_cot)) ndvi_cot else NULL,
    et            = if (!is.null(et_cot)) et_cot else NULL,
    county_v = county_v,
    template      = template,
    crop_mask   = crop_mask,
    eligibility   = eligible_mask,
    diagnostics   = diag
  ))
  
}

# test_gets_weight <- get_weights_for_county_year(county_sf, county_name, sensor = "LANDSAT",year,
#                                                 season_start, season_end,
#                                                 cdl_rast = cdl_THP, CRS_TARGET,crop_code = crop_code, # CDL code for cotton
#                                                 a_ndvi = 0.5, a_et = 0.5,
#                                                 NASS_harvest_area_ha = NASS_harvest_area_ha,
#                                                 thresh = 0.80,
#                                                 area_buffer_frac = 0.03,
#                                                 pre2008_mask_strategy = "auto")
