
# 01_ndvi_et_season.R
# 
# Purpose:
#   Helper functions to build seasonal NDVI and ET images (Google Earth Engine),
#   for a given county-year and crop season.
#
# Key outputs:
#   - eemetric_et_season_image(): list(et = ee$Image, scale_m = numeric)
#   - ndvi_season_image(): list(ndvi = ee$Image, scale_m = numeric)
#
# Dependencies:
#   - rgee, sf, dplyr, terra
#
# Notes:
#   - Requires ee_Initialize() and valid Earth Engine + OpenET access.



# if (!requireNamespace("rgee", quietly = TRUE)) {
#   stop("NAYD::<fn>() requires 'rgee'. Install it first.", call. = FALSE)
# }


# --- Cloud/QA masks ----
mask_landsat_c2_sr <- function(img) {
  qa <- img$select("QA_PIXEL")
  # bits: 3=cloud shadow, 4=snow, 5=cloud, 7=cirrus
  cloudShadowBitMask <- bitwShiftL(1, 3)
  snowBitMask        <- bitwShiftL(1, 4)
  cloudBitMask       <- bitwShiftL(1, 5)
  cirrusBitMask      <- bitwShiftL(1, 7)
  
  mask <- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    And(qa$bitwiseAnd(snowBitMask)$eq(0))$
    And(qa$bitwiseAnd(cloudBitMask)$eq(0))$
    And(qa$bitwiseAnd(cirrusBitMask)$eq(0))
  
  # scale factors for C2 SR:
  # reflectance = DN*0.0000275 - 0.2  (per USGS)
  img$updateMask(mask)$
    select(c("SR_B.*"))$
    multiply(0.0000275)$add(-0.2)
}

mask_s2_sr <- function(img) {
  qa  <- img$select("QA60")
  # bits: 10=clouds, 11=cirrus
  cloudBitMask  <- bitwShiftL(1, 10)
  cirrusBitMask <- bitwShiftL(1, 11)
  mask <- qa$bitwiseAnd(cloudBitMask)$eq(0)$
    And(qa$bitwiseAnd(cirrusBitMask)$eq(0))
  # Scale: S2 SR is in 1e-4 reflectance
  img$updateMask(mask)$
    select(c("B4","B8","B5","B6","B7","B8A"))$
    multiply(1e-4)
}

# --- NDVI calculators (all return an Image with band named "ndvi" in [0,1]) ---
landsat_ndvi <- function(start, end, geom) {
  # Merge L5/L7/L8/L9 SR collections (C2 L2)
  col_names <- c(
    "LANDSAT/LT05/C02/T1_L2",
    "LANDSAT/LE07/C02/T1_L2",
    "LANDSAT/LC08/C02/T1_L2",
    "LANDSAT/LC09/C02/T1_L2"
  )
  col <- ee$ImageCollection(col_names[[1]])
  for (i in 2:length(col_names)) col <- col$merge(ee$ImageCollection(col_names[[i]]))
  col <- col$
    filterDate(start, end)$
    filterBounds(geom)$
    map(ee_utils_pyfunc(mask_landsat_c2_sr))$
    map(ee_utils_pyfunc(function(img) {
      # Red=SR_B4 (L8/9) or SR_B3 (L5/7); NIR=SR_B5 (L8/9) or SR_B4 (L5/7)
      # We harmonize by renaming if present:
      # Try modern band names first (L8/9)
      red <- img$select("SR_B4")$rename("red")
      nir <- img$select("SR_B5")$rename("nir")
      # If those are empty, fall back to L5/7
      red <- ee$Algorithms$If(red$bandNames()$size()$gt(0), red, img$select("SR_B3")$rename("red"))
      nir <- ee$Algorithms$If(nir$bandNames()$size()$gt(0), nir, img$select("SR_B4")$rename("nir"))
      img2 <- ee$Image$cat(list(red, nir))$toFloat()
      ndvi <- img2$normalizedDifference(c("nir","red"))$rename("ndvi")
      ndvi$copyProperties(img2, img2$propertyNames())
    }))
  # Seasonal composite (p90 is robust)
  col$select("ndvi")$reduce(ee$Reducer$percentile(list(90)))$rename("ndvi")$
    clip(geom)$
    max(0)$min(1)
}

s2_ndvi <- function(start, end, geom) {
  col <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(start, end)$
    filterBounds(geom)$
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 70))$ # generous; masking refines it
    map(ee_utils_pyfunc(mask_s2_sr))$
    map(ee_utils_pyfunc(function(img) {
      ndvi <- img$normalizedDifference(c("B8","B4"))$rename("ndvi")
      ndvi$copyProperties(img, img$propertyNames())
    }))$
    select("ndvi")
  
  has_data <- col$size()$gt(0)
  
  ee$Image(
    ee$Algorithms$If(
      has_data,
      col$reduce(ee$Reducer$percentile(list(90)))$
        rename("ndvi")$
        clip(geom)$
        max(0)$min(1),
      ee$Image$constant(0)$toFloat()$rename("ndvi")$clip(geom)
    )
  )
}

# choose a source per year
get_ndvi_ic <- function(start, end, geom, yr) {
  if (yr >= 2023) {
    # Aqua MODIS v6.1 (250 m)
    ee$ImageCollection("MODIS/061/MYD13Q1")$
      filterDate(start, end)$filterBounds(geom)$select("NDVI")
  } else {
    # Terra MODIS v6.1 (250 m)
    ee$ImageCollection("MODIS/061/MOD13Q1")$
      filterDate(start, end)$filterBounds(geom)$select("NDVI")
  }
}

modis_ndvi <- function(start, end, geom, yr) {
  # MOD13Q1 Terra 16-day, 250 m. NDVI scaled by 1e-4.
  ic <- get_ndvi_ic(start, end, geom, yr)
  # guard against empty collections
  has_data <- ee$Number(ic$size())$gt(0)
  ee$Image(
    ee$Algorithms$If(
      has_data,
      ic$median()$rename("NDVI")$multiply(0.0001),  # MODIS NDVI scale factor
      ee$Image(0)$toFloat()$rename("NDVI")          # fallback to avoid 0-band errors
    )
  )
}

# --- Dispatcher that returns (ndvi_img, reduce_scale_m) ---
ndvi_season_image <- function(sensor = c("MODIS","LANDSAT","SENTINEL2"),
                              start, end, geom, yr = NULL) {
  sensor <- match.arg(sensor)
  if (sensor == "MODIS") {
    list(
      ndvi = modis_ndvi(start, end, geom,yr),
      scale_m = 250
    )
  } else if (sensor == "LANDSAT") {
    list(
      ndvi = landsat_ndvi(start, end, geom),
      scale_m = 30
    )
  } else {
    list(
      ndvi = s2_ndvi(start, end, geom),
      scale_m = 10
    )
  }
}


# Seasonal OpenET eeMETRIC ET (mm) at 30 m
eemetric_et_season_image <- function(start, end, geom) {
  # start, end: ee$Date
  # geom: ee$Geometry (AOI)
  
  col <- ee$ImageCollection("OpenET/EEMETRIC/CONUS/GRIDMET/MONTHLY/v2_0")$
    filterDate(start, end)$
    filterBounds(geom)$
    select("et")$
    map(ee_utils_pyfunc(function(img) {
      # Apply the per-image scale factor stored in image properties
      scale <- ee$Number(img$get("scale_factor_et"))
      img$select("et")$
        multiply(scale)$        # scale ET to mm
        rename("et")$
        copyProperties(img, img$propertyNames())
    }))
  
  # Sum monthly ET over the season
  et_season <- col$sum()$rename("ET")
  
  # eeMETRIC is 30 m
  list(
    et       = et_season,
    scale_m  = 30L
  )
}
