utils::globalVariables(c(
  "patches", "area", "w_norm", "area_ha", "w_sum"
))
#' @importFrom rlang .data
NULL

#' Normalize county names
#'
#' Utility to normalize county names so that joins between NASS tables,
#' shapefiles, and validation data are consistent.
#'
#' Typical usage is to lower case, trim whitespace, and collapse
#' punctuation / multiple spaces.
#'
#' @param x Character vector of county names.
#'
#' @return Character vector of normalized county names.
#' @export
norm_county <- function(x) {
  x |>
    toupper() |>
    str_replace("\\s+COUNTY$", "") |>
    str_trim()
}

#' @keywords internal
norm_sensor <- function(x_cot, new_min = 0.01, new_max = 1) {
  r <- terra::global(x_cot, c("min", "max"), na.rm = TRUE)
  x_min <- as.numeric(r[1, "min"])
  x_max <- as.numeric(r[1, "max"])
  
  if (!is.finite(x_min) || !is.finite(x_max) || x_max <= x_min) {
    x_norm <- x_cot * NA
    names(x_norm) <- "weight"
    return(x_norm)
  }
  
  x_s <- (x_cot - x_min) / (x_max - x_min)
  x_norm <- x_s * (new_max - new_min) + new_min
  names(x_norm) <- "weight"
  x_norm
}

#' @keywords internal
fmt_ha <- function(x) {
  if (is.null(x) || is.na(x) || !is.finite(x)) return("NA")
  format(round(as.numeric(x), 1), big.mark = ",")
}

#' @keywords internal
fmt_pct <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%.1f%%", 100 * x)
}

#' @keywords internal
.area_ha_bin <- function(binmask) {
  if (is.null(binmask)) return(NA_real_)
  a <- terra::cellSize(binmask, unit = "m")
  a <- terra::mask(a, binmask)
  as.numeric(terra::global(a, "sum", na.rm = TRUE)) / 1e4
}
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
#' @keywords internal
.sanitize_id <- function(x) {
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

#' CDL mosaic path for a given year
#'
#' Returns the expected CDL file path (for Texas, state FIPS 48) and
#' optionally checks if it exists on disk.
#'
#' For years before 2008, it returns the 2008 mosaic path, which is then
#' typically used as a static mask.
#'
#' @param year Integer year.
#' @param cdl_dir Directory that contains CDL mosaics named
#'   like \code{"CDL_<year>_48.tif"}.
#' @param check_exists Logical; if TRUE, returns \code{NA_character_}
#'   when the file does not exist.
#'
#' @return Character file path, or \code{NA_character_} if
#'   \code{check_exists = TRUE} and the file is missing.
#' @export
cdl_for_processing <- function(year,
                               cdl_dir      = getOption("NAYD.cdl_dir", "."),
                               check_exists = TRUE) {
  y_eff <- if (year < 2008L) 2008L else as.integer(year)
  fp <- file.path(cdl_dir, sprintf("CDL_%d_48.tif", y_eff))
  if (check_exists && !file.exists(fp)) {
    return(NA_character_)
  }
  fp
}
#' @keywords internal
rs_mode_for_year <- function(y) {
  if (y >= 2008) "NDVI_ET"
  else if (y > 1999) "NDVI_ET_2008mask"
  else "NDVI_ONLY_2008mask"
}

#' Internal helper: build cache paths for NAYD rasters
#'
#' Construct a standardized path in \code{cache_dir/kind/CROP_NAME/}
#' for NDVI, ET, or weight rasters.
#'
#' @param cache_dir Directory where cache files are stored.
#' @param kind One of \code{"ndvi"}, \code{"et"}, \code{"weight"}.
#' @param county County name (character).
#' @param year Integer year.
#' @param crop_name Character crop name.
#' @param sensor Optional sensor tag (e.g. \code{"LANDSAT"}).
#'
#' @return A character string with the full path to the cache file.
#' @keywords internal
nayd_cache_path <- function(cache_dir, kind, county, year, crop_name, sensor = NULL) {
  stopifnot(kind %in% c("ndvi", "et", "weight"))
  
  subdir <- file.path(cache_dir, kind, .sanitize_id(crop_name))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  sensor_tag <- if (!is.null(sensor) && nzchar(sensor)) paste0("_", sensor) else ""
  
  file.path(
    subdir,
    sprintf(
      "%s_%s_%d%s.tif",
      kind,
      .sanitize_id(county),
      as.integer(year),
      sensor_tag
    )
  )
}

#' @rdname nayd_cache_path
#' @keywords internal
.cache_path <- function(kind, county, year, crop_name, sensor = NULL,
                        cache_dir = getOption("NAYD.cache_dir", "cache")) {
  nayd_cache_path(cache_dir, kind, county, year, crop_name, sensor)
}

# Global variable declarations to satisfy R CMD check
# (used in dplyr pipelines and rgee helpers)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".data",
    "patches",
    "area",
    "area_ha",
    "w_norm",
    "w_sum",
    "ee",
    "ee_utils_pyfunc"
  ))
}