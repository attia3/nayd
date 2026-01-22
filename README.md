NASS-anchored-Yield-Disaggregation-Workflow (NAYD)
This repository provides an R workflow to disaggregate official, coarse-scale yield statistics (e.g. USDA NASS county yields) to field-scale “harvest units” using a combination of: land-cover / crop-type rasters (e.g. USDA Cropland Data Layer), remote-sensing metrics (NDVI, ET) from Google Earth Engine, and reported area and yield for each administrative unit. The code can be applied to any crop and region where equivalent inputs are available. In the documentation and examples we demonstrate the workflow using winter wheat and cotton in Texas, but the methods are not limited to these crops or locations.

## Features
Builds county- and crop-specific masks from USDA Cropland Data Layer (CDL) with consistency checks against NASS planted/harvested areas.
Derives seasonal NDVI and ET composites from Landsat / OpenET.
Computes pixel-level weights combining NDVI and ET.
Segments fields and large pseudo-fields into harvest units with area constraints.
Disaggregates NASS county yield across units using a beta-weighting scheme while preserving county means.
Clusters harvest units into ~5–10 km production zones.
Validates segmented yields against independent variety-testing trials.
Repository structure

R/01_ndvi_et_season.R – NDVI/ET seasonal composites (rgee/OpenET).

R/02_build_spatial_weights.R – CDL mask + NDVI/ET weight generator.

R/03_segment_fields_pseudo_units.R – field segmentation and tiling.

R/04_disagg_units_yield_safe.R – NASS-consistent yield assignment to units.

R/05_cluster_units_5km.R – cluster units into production zones.

R/06_process_county_year.R – wrapper for one county-year.

R/07_compute_weights_batch.R – batch weights for many counties/years.

R/08_run_disagg_clustering.R – main loop + segmentation + export.

## Requirements
R >= 4.3
Packages: sf, terra, dplyr, data.table, ggplot2, rgee, reticulate, readxl, stringr, purrr, viridis, patchwork, etc.
Google Earth Engine account and rgee configured.
Access to:
Texas county boundary shapefile (EPSG:4326)
CDL mosaics (30 m) for 2008–present, in EPSG:5070
NASS county-level planted/harvested area and yield tables
## Main features
Crop-specific masks from CDL

Builds county- and crop-specific masks from CDL
Applies consistency checks against NASS planted / harvested area
Handles pre-2008 years via multi-year unions / intersections when CDL is noisy
Remote-sensing-based weights

Derives seasonal NDVI and ET composites from Landsat / Sentinel-2 via Google Earth Engine / OpenET
Normalizes NDVI and ET to [0, 1] and combines them into a per-pixel weight surface
Field segmentation / harvest units

Segments the crop mask into contiguous “fields”
Splits overly large fields into pseudo-fields using a regular grid
Enforces minimum / maximum field-size constraints (e.g. 5–1000 ha)
NASS-consistent yield disaggregation

Scales field weights with a β-exponent to control contrast between high- and low-weight pixels
Allocates NASS county yield to units in a way that:
preserves the county mean yield (within numerical noise)
respects crop-specific min / max yield caps
avoids pathological extreme yields
Clustering into production zones

Clusters nearby harvest units (e.g. ~5–10 km radius) into production zones
Returns cluster-level mean yield and area for mapping and downstream modelling
Validation against variety trials

Matches production-zone yields to independent variety-trial plots
Supports evaluation across multiple radii and β values
Repository structure
A typical workflow uses the scripts in R/ in this order:

R/01_ndvi_et_season.R
Build seasonal NDVI and ET composites (via rgee + Google Earth Engine / OpenET).

R/02_build_spatial_weights.R

Construct CDL masks
Apply NASS area consistency checks
Compute NDVI+ET pixel weights and normalize to 1 within each county / crop / year.
R/03_segment_fields_pseudo_units.R
Segment the weight raster into fields and pseudo-fields and compute per-unit weight and area.

R/04_disagg_units_yield_safe.R
Allocate NASS county yield to units with:

β-weighting
crop-specific min / max yield caps
a final recentering step to recover the county mean.
R/05_cluster_units_prod_zones.R
Cluster harvest units into production zones (e.g. 5–10 km, complete-linkage).

R/06_process_county_year.R
Wrapper that:

takes a pre-computed weight raster for one county-year,
segments to units,
disaggregates yield,
clusters units,
returns diagnostic tables and plots.
R/07_compute_weights_batch.R
Convenience loop to build weights across many counties / years, caching intermediate rasters.

R/08_run_disagg_clustering.R
Main driver script: loops over counties × years, calls process_county_year(), and writes out cluster-level yield tables + diagnostic figures.

## Demo data (no-Earth-Engine example)
The repository ships a small demo dataset under inst/extdata/ to allow testing the workflow without running Google Earth Engine:

inst/extdata/nass_demo.csv
Small NASS-style table with county, year, crop, harvested area (ha) and yield (kg/ha) for a few Texas counties.

inst/extdata/tx_counties_demo.gpkg
A simplified Texas counties layer in EPSG:4326, containing only the counties used in the demo (e.g. Deaf Smith, Gaines).

inst/extdata/WH_weight_DeafSmith_16.tif
Pre-computed normalised weight raster for winter wheat in Deaf Smith County, TX, for a representative season.
Each pixel value represents the relative contribution to county-level wheat production (sum = 1 over eligible pixels).

inst/extdata/CO_weight_Gaines_21.tif
Pre-computed weight raster for cotton in Gaines County, TX, for a representative season.

inst/extdata/CDL_DeafSmith_16.tif
A CDL excerpt for Deaf Smith County (crop codes and attributes as in the original CDL). Used in the demo to illustrate the crop mask and weight surface.

inst/extdata/WH_dis-yield.pdf
Example figure: disaggregated winter-wheat yield in Deaf Smith County with histogram and NASS comparison.

inst/extdata/CO_dis-yield.pdf
Example figure: disaggregated cotton yield in Gaines County with histogram and NASS comparison.

These files are only for illustration and unit tests; they are not intended as an authoritative dataset.
## Quick start

```r
# 0. Install dependencies
# install.packages(c("sf","terra","dplyr","rgee", ...))

# 1. Initialize Earth Engine
# Optional: manually point reticulate to the Python env used by rgee
Sys.setenv(
  RETICULATE_PYTHON = "~/AppData/Local/r-miniconda/envs/rgee_py311/python.exe"
)
library(reticulate)
py_config()

library(rgee)
# This installs a miniconda env and earthengine-api into it
rgee::ee_install()
# One-time authentication (will open browser)
rgee::ee_Authenticate()
# Initialize in your R session (usually once per session)
rgee::ee_Initialize(drive = TRUE)

# 2. Prepare NASS tables, county shapefile, and CDL mosaics
# (see data/README_data.md)

# 3. Compute weights (optional, precompute)
source("R/02_build_spatial_weights.R")
source("R/07_compute_weights_batch.R")
compute_weights_batch(...)

# 4. Run disaggregation
source("R/03_segment_fields_pseudo_units.R")
source("R/04_disagg_units_yield_safe.R")
source("R/05_cluster_units_prod_zones.R")
source("R/06_process_county_year.R")
source("R/08_run_disagg_clustering.R")

library(NAYD)

nass_cotton <- read.csv(system.file("extdata", "nass_demo_cotton.csv", package = "NAYD"))
tx_cnties   <- sf::st_read(system.file("extdata", "tx_counties_demo.gpkg", package = "NAYD"))

compute_weights_batch(
  nass_df     = nass_cotton,
  years       = 2018:2020,
  counties_sf = tx_cnties,
  crop_name   = "Cotton",
  crop_code   = 2,
  sensor      = "Landsat",
  season_start = "06-01",
  season_end   = "09-30"
)
