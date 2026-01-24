NASS-anchored-Yield-Disaggregation-Workflow (NAYD)

library(remotes)
or
library(devtools)

install_github("attia3/NASS-anchored-Yield-Disaggregation-Workflow-NAYD")

This repository provides an R workflow to disaggregate official, coarse-scale yield statistics (e.g. USDA NASS county yields) to field-scale “harvest units” using a combination of: land-cover / crop-type rasters (e.g. USDA Cropland Data Layer), remote-sensing metrics (NDVI, ET) from Google Earth Engine, and reported area and yield for each administrative unit. The code can be applied to any crop and region where equivalent inputs are available. In the documentation and examples we demonstrate the workflow using winter wheat and cotton in Texas, but the methods are not limited to these crops or locations.

\## Features
Builds county- and crop-specific masks from USDA Cropland Data Layer (CDL) with consistency checks against NASS planted/harvested areas.

Derives seasonal NDVI and ET composites from Landsat / OpenET.

Computes pixel-level weights combining NDVI and ET.

Segments fields and large pseudo-fields into harvest units with area constraints.

Disaggregates NASS county yield across units using a beta-weighting scheme while preserving county means.

Clusters harvest units into ~5–10 km production zones.

Validates segmented yields against independent variety-testing trials.


Features

Builds county- and crop-specific masks from USDA Cropland Data Layer (CDL) with consistency checks against NASS planted/harvested areas. Derives seasonal NDVI and ET composites from Landsat / OpenET. Computes pixel-level weights combining NDVI and ET. Segments fields and large pseudo-fields into harvest units with area constraints. Disaggregates NASS county yield across units using a beta-weighting scheme while preserving county means. Clusters harvest units into ~5–10 km production zones. Validates segmented yields against independent variety-testing trials. Repository structure


R/01\_ndvi\_et\_season.R – NDVI/ET seasonal composites (rgee/OpenET).


R/02\_build\_spatial\_weights.R – CDL mask + NDVI/ET weight generator.


R/03\_segment\_fields\_pseudo\_units.R – field segmentation and tiling.

R/04\_disagg\_units\_yield\_safe.R – NASS-consistent yield assignment to units.

R/05\_cluster\_units\_5km.R – cluster units into production zones.

R/06\_process\_county\_year.R – wrapper for one county-year.

R/07\_compute\_weights\_batch.R – batch weights for many counties/years.

R/08\_run\_disagg\_clustering.R – main loop + segmentation + export.

\## Requirements

R >= 4.3

Packages: sf, terra, dplyr, data.table, ggplot2, rgee, reticulate, readxl, stringr, purrr, viridis, patchwork, etc.

Google Earth Engine account and rgee configured.

Access to:

Texas county boundary shapefile (EPSG:4326)

CDL mosaics (30 m) for 2008–present, in EPSG:5070

NASS county-level planted/harvested area and yield tables

\## Main features

Requirements

R >= 4.3 Packages: sf, terra, dplyr, data.table, ggplot2, rgee, reticulate, readxl, stringr, purrr, viridis, patchwork, etc. Google Earth Engine account and rgee configured. Access to: Texas county boundary shapefile (EPSG:4326) CDL mosaics (30 m) for 2008–present, in EPSG:5070 NASS county-level planted/harvested area and yield tables

Main features

Crop-specific masks from CDL

Builds county- and crop-specific masks from CDL Applies consistency checks against NASS planted / harvested area Handles pre-2008 years via multi-year unions / intersections when CDL is noisy Remote-sensing-based weights

Derives seasonal NDVI and ET composites from Landsat / Sentinel-2 via Google Earth Engine / OpenET Normalizes NDVI and ET to \[0, 1] and combines them into a per-pixel weight surface Field segmentation / harvest units


Segments the crop mask into contiguous “fields” Splits overly large fields into pseudo-fields using a regular grid Enforces minimum / maximum field-size constraints (e.g. 5–1000 ha) NASS-consistent yield disaggregation


Scales field weights with a β-exponent to control contrast between high- and low-weight pixels Allocates NASS county yield to units in a way that: preserves the county mean yield (within numerical noise) respects crop-specific min / max yield caps avoids pathological extreme yields Clustering into production zones

Clusters nearby harvest units (e.g. ~5–10 km radius) into production zones Returns cluster-level mean yield and area for mapping and downstream modelling Validation against variety trials

Matches production-zone yields to independent variety-trial plots Supports evaluation across multiple radii and β values Repository structure A typical workflow uses the scripts in R/ in this order:

R/01\_ndvi\_et\_season.R Build seasonal NDVI and ET composites (via rgee + Google Earth Engine / OpenET).


R/02\_build\_spatial\_weights.R


Construct CDL masks Apply NASS area consistency checks Compute NDVI+ET pixel weights and normalize to 1 within each county / crop / year. R/03\_segment\_fields\_pseudo\_units.R Segment the weight raster into fields and pseudo-fields and compute per-unit weight and area.


R/04\_disagg\_units\_yield\_safe.R Allocate NASS county yield to units with:

β-weighting crop-specific min / max yield caps a final recentering step to recover the county mean. R/05\_cluster\_units\_prod\_zones.R Cluster harvest units into production zones (e.g. 5–10 km, complete-linkage).

R/06\_process\_county\_year.R Wrapper that:

takes a pre-computed weight raster for one county-year, segments to units, disaggregates yield, clusters units, returns diagnostic tables and plots. R/07\_compute\_weights\_batch.R Convenience loop to build weights across many counties / years, caching intermediate rasters.

R/08\_run\_disagg\_clustering.R Main driver script: loops over counties × years, calls process\_county\_year(), and writes out cluster-level yield tables + diagnostic figures.

\## Demo data (no-Earth-Engine example)

The repository ships a small demo dataset under inst/extdata/ to allow testing the workflow without running Google Earth Engine:

inst/extdata/nass\_demo.csv Small NASS-style table with county, year, crop, harvested area (ha) and yield (kg/ha) for a few Texas counties.

inst/extdata/tx\_counties\_demo.gpkg A simplified Texas High Plains counties layer in EPSG:4326, containing the counties used in the demo (e.g. Deaf Smith, Gaines).

inst/extdata/WH\_weight\_DeafSmith\_16.tif Pre-computed normalised weight raster for winter wheat in Deaf Smith County, TX, for a representative season. Each pixel value represents the relative contribution to county-level wheat production (sum = 1 over eligible pixels).

inst/extdata/CO\_weight\_Gaines\_21.tif Pre-computed weight raster for cotton in Gaines County, TX, for a representative season.

inst/extdata/CDL\_DeafSmith\_16.tif A CDL excerpt for Deaf Smith County (crop codes and attributes as in the original CDL). Used in the demo to illustrate the crop mask and weight surface.

inst/extdata/WH\_dis-yield.pdf Example figure: disaggregated winter-wheat yield in Deaf Smith County with histogram and NASS comparison.

inst/extdata/CO\_dis-yield.pdf Example figure: disaggregated cotton yield in Gaines County with histogram and NASS comparison.

These files are only for illustration and unit tests; they are not intended as an authoritative dataset.

\## Quick start

0\. Install dependencies

install.packages(c("sf","terra","dplyr","rgee", ...))

1\. Initialize Earth Engine

&nbsp;Optional: manually point reticulate to the Python env used by rgee

Sys.setenv(

&nbsp; RETICULATE\_PYTHON = "~/AppData/Local/r-miniconda/envs/rgee\_py311/python.exe"

)

library(reticulate)

py\_config()

library(rgee)

This installs a miniconda env and earthengine-api into it

rgee::ee\_install()

One-time authentication (will open browser)

rgee::ee\_Authenticate()

Initialize in your R session (usually once per session)

rgee::ee\_Initialize(drive = TRUE)

2\. Prepare NASS tables, county shapefile, and CDL mosaics

(see data/README\_data.md)

3\. Compute weights (optional, precompute)

source("R/02\_build\_spatial\_weights.R")

source("R/07\_compute\_weights\_batch.R")

compute\_weights\_batch(...)

4\. Run disaggregation

source("R/03\_segment\_fields\_pseudo\_units.R")

source("R/04\_disagg\_units\_yield\_safe.R")

source("R/05\_cluster\_units\_prod\_zones.R")

source("R/06\_process\_county\_year.R")

source("R/08\_run\_disagg\_clustering.R")

#Working EX for cotton in Gaines County in 2021 in TX 
#Step 1 building weight layer 

library(NAYD)
library(sf)
library(dplyr)
library(purrr)
library(units)
library(stats)  # hclust, cutree
library(terra)
library(sf)
library(mapview)
library(dplyr)
library(data.table)
library(stringr)
library(terra)
library(scales)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(patchwork)
library(stringr)
library(rlang)
library(viridisLite) # palettes
library(cowplot)
library(viridis)
library(readxl)

# 
#- the correct order here 
Sys.setenv(
  RETICULATE_PYTHON = "C:/Users/ahmed.attia/AppData/Local/r-miniconda/envs/rgee311/python.exe"
) #update to your local path
library(reticulate)
py_config()
library(rgee)
ee_Authenticate()
ee_Initialize()

nass_csv <- "C:/Users/ahmed.attia/OneDrive - Texas A&M AgriLife/ML_study/SpatialData/CSV/TX_CO_PLAN_HARV_Y.csv"  # <-- update
nass_df <- read.csv(nass_csv)

season_start <- "06-01"
season_end   <- "09-30"
crop_code <- c(2)
crop_name <- "Cotton"
CRS_TARGET  <- "EPSG:3857"

THP_shp <- st_read(system.file("extdata", "tx_counties_demo.gpkg", package = "NAYD"), quiet = TRUE) |> st_make_valid() |> st_transform(4326)
county_sf <- THP_shp[THP_shp$COUNTY == "Gaines",]

#the CDL_2021_48.tif is clipped CDL layer for Gaines county in 2021. 

get_weights_for_county_year(county_sf = county_sf,
                            county_name = "Gaines", sensor="LANDSAT", year=2021, 
                            season_start, 
                            season_end, CRS_TARGET, crop_code, crop_name, 
                            a_ndvi=0.7, a_et=0.3, 
                            NASS_planted_area_ha = nass_df$area_planted_ha[(nass_df$year %in% 2021 & nass_df$county == norm_county(county_name))],
                            NASS_harvest_area_ha = nass_df$area_harvested_ha[(nass_df$year %in% 2021 & nass_df$county == norm_county(county_name))], 
                            lower_thresh = 0.6, upper_thresh = 2, area_buffer_frac = 0.05, 
                            cdl_dir="C:/NAYD/inst/extdata", cache_dir="Path_to_outPut_dir")



