# usda-nass-yield-disaggregating

This repository provides an R workflow to disaggregate official, coarse-scale yield statistics (e.g. USDA NASS county yields) to field-scale “harvest units” using a combination of:
land-cover / crop-type rasters (e.g. USDA Cropland Data Layer),
remote-sensing metrics (NDVI, ET) from Google Earth Engine, and
reported area and yield for each administrative unit.
The code can be applied to any crop and region where equivalent inputs are available.
In the documentation and examples we demonstrate the workflow using winter wheat and cotton in Texas, but the methods are not limited to these crops or locations.

## Features

- Builds county- and crop-specific masks from USDA Cropland Data Layer (CDL)
  with consistency checks against NASS planted/harvested areas.
- Derives seasonal NDVI and ET composites from Landsat / OpenET.
- Computes pixel-level weights combining NDVI and ET.
- Segments fields and large pseudo-fields into harvest units with area constraints.
- Disaggregates NASS county yield across units using a beta-weighting scheme
  while preserving county means.
- Clusters harvest units into ~5–10 km production zones.
- Validates segmented yields against independent variety-testing trials.

## Repository structure

- `R/01_ndvi_et_season.R` – NDVI/ET seasonal composites (rgee/OpenET).
- `R/02_get_weights_county_year.R` – CDL mask + NDVI/ET weight generator.
- `R/03_segment_fields_pseudo_units.R` – field segmentation and tiling.
- `R/04_disagg_units_yield_safe.R` – NASS-consistent yield assignment to units.
- `R/05_cluster_units_5km.R` – cluster units into production zones.
- `R/06_process_county_year.R` – wrapper for one county-year.
- `R/07_compute_weights_batch.R` – batch weights for many counties/years.
- `R/08_run_disagg_clustering.R` – main loop + segmentation + export.

## Requirements

- R >= 4.3
- Packages: `sf`, `terra`, `dplyr`, `data.table`, `ggplot2`, `rgee`,
  `reticulate`, `readxl`, `stringr`, `purrr`, `viridis`, `patchwork`, etc.
- Google Earth Engine account and `rgee` configured.
- Access to:
  - Texas county boundary shapefile (EPSG:4326)
  - CDL mosaics (30 m) for 2008–present, in EPSG:5070
  - NASS county-level planted/harvested area and yield tables

## Quick start

```r
# 0. Install dependencies
# install.packages(c("sf","terra","dplyr","rgee", ...))

# 1. Initialize Earth Engine
rgee::ee_Initialize(drive = TRUE)

# 2. Prepare NASS tables, county shapefile, and CDL mosaics
# (see data/README_data.md)

# 3. Compute weights (optional, precompute)
source("R/02_get_weights_county_year.R")
source("R/07_compute_weights_batch.R")
compute_weights_batch(...)

# 4. Run disaggregation
source("R/03_segment_fields_pseudo_units.R")
source("R/04_disagg_units_yield_safe.R")
source("R/05_cluster_units_5km.R")
source("R/06_process_county_year.R")
source("R/08_run_disagg_clustering.R")
