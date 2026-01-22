
source("R_01_ndvi_et_season.R")
source("R_02_build_spatial_weights.R")

# --- Cotton
nass_csv <- "NASS_Cotton.csv"  # <-- update
nass_df <- read.csv(nass_csv)
season_start <- "06-01"
season_end   <- "09-30"
crop_code <- c(2)
crop_name <- "Cotton"
CRS_TARGET  <- "EPSG:3857"
# --- Winter wheat
nass_csv <- "NASS_wheat.csv"  # <-- update
nass_df <- read.csv(nass_csv)
season_start <- "03-01"  
season_end   <- "05-30"
crop_code   <- c(24)  
crop_name <- "Winter_wheat"
CRS_TARGET  <- "EPSG:3857"

years <- c(2018:2024)  #desired years 

for (yr in years) {
  
  cdl_fp_year <- cdl_for_processing(yr)
  if (is.na(cdl_fp_year)) { 
    cat("Year", yr, ": no CDL mosaic\n"); 
    next 
  }
  cdl_TX <- terra::rast(cdl_fp_year)
  
  sub_yr <- var_test[var_test$year == yr,]
  cnties <- unique(sub_yr$county)
  
  for (c in cnties) {
    cname     <- norm_county(c)
    county_sf <- test_counties_sf2[test_counties_sf2$COUNTY == c, ]
    
    out_fp <- .cache_path("weight", cname, yr,norm_county(crop_name),sensor)
    if (file.exists(out_fp)){
      cat(sprintf(" ✓ skip %s_%s_%d_%s_%s already exists\n",
                  "weight",cname,yr,norm_county(crop_name),sensor))
      next
    }
    
    # --- lookup NASS harvested area for this county-year ---
    row_i <- nass_df[nass_df$county == cname & nass_df$year == yr, ]
    
    if (nrow(row_i) == 0L) {
      # no data for this county-year → handle skip
      cat(sprintf( " -SKIP %s %d : no NASS row\n", cname, yr))
      next 
    } else {
      # keep only the row with largest harvested area
      row_i <- row_i[which.max(row_i$area_harvested_ha), , drop = FALSE]
    }
    
    NASS_planed_area_ha <- row_i$area_planted_ha
    NASS_harvest_area_ha <- row_i$area_harvested_ha
    
    # --- skip if no valid harvested area ---
    if (!is.finite(NASS_harvest_area_ha) || is.na(NASS_harvest_area_ha) || NASS_harvest_area_ha <= 0) {
      cat(sprintf("  - SKIP %s %d: no NASS harvested area\n", cname, yr))
      next
    }
    
    # --- compute weights ---
    wt <- tryCatch(
      get_weights_for_county_year(
        county_sf       = county_sf,
        county_name     = cname,
        sensor          = sensor,
        year            = yr,
        season_start    = season_start,
        season_end      = season_end,
        TX_shp        = TX_shp,
        CRS_TARGET      = CRS_TARGET,
        crop_code       = crop_code,
        crop            = crop_name,
        a_ndvi          = 0.70,
        a_et            = 0.30,
        NASS_planed_area_ha = NASS_planed_area_ha,
        NASS_harvest_area_ha  = NASS_harvest_area_ha,
        lower_thresh = 0.5,
        upper_thresh = 3.50,      # max ratio of NASS planted area 
        area_buffer_frac = 0.2,  
        pre2008_mask_strategy = "auto"
      ),
      error = function(e) {
        cat(sprintf("  - ERROR %s %d: %s\n", cname, yr, conditionMessage(e)))
        NULL
      }
    )
    
    if (is.null(wt) || !("weight" %in% names(wt))) {
      cat(sprintf("  - SKIP %s %d: weights NULL\n", cname, yr))
      next
    }
    
    wsum <- as.numeric(terra::global(wt$weight, "sum", na.rm = TRUE))
    cat(sprintf("  ✓ %s %d (weights cached, sum = %.6f, NASS_harvest_area_ha = %.1f)\n",
                cname, yr, wsum, NASS_harvest_area_ha))
  }
}

