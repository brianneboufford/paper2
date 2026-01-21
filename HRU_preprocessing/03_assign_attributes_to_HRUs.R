# ------------------------------------------------------------------------------
# extract attributes from new HRU polygons
# 
# author: Brianne Boufford 
# date: November 10th, 2025
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# libraries 
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(units)

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
# working dir 
project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

# path to new hrus 
hru_path <- file.path(".", "data", "HRU_delineation", "interm_unions", "6_all_data_filled_1ha.shp")
hru_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "HRU_no_slope_aspect_Q.shp")

# read in old HRUs for attribute column names
og_hru_path <- file.path(".", "data", "src", "HRU_paper1", "HRU_Poly_attributes.shp")

# path to other continuous data (DEM, slope, aspect)
dem_path <- file.path(".", "data", "src", "lidar_derived", "dem.tif")
slp_path <- file.path(".", "data", "src", "lidar_derived", "slope.tif")
asp_path <- file.path(".", "data", "src", "lidar_derived", "aspect.tif")

# output path
outpath <- file.path(".", "data", "HRU_delineation")

# ------------------------------------------------------------------------------
# read data
# ------------------------------------------------------------------------------
hrus <- st_read(hru_path)
og_hrus <- st_read(og_hru_path)
dem <- rast(dem_path)
slp <- rast(slp_path)
asp <- rast(asp_path)

atts <- names(og_hrus)
# --> "ID"       "AREA" X    "ELEVATI" X "LATITUD"  "LONGITU"  "BASIN_I" X "LAND_US"  X
# "VEG_CLA" X "SOIL_PR" X "AQUIFER" X "TERRAIN"X  "SLOPE"X    "ASPECT"  X "label" X   "geometry"X

hrus_new <- hrus %>% 
  select("LAND_US", "VEG_CLA", "SOIL_PR", "AQUIFER", "TERRAIN", "BASIN_I",
         "ZONE", "fire_yr", "harv_yr")
hrus_new$AREA <- st_area(hrus_new)
hrus_new$ELEVATI <- terra::extract(dem, vect(hrus_new), fun = mean, ID = FALSE, na.rm = TRUE)$dem
hrus_new$SLOPE <- terra::extract(slp, vect(hrus_new), fun = mean, ID = FALSE, na.rm = TRUE)$slope

# convert aspect to radians to get HRU-based average to avoid circularity of North 
# facing aspects (average of 0 and 360 is 180 (South) even though its true north)
circ_mean <- function(a){
  a <- a*pi/180
  mean_x <- mean(cos(a), na.rm = TRUE)
  mean_y <- mean(sin(a), na.rm = TRUE)
  ang <- atan2(mean_y, mean_x) * 180/pi
  if(is.na(ang)){
    ang <- NA
  } else if (ang < 0){ 
    ang <- ang + 360
  }
  return(ang)
}

circ_asp <- terra::extract(asp, hrus_new, fun=circ_mean)
hrus_new$ASPECT <- circ_asp$aspect

# hrus with missing topography data 
hrus_missing <- hrus_new[is.na(hrus_new$ELEVATI), ]
hrus_missing <- st_buffer(hrus_missing, dist = 500)
hrus_missing$ELEVATI <- terra::extract(dem, vect(hrus_missing), fun = mean, ID = FALSE, na.rm = TRUE)$dem
hrus_missing$SLOPE <- terra::extract(slp, vect(hrus_missing), fun = mean, ID = FALSE, na.rm = TRUE)$slope
hrus_missing$ASPECT <- terra::extract(asp, vect(hrus_missing), fun = mean, ID = FALSE, na.rm = TRUE)$aspect

# replace in hru dataset 
hrus_new[is.na(hrus_new$ELEVATI), "SLOPE"] <- hrus_missing$SLOPE
hrus_new[is.na(hrus_new$ELEVATI), "ASPECT"] <- hrus_missing$ASPECT
hrus_new[is.na(hrus_new$ELEVATI), "ELEVATI"] <- hrus_missing$ELEVATI

# clean other data from OG dataset 
hrus_new$LAND_US[is.na(hrus_new$LAND_US)] <- "FOREST"
hrus_new$VEG_CLA[is.na(hrus_new$VEG_CLA)] <- "FOREST"
hrus_new$BASIN_I[is.na(hrus_new$BASIN_I)] <- 1
hrus_new$SOIL_PR[is.na(hrus_new$SOIL_PR)] <- "SOIL"

hrus_new$AREA <- drop_units(hrus_new$AREA)/(1000*1000) # in km2
hrus_new$ID <- 1:length(hrus_new$LAND_US)
hrus_new$LONGITU <- st_coordinates(st_transform(st_centroid(hrus_new$geometry), crs=4326))[, "X"]
hrus_new$LATITUD <- st_coordinates(st_transform(st_centroid(hrus_new$geometry), crs=4326))[, "Y"]

hrus_new$AREA <- round(hrus_new$AREA, digits = 2)
hrus_new$ELEVATI <- round(hrus_new$ELEVATI, digits = 0)
hrus_new$SLOPE <- round(hrus_new$SLOPE, digits = 0)
hrus_new$ASPECT <- round(hrus_new$ASPECT, digits = 0)

# make forest distinct by BEC Zone
hrus_new$LAND_US[hrus_new$LAND_US == "FOREST" & hrus_new$ZONE == "MS"] <- "FOREST_MS"
hrus_new$LAND_US[hrus_new$LAND_US == "FOREST" & hrus_new$ZONE == "IDF"] <- "FOREST_IDF"
hrus_new$LAND_US[hrus_new$LAND_US == "FOREST" & hrus_new$ZONE == "ESSF"] <- "FOREST_ESSF"

# apply to VEG CLASS too
hrus_new$VEG_CLA <- hrus_new$LAND_US

# reorder columns 
hrus_new <- hrus_new[, c("ID", "AREA", "ELEVATI", "LATITUD", "LONGITU", "BASIN_I", 
                         "LAND_US", "VEG_CLA", "SOIL_PR", "AQUIFER", "TERRAIN", 
                         "SLOPE", "ASPECT", "ZONE", "fire_yr", "harv_yr", "geometry")]

# write to new shapefile
st_write(hrus_new, 
         file.path(outpath, "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp"),
         append = FALSE)
