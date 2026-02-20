# ------------------------------------------------------------------------------
# find areas to distrub
# 
# adapted from 005_LAI_seasonal_curves_figure
# Dec 1, 2025
# adpated Feb 19 for new 
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# library 
# ------------------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)
library(units)
library(docstring)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# set wd 
setwd(file.path("C:","Users", "blbouf", "Sync", "Paper2"))
output_path <- file.path(".", "data", "disturbance_recovery_scenarios")

allHRUs_path <- file.path(".", "data", "LCC_HRU_files", "Feb19", "HRUs_1923_2023") # was nov 16
hru_2023_path <- file.path(allHRUs_path, "HRU_2023.shp")
hru_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp")

# path to other continuous data (DEM, slope, aspect)
dem_path <- file.path(".", "data", "src", "lidar_derived", "dem.tif")
# slp_path <- file.path(".", "data", "src", "lidar_derived", "slope.tif")
# asp_path <- file.path(".", "data", "src", "lidar_derived", "aspect.tif")

# ------------------------------------------------------------------------------
# data 
# ------------------------------------------------------------------------------

dem <- rast(dem_path)
hrus <- st_read(hru_path)


# disturb files 
# --------------------------------------------------------------30
# high 30
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 30,
                 elevation_q50 = "high",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# low 30
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 30,
                 elevation_q50 = "low",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# --------------------------------------------------------------15
# high 15
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 15,
                 elevation_q50 = "high",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# low 15
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 15,
                 elevation_q50 = "low",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# --------------------------------------------------------------20
# high 20
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 20,
                 elevation_q50 = "high",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# low 20
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 20,
                 elevation_q50 = "low",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# --------------------------------------------------------------10
# high 10
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 10,
                 elevation_q50 = "high",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# low 10
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 10,
                 elevation_q50 = "low",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# --------------------------------------------------------------40
# high 40
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 40,
                 elevation_q50 = "high",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)
# low 40
disturb_polygons(dem = dem,
                 hrus = hrus,
                 area_percent = 40,
                 elevation_q50 = "low",
                 output_path = output_path,
                 write_dem_binary = FALSE, 
                 write_output = TRUE)

# ------------------------------------------------------------------------------
# create recovery files 
# ------------------------------------------------------------------------------

disturbance_files <- list.files(output_path, 
                                pattern = "0yrs",
                                full.names = TRUE) 

disturbance_shp_files <- disturbance_files[grepl(pattern = "0yrs.shp$", disturbance_files)]

recovery_years <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50) # originally was 5, 10, 20, ... every 10 yrs.. 

lapply(disturbance_shp_files, 
       recover_disturbed_polygons,
       recovery_years
)

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

disturb_polygons <- function(dem, 
                             hrus, 
                             area_percent = 30,
                             elevation_q50, 
                             output_path, 
                             write_dem_binary = FALSE,
                             write_output = TRUE){
  
  #' @description
  #' Disturb X% of forested area in watershed either in the high elevation 
  #' or low elevation areas 
  #' 
  #' @param dem dem raster of study area 
  #' @param hrus sf object of hrus
  #' @param area_percent percent area of watershed to disturb (out of 100)
  #' @param elevation_q50 string "high" or "low" representing region to disturb based on 50th percentile area
  #' @param output_path output path 
  #' @param write_dem_binary TRUE/FALSE
  #' @param write_output TRUE/FALSE
  #' @return returns hru sf object with disturbed polygons set with an age of 0
  # ----------------------------------------------------------------------------
  # compute 50th percentile
  # ----------------------------------------------------------------------------
  
  med_elev <- quantile(values(dem), probs = 0.5, na.rm = TRUE)
  dem_binary <- dem
  
  if (elevation_q50 == "high"){ # high elevation 
    
    dem_binary[dem > med_elev] <- 1
    dem_binary[dem <= med_elev] <- 0
    
  } else if (elevation_q50 == "low"){ # low elevation 
    
    dem_binary[dem < med_elev] <- 1
    dem_binary[dem >= med_elev] <- 0
    
  } else {
    return("please designate elevation_q50 as 'high' or 'low'")
  }
  
  # write out binary layer if applicable
  if (write_dem_binary == TRUE){
    writeRaster(dem_binary,
                file.path(output_path, "binary_DEM_50pct.tif"), 
                overwrite=TRUE)
  }
  
  # ----------------------------------------------------------------------------
  # sample the binary DEM inside each polygon 
  # ----------------------------------------------------------------------------
  
  # extract = returns summary for each polygon by default
  vals <- terra::extract(dem_binary, hrus, fun = median)
  
  # add mean high-elevation fraction per polygon
  hrus$p50 <- vals[,2] 
  
  hru_p50_missing <- hrus$ELEVATI[is.na(hrus$p50)]
  
  # deal with hrus that have missing elev class 
  if (elevation_q50 == "high"){
    hru_p50_missing[hru_p50_missing > med_elev] <- 1 # high elevation
    hru_p50_missing[hru_p50_missing <= med_elev] <- 0 # low elevation
  } else if (elevation_q50 == "low"){
    hru_p50_missing[hru_p50_missing < med_elev] <- 1 # low elevation
    hru_p50_missing[hru_p50_missing >= med_elev] <- 0 # high elevation 
  }
  
  hrus$p50[is.na(hrus$p50)] <- hru_p50_missing 
  
  # ----------------------------------------------------------------------------
  # weighted random selection of polygons 
  # ----------------------------------------------------------------------------
  
  # Define weights
  hrus$w <- hrus$p50
  
  hrus <- hrus[hrus$VEG_CLA %in% c("FOREST_MS", "FOREST_IDF", "FOREST_ESSF"), ]
  
  # ----------------------------------------------------------------------------
  # select polygons
  # ----------------------------------------------------------------------------
  
  total_area <- sum(st_area(hrus)) %>% drop_units()
  target_area <- area_percent/100 * total_area
  
  selected <- c()
  selected_area <- 0
  
  set.seed(1117)  # reproducible
  
  while(selected_area < target_area) {
    
    i <- sample(1:nrow(hrus), size = 1, prob = hrus$w)
    
    if(!(i %in% selected)) {
      selected <- c(selected, i)
      selected_area <- sum(st_area(hrus[selected, ])) %>% drop_units()
    }
  }
  
  selected_polys <- hrus[selected, ] 
  selected_polys$age <- 0

  if (write_output == TRUE){
    st_write(selected_polys,
             file.path(output_path, paste0("p50_", elevation_q50, "_", area_percent, "_", "0yrs.shp")),
             append = FALSE)
  }

  return(paste0(elevation_q50, "\npercent area: ", area_percent))
  
}

recover_disturbed_polygons <- function(path_to_disturbed_polygons, 
                                       recovery_years){
  #' @description
    #' function to take hru layers with new disturbance and to set age based on additional recovery years
    #' 
    #' @param path_to_disturbed_polygons full file path to disturbed HRU shapefiles
    #' @param recovery_years list of years to simulate recovery
    #' 
  
  hrus <- st_read(path_to_disturbed_polygons)
  
  for (a in 1:length(recovery_years)){
    
    age <- recovery_years[a]
    
    recovering_hrus <- hrus
    recovering_hrus$age[hrus$age == 0] <- age
    
    new_filename <- stringr::str_replace(path_to_disturbed_polygons, 
                                         "0yrs", 
                                         paste0(age, "yrs"))
    st_write(recovering_hrus,
             file.path(new_filename), 
             append = FALSE)
  }
  
  return("done!")
}
                                       



