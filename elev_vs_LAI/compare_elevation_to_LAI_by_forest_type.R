# ------------------------------------------------------------------------------
# script to clip fCOVER to study area and join into a single raster stack 
# 
# date: november 25, 2025
# author: Brianne Boufford 
# modified from scale_lai_stack_rasters.R
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# load 
library(dplyr)
library(terra)
library(sf)
library(stringr)
library(ggplot2)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
data_path <- file.path(".", "data")
lai_path <- file.path(data_path, "src", "LAI", "adjusted_hls_lai.tif")
bec_path <- file.path(data_path, "src", "BEC", "bec_zones_tc.shp")
dem_path <- file.path(data_path, "src", "lidar_derived", "dem.tif")
catchment_path <- file.path(data_path, "src", "TC_catchment", "catchments.shp")
p1_hru_path <- file.path(data_path, "src", "HRU_paper1", "HRU_updated_oct7", "HRUs_1923_2023")

# ------------------------------------------------------------------------------
# data 
# ------------------------------------------------------------------------------
lai <- rast(lai_path)
bec <- st_read(bec_path)
dem <- rast(dem_path)
c <- st_read(catchment_path)
p1_files <- list.files(p1_hru_path)

# rasterize bec zones
bec_rast <- rasterize(bec, lai, field = "ZONE")

# transform catchment 
c_proj <- st_transform(c, crs = crs(lai))

# resample and mask dem 
dem_resamp <- resample(dem, lai)
dem_masked <- mask(dem_resamp, c_proj)

# only get LAI files that are from 2014-2022
hru_files <- grep("20(1[4-9]|2[0-2]).shp$", p1_files, value = TRUE)

# subset LAI files to years 2014-2022 
pattern <- "20(1[4-9]|2[0-2])_07$"
peak_lai_idx <- grepl(pattern, names(lai))
peak_lai <- lai[[peak_lai_idx]]

yrs <- 2014:2022 

results <- lapply(yrs, 
                  summarize_lai_by_elevation, 
                  dem = dem_resamp,
                  peak_lai= peak_lai,
                  hru_files = hru_files,
                  hru_files_base_path = p1_hru_path) %>%
  do.call(rbind, .)

results_clean <- na.omit(results) %>%
  unique()

results_grouped <- results_clean %>%
  group_by(elevation) %>%
  summarize(lai_av = mean(LAI),
            forest_type = first(forest_type))

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------

ggplot(results_grouped, aes(x = lai_av, y = elevation, color = forest_type)) +
  geom_point(alpha = 0.6) +                 # scatter points
  #geom_smooth(method = "lm", se = TRUE) +   # fitted regression line
  labs(x = "LAI",
    y = "Elevation (m)",
    title = "Elevation vs LAI by Forest Type"
  ) +
  theme_bw()

# ------------------------------------------------------------------------------
# function to get elevation and LAI for each mature forest HRU  
# ------------------------------------------------------------------------------

summarize_lai_by_elevation <- function(yr, 
                                       dem,
                                       peak_lai, 
                                       hru_files, 
                                       hru_files_base_path){
  
  peak_lai_i <- grepl(yr, names(peak_lai))
  p_lai <- peak_lai[[peak_lai_i]]
  
  hru_i <- hru_files[grepl(yr, hru_files)]
  hru_path <- file.path(hru_files_base_path, hru_i)
  hru <- st_read(hru_path)
  
  hru_forested <- hru[hru$VEG_CLA == "mature", ]
  hru_lai <- extract(p_lai, hru_forested, fun = mean)
  hru_elev <- extract(dem_masked, hru_forested, fun = mean)
  
  results <- cbind(hru_elev, hru_lai, hru_forested$frst_cl) %>% 
    select(-c("ID")) %>% 
    select(-c("ID"))
  
  names(results) <- c("elevation", "LAI", "forest_type")
  
  return(results)
  
}
