# -----------------------------------------------------------------------------_
# Join CHM and canopy cover (pzabove2m) to 2015 catchment data 
# 2014 age data has quite a few inconsistencies in the ESSF forested area 
# 
# Date created: Feburary 9, 2025
#
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(smoothr)
library(stringr)
library(ggplot2)
library(future)
library(future.apply)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

output_prod_path <- file.path(".", "data", "forest_params_by_age")
write_interm_data <- TRUE

# paths
if(TRUE){
  
  #fc_path <- file.path(".", "data", "src", "ntems", "percentage_first_returns_above_2m")
  h_path <- file.path(".", "data", "src", "lidar_derived", "chm.tif")
  
  # forest cover
  fc_path <- file.path(".", "data", "src", "lidar_derived", "standard_metrics.tif")
  
  # dem 
  dem_path <- file.path(".", "data", "src", "lidar_derived", "dem.tif")
  
  # bec_path 
  bec_path <- file.path(".", "data", "src", "BEC", "bec_zones_tc.shp")
  
  # hru path
  hru_path <- file.path(".", "data", "LCC_HRU_files", "Nov10", "HRUs_1923_2023", "HRU_2023.shp")
  
  # last change year polygon later 
  lcy_poly_path <- file.path(".", "data", "src", "ntems", "disturbance_year_polygon.shp")
  
}

if(TRUE){
  
  # resample everything to 20m --> same resolution as the std metrics 
  
  # forest cover
  std_metrics <- rast(fc_path)
  fc <- std_metrics$pzabove2
  
  # height
  h_files <- rast(h_path) %>% 
    terra::resample(fc, method = "average")
  
  # dem
  dem <- rast(dem_path) %>% 
    terra::resample(fc, method = "average")
  
  age_combined <- rast(file.path(".", "data", "age", "age_combined.tif")) %>%
    terra::resample(fc, method='mode')
  
  age_combined <- age_combined$year_2015
  
  zone_raster <- rast(file.path(".", "data", "src", "BEC", "zone_raster.tif")) %>%
    terra::resample(fc, method='mode')
  
  # skipping for now
  hru_2023 <- st_read(hru_path) %>%
    st_transform(., crs = crs(fc))
  
  hru_rast <- rasterize(hru_2023, age_combined, field = "frst_cl")
  
  lcy_poly <- st_read(lcy_poly_path)
  lcy_poly$fid[lcy_poly$Last_CY == 0] <- 3333
  lcy_rast <- rasterize(lcy_poly, fc, field = "fid")
  
  age_i <- age_combined
  
  lai_grp <- dem
  lai_grp[dem >= 1560] <- 1
  lai_grp[dem < 1560] <- 0
  
  h_i <- h_files %>% 
    terra::project(., fc, method='average')
  
  fc_i <- fc %>% 
    terra::project(., fc, method='average')
  
  data_i <- c(age_i, zone_raster, hru_rast, dem, lcy_rast, lai_grp) 
  data_ii <- c(h_i, fc_i)
  
  if(write_interm_data){
    
    data_all <- c(data_i, data_ii)
    writeRaster(data_all, 
                file.path(".", "data", "interm", "ALS_fparams_joined", "ALS_fparams_stack_20m.tif"))
  }

  df_1 <- as.data.frame(data_i, xy=TRUE) %>%
    na.omit() 
  df_2 <- as.data.frame(data_ii, xy=TRUE)
  
  df_i <- merge(df_1, df_2, by=c("x", "y"))
  
  names(df_i) <- c("x", "y", "age", "id", "frst_cl", "dem", "lcy_id", "lai_grp", "height", "fc")
  
  df_i <- df_i[df_i$frst_cl %in% c("ESSF", "MS", "IDF"), ]
  
  # write fully merged data 
  write.csv(df_i,
            file.path(output_prod_path, "fc_height_2015_feb11.csv"),
            row.names = FALSE)
}

# ------------------------------------------------------------------------------
# sample n = 1000 from each age 
# ------------------------------------------------------------------------------

set.seed(1113)

df_i <- read.csv(file.path(output_prod_path, "fc_height_2015_feb11.csv"))


# sampled_data <- df_i %>% 
#   group_by(age, lai_grp) %>%
#   slice_sample(n = 1000) %>%
#   ungroup()
# 
# sampled_forested_rse <- df_i %>%
#   group_by(age, lai_grp) %>%
#   slice_sample(n = 1000) %>%
#   ungroup() %>%
#   filter(age < 100) %>%
#   group_by(age, lai_grp) %>% 
#   summarize(rse_fcover = AnglerCreelSurveySimulation::calculate_rse(na.omit(fc)),
#             n_fc = length(fc), 
#             rse_height = AnglerCreelSurveySimulation::calculate_rse(na.omit(height)),
#             n_h = length(height), 
#   )

# summarize the number of samples per age in each group (n)
# sampled_frst_n <- sampled_data %>%
#   group_by(age, lai_grp) %>%
#   summarize(n = length(age))

########################################################
# write out 1000 values for each yaer 
# write.csv(sampled_data,
#           file.path(output_prod_path, "sampled_fc_height_2015_feb10.csv"),
#           row.names = FALSE)
