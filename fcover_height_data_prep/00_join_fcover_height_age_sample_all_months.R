# -----------------------------------------------------------------------------_
# script to sample NTEMS FCOVER and veg height by age - only forested pixels
# keep all winter LAI as monthly for now, calculate median later in the process  ****
#
#
# adapated from 03_prep_sampled_fcover_age_data.R and then adapted from 
# 00_join_fcvoer_age_and_sample.R
# get fCOVER by age and pull out alpine, wetland, and shrub classes 
# 
# february 13, 2025
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

# paths
if(TRUE){
  
  # ntems fcover and height
  fc_ntems_path <- file.path(".", "data", "src", "ntems", "percentage_first_returns_above_2m")
  h_ntems_path <- file.path(".", "data", "src", "ntems", "elev_p95")
  
  # adjusted fcover and LAI from leaf data 
  fc_leaf_path <- file.path(".", "data", "fCOVER_clean", "fCOVER_stacked", "fcover_stacked_TC_area.tif")
  lai_leaf_path <- file.path(".", "data", "LAI", "adjusted_hls_lai.tif")
  
  # chm height
  h_chm_path <- file.path(".", "data", "src", "lidar_derived", "chm.tif")
  
  # dem 
  dem_path <- file.path(".", "data", "src", "lidar_derived", "dem.tif")
  
  # bec_path 
  bec_path <- file.path(".", "data", "src", "BEC", "bec_zones_tc.shp")
  
  # hru path
  hru_path <- file.path(".", "data", "LCC_HRU_files", "Nov10", "HRUs_1923_2023", "HRU_2023.shp")
  
  # last change year polygon later 
  lcy_poly_path <- file.path(".", "data", "src", "ntems", "disturbance_year_polygon.shp")
  
  # 2014 chm 
  chm_path <- file.path(".", "data", "src", "lidar_derived", "chm.tif")
  
  # 2014 standard metrics 
  std_metrics_path <- file.path(".", "data", "src", "lidar_derived", "standard_metrics.tif")
  
  
}

if(TRUE){
  
  # ntems fcover
  fc_ntems_files <- list.files(fc_ntems_path, pattern = ".dat$", full.names = TRUE) # list raster files 
  fc_ntems_files <- str_subset(fc_ntems_files, "201[5-9]|202[0-1]") # subset to 2014:2021
  
  fc_eg <- rast(fc_ntems_files[1]) # get an example 
  
  # ntems height
  h_ntems_files <- list.files(h_ntems_path, pattern = ".dat$", full.names = TRUE) # list raster files  
  h_ntems_files <- str_subset(h_ntems_files, "201[5-9]|202[0-1]") # subset to 2014:2021
  
  # leaf LAI 
  leaf <- rast(lai_leaf_path)
  leaf[leaf == 0] <- NA
  leaf <- leaf[[-c(1:12)]]
  
  # leaf fcover 
  fcov <- rast(fc_leaf_path)
  
  # subset to july aka peak fcover
  peak_fcov <- fcov[[grepl("_07", names(fcov))]]
  peak_fcov <- peak_fcov[[1:(nlyr(peak_fcov) - 2)]]
  
  # unique years
  years <- 2015:2021
  
  age_combined <- rast(file.path(".", "data", "age", "age_combined.tif"))
  age_combined <- age_combined[[-1]]
  
  # dem
  dem <- rast(dem_path) %>% 
    terra::resample(age_combined, method = "average")
  
  zone_raster <- rast(file.path(".", "data", "src", "BEC", "zone_raster.tif")) %>%
    terra::resample(age_combined, method='mode')
  
  # skipping for now
  hru_2023 <- st_read(hru_path) %>%
    st_transform(., crs = crs(age_combined))
  
  chm <- rast(chm_path) %>% 
    terra::resample(age_combined, method='max') # taking max height because I want to remove edge pixels where LAI is too coarse to see forest next to 
  # disturbed or regenerating stand, so If I take the max then I can filter out these pixels, If I use median or average it will just smooth them over 
  # which I don't want bc age will be inaccurate for those conditions 
  
  fc <- rast(std_metrics_path) %>% 
    subset("pzabove2") %>% 
    terra::resample(age_combined, method='max')
  
  hru_rast <- rasterize(hru_2023, age_combined, field = "frst_cl")
  
  lcy_poly <- st_read(lcy_poly_path)
  lcy_poly$fid[lcy_poly$Last_CY == 0] <- 3333
  lcy_rast <- rasterize(lcy_poly, age_combined, field = "fid")
  
  # subset monthly LAI data ------------------------------------------
    
  # may 
  may_lai <- leaf[[grepl("_05", names(leaf))]]
  may_lai <- may_lai[[1:(nlyr(may_lai) - 2)]]
  
  # june 
  june_lai <- leaf[[grepl("_06", names(leaf))]]
  june_lai <- june_lai[[1:(nlyr(june_lai) - 2)]]
  
  # july 
  july_lai <- leaf[[grepl("_07", names(leaf))]]
  july_lai <- july_lai[[1:(nlyr(july_lai) - 2)]]
  
  # august 
  aug_lai <- leaf[[grepl("_08", names(leaf))]]
  aug_lai <- aug_lai[[1:(nlyr(aug_lai) - 2)]]
  
  # sept
  sept_lai <- leaf[[grepl("_09", names(leaf))]]
  sept_lai <- sept_lai[[1:(nlyr(sept_lai) - 2)]]
  
  # oct
  oct_lai <- leaf[[grepl("_10", names(leaf))]]
  oct_lai <- oct_lai[[1:(nlyr(oct_lai) - 2)]]
  
  # april
  apr_lai <- leaf[[grepl("_04", names(leaf))]]
  apr_lai <- apr_lai[[1:(nlyr(apr_lai) - 2)]]

  # nov 
  nov_lai <- leaf[[grepl("_11", names(leaf))]]
  nov_lai <- nov_lai[[1:(nlyr(nov_lai) - 2)]]
  
  # dec
  dec_lai <- leaf[[grepl("_12", names(leaf))]]
  dec_lai <- dec_lai[[1:(nlyr(dec_lai) - 2)]]
  
  # jan 
  jan_lai <- leaf[[grepl("_01", names(leaf))]]
  jan_lai <- jan_lai[[1:(nlyr(jan_lai) - 2)]]
  
  # feb  
  feb_lai <- leaf[[grepl("_02", names(leaf))]]
  feb_lai <- feb_lai[[1:(nlyr(feb_lai) - 2)]]
  
  # march
  mar_lai <- leaf[[grepl("_03", names(leaf))]]
  mar_lai <- mar_lai[[1:(nlyr(mar_lai) - 2)]]
  
}

# ---------------------
# function to combine all the data into a single dataframe for further analysis 
join_fcover_age <- function(age_combined, 
                            zone_raster, 
                            hru_rast, # hru rast
                            dem, # elevation 
                            
                            jan_lai, 
                            feb_lai,
                            mar_lai,
                            apr_lai, 
                            may_lai, 
                            june_lai, 
                            july_lai, 
                            aug_lai, 
                            sept_lai, 
                            oct_lai,
                            nov_lai, 
                            dec_lai,
                            
                            lcy_rast, 
                            fc, # als forest cover 
                            chm, # als chm (height)
                            
                            fc_ntems_files, # ntems file paths
                            h_ntems_files,  # ntems file paths 
                            
                            peak_fcov, # peak fcov from leaf 
                            output_prod_path,
                            years = years 
                            # h_chm # height from chm # remove for now
){
  
  # function to join age, fcover, and fcover mode data together 
  #   
  # inputs:
  #   @age_combined: raster with merges NTEMS layers for each year from 2015:2020
  #   @ zone_raster : static BEC ZONE raster with id
  #   @ hru_rast : raster with frst_cl attribute from hru polygon layer 
  #   @ dem : dem 
  #   @ [ ]_lai: monthly (or winter) lai from leaf
  #   @ lcy_rast : raster with last change year ID : unqique ID for each polygon 
  #   @ fc_ntems_files: ntems forest cover paths 
  #   @ h_ntems_files: ntems height file paths
  #   @ peak_fcov: july fcov from leaf 
  #   @ h_chm: height from CHM 
  # 
  # outputs 
  #   merged dataframe with age, x, y, and monthly fCOVER values 
  # 
  for (i in 1:nlyr(age_combined)){
    
    age_i <- age_combined[[i]]
    
    jan_i <- jan_lai[[i]]
    feb_i <- feb_lai[[i]]
    mar_i <- mar_lai[[i]]
    apr_i <- apr_lai[[i]]
    may_i <- may_lai[[i]]
    june_i <- june_lai[[i]]
    jul_i <- july_lai[[i]]
    aug_i <- aug_lai[[i]]
    sept_i <- sept_lai[[i]]
    oct_i <- oct_lai[[i]]
    nov_i <- nov_lai[[i]]
    dec_i <- dec_lai[[i]]
    
    fc_i <- rast(fc_ntems_files[i])
    fc_i <- fc_i/100 
    
    fc_i <- fc_i %>%
      terra::project(., apr_i, method='average')
    
    h_i <- rast(h_ntems_files[i])
    h_i <- h_i/1000 
    h_i <- h_i %>% 
      terra::project(., apr_i, method='average')
    
    # peak fcover from leaf 
    peak_fc_i <- peak_fcov[[1]] %>% 
      terra::project(., apr_i, method='average')
    
    data_i <- c(age_i, zone_raster, hru_rast, dem, lcy_rast, chm, fc)
    data_i$lai_grp <- 0 
    data_i$lai_grp[data_i$dem >= 1560] <- 1
    
    data_ii <- c(jan_i, feb_i, mar_i, apr_i, may_i, june_i, jul_i, 
                 aug_i, sept_i, oct_i, nov_i, dec_i, fc_i, h_i, peak_fc_i) %>% 
      terra::project(., data_i, method = "average")
    
    data_rast <- c(data_i, data_ii) %>% 
      writeRaster(., 
                  file.path(output_prod_path, paste0("veg_data_", 2014+i, ".tif")),
                  overwrite = TRUE)
    
    df_1 <- as.data.frame(data_i, xy=TRUE) %>%
      na.omit() 
    df_2 <- as.data.frame(data_ii, xy=TRUE)
    
    df_i <- merge(df_1, df_2, by=c("x", "y"))
    df_i$year <- years[i]
    
    names(df_i) <- c("x", "y", "age", "id", "frst_cl", "dem", "lcy_id", "2014_height", "2014_fc", "lai_grp",
                     "jan_lai", "feb_lai", "mar_lai", "april_lai", "may_lai", "june_lai", 
                     "july_lai", "aug_lai", "sept_lai", "oct_lai", "nov_i", "dec_i", "ntems_fcover", "ntems_height",
                     "peak_fcover_leaf", "year")
    
    if (i == 1){
      full_df <- df_i
    } else {
      full_df <- rbind(full_df, df_i)
    }
    
  }
  return(full_df)
}

# join data together 

full_df <- join_fcover_age(age_combined = age_combined, 
                           zone_raster = zone_raster,
                           hru_rast = hru_rast,
                           dem = dem,
                           jan_lai = jan_lai,
                           feb_lai = feb_lai, 
                           mar_lai = mar_lai,
                           apr_lai = apr_lai, 
                           may_lai = may_lai, 
                           june_lai = june_lai, 
                           july_lai = july_lai, 
                           aug_lai = aug_lai, 
                           sept_lai = sept_lai, 
                           oct_lai = oct_lai,
                           nov_lai = nov_lai,
                           dec_lai = dec_lai, 
                           lcy_rast = lcy_rast,
                           chm = chm, 
                           fc = fc,
                           fc_ntems_files = fc_ntems_files,
                           h_ntems_files = h_ntems_files, 
                           peak_fcov = peak_fcov, 
                           output_prod_path = output_prod_path, 
                           years = years)


# write fully merged data 
write.csv(full_df,
          file.path(output_prod_path, "veg_params_2015_2021_allmonths_feb19.csv"), # much smaller file size because downsampled to 
          row.names = FALSE)
# ------------------------------------------------------------------------------
# try sampling with limiting by same disturbance 
# ------------------------------------------------------------------------------

set.seed(1113)

# sample 25% of the rows in each df 
# --- df meant to be for each disturbance polygon
sample_25p <- function(df) {
  sampled_df <- df %>% sample_frac(0.25)
  return(sampled_df)
}

# ---------------------------------------------------
# sample full dataset 
# ---------------------------------------------------
# Apply the function per age group
sampled_df <- full_df %>%
  group_by(lcy_id) %>%
  group_modify(~ sample_25p(.x)) %>%
  ungroup()

samp_n <- sampled_df %>%
  group_by(age, lai_grp) %>%
  summarize(n = length(x))

# filtered to only forested < 100 years bc this is what we care about 
# RSE <0.07 for all and for majority < 0.05 
# originally just sampled 1000 points per age but missing a lot of data for IDF this way 
# because unequal area 
# so switching to trying to sample 1000 points, with same < 25% of one disturbance,
# in each BEC zone for each age, forested pixels only 
# sample 1000 points per year 

# ------------------------------------------------------------------------------
# sample by BEC zone
# ------------------------------------------------------------------------------
# 
# sampled_forested_rse <- sampled_df %>%
#   group_by(age, frst_cl) %>%
#   slice_sample(n = 1000) %>%
#   ungroup() %>%
#   filter(frst_cl %in% c("MS", "ESSF", "IDF")) %>%
#   filter(age < 100) %>%
#   group_by(age, frst_cl) %>% 
#   summarize(rse_fcover = AnglerCreelSurveySimulation::calculate_rse(na.omit(peak_fcover_leaf)),
#             rse_height = AnglerCreelSurveySimulation::calculate_rse(na.omit(ntems_height))
#   )
# 
# sampled_data <- full_df %>% # was sampled_df 
#   group_by(age, frst_cl) %>%
#   slice_sample(n = 1000) %>%
#   ungroup()
# 
# #p <- vect(sampled_data, geom = c("x", "y")) %>%
# #  st_as_sf()
# #st_crs(p) <- crs(leaf)
# 
# sampled_n <- sampled_data %>%
#   group_by(age, frst_cl) %>%
#   summarize(n = length(age))
# 
# 
# ########################################################
# # write out 1000 values for each yaer 
# write.csv(sampled_data,
#           file.path(output_prod_path, "sampled_veg_params_by_BECZONE_2015_2021_feb10.csv"),
#           row.names = FALSE)

# ------------------------------------------------------------------------------
# sample by elevation group
# ------------------------------------------------------------------------------

# subset to forested area
sampled_df <- full_df # rm if switching back
sampled_frst_df <- sampled_df[sampled_df$frst_cl %in% c("MS", "ESSF", "IDF"), ]

# print percent of TC that is forested as a check
# print(length(sampled_frst_df$lcy_id)/length(sampled_df$lcy_id)*100)

# check RSE for each age and LAI elevation group (aiming for < 0.05)
sampled_forested_rse <- sampled_frst_df %>%
  group_by(age, lai_grp) %>%
  slice_sample(n = 10000) %>%
  ungroup() %>%
  filter(age < 100) %>%
  group_by(age, lai_grp) %>% 
  summarize(rse_ntems_fcover = AnglerCreelSurveySimulation::calculate_rse(na.omit(ntems_fcover)),
            rse_ntems_height = AnglerCreelSurveySimulation::calculate_rse(na.omit(ntems_height)), 
            rse_peak_fcover_leaf = AnglerCreelSurveySimulation::calculate_rse(na.omit(peak_fcover_leaf)), 
            rse_july_lai = AnglerCreelSurveySimulation::calculate_rse(na.omit(july_lai))
  )

# sample 1000 points for each age in each lai group
sampled_frst_data <- sampled_frst_df %>% 
  group_by(age, lai_grp) %>%
  slice_sample(n = 10000) %>%
  ungroup()

# summarize the number of samples per age in each group (n)
sampled_frst_n <- sampled_frst_data %>%
  group_by(age, lai_grp) %>%
  summarize(n = length(age))


########################################################
# write out sampled forested data (10000 points per class)
write.csv(sampled_frst_data,
          file.path(output_prod_path, "sampled_veg_params_byLAI_GRP_2015_2021_feb19.csv"),
          row.names = FALSE)


# ------------------------------------------------------------------------------
# write out the other classes (WET_LAND, ALPINE, SHRUB) because I want this data too 
# to check LAI values 
# not going to downsample through because age doesn't matter - no recovery traj 

# subset to non forested classes
sampled_non_frst_df <- sampled_df[!sampled_df$frst_cl %in% c("MS", "ESSF", "IDF"), ]

# write 
write.csv(sampled_non_frst_df,
          file.path(output_prod_path, "sampled_veg_params_nonfrst_2015_2021_feb19.csv"),
          row.names = FALSE)
