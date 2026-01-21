# -----------------------------------------------------------------------------_
# script to sample shifted fCOVER data by age and write sampled data sets 
#
#
# adapated from 00_prep_LAI_age_layers.R which is adapted from 
# adapted from 01_prep_bmu_recovery_layers 
# get fCOVER by age and pull out alpine, wetland, and shrub classes 
# 
# november 18, 2025
# adapted form june 3, 2025 version 
#
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)
library(future)
library(future.apply)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

output_prod_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "fCOVER_age")

# paths
if(TRUE){
  
  leaf_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "adjusted_hls_cc.tif")
  
  # dem 
  dem_path <- file.path(".", "data", "src", "lidar_derived", "dem.tif")
  
  # slope 
  slope_path <- file.path(".", "data", "src", "lidar_derived", "slope.tif")
  
  # aspect 
  aspect_path <- file.path(".", "data", "src", "lidar_derived", "aspect.tif")
  
  # bec_path 
  bec_path <- file.path(".", "data", "src", "BEC", "bec_zones_tc.shp")
  
  # hru path
  hru_path <- file.path(".", "data", "LCC_HRU_files", "Nov10", "HRUs_1923_2023", "HRU_2023.shp")
  
  # last change year polygon later 
  lcy_poly_path <- file.path(".", "data", "src", "ntems", "disturbance_year_polygon.shp")
  
}

if(TRUE){
  
  # fcover
  leaf <- rast(leaf_path)       
  leaf[leaf == 0] <- NA
  
  # dem
  dem <- rast(dem_path) %>%
    terra::resample(leaf, method='average') 
  
  # slope 
  slope_rast <- rast(slope_path) %>%
    terra::resample(leaf, method='average')
  
  # aspect
  aspect_rast <- rast(aspect_path) %>%
    terra::resample(leaf, method='average')
  
  # winter 
  winter <- leaf[[grepl("_12|_11|_01|_02|_03", names(leaf))]]
  
  # get winter mean for each year
  layer_info <- data.frame(
    name = names(winter),
    year = sub(".*_(\\d{4})_\\d{2}$", "\\1", names(winter))
  )
  
  
  # Unique years
  years <- unique(layer_info$year)
  
  age_dist <- rast(file.path(".", "data", "age", "age_disturbed.tif")) %>%
    terra::resample(leaf, method='mode')
  age_undist <- rast(file.path(".", "data", "age", "age_undisturbed.tif")) %>%
    terra::resample(leaf, method='mode')
  age_combined <- rast(file.path(".", "data", "age", "age_combined.tif")) %>%
    terra::resample(leaf, method='mode')
  zone_raster <- rast(file.path(".", "data", "src", "BEC", "zone_raster.tif")) %>%
    terra::resample(leaf, method='mode')
  
  # skipping for now
  hru_2023 <- st_read(hru_path) %>%
   st_transform(., crs = crs(age_combined))

  hru_rast <- rasterize(hru_2023, age_combined, field = "frst_cl")
  
  # end skipping 
  
  lcy_poly <- st_read(lcy_poly_path)
  lcy_poly$fid[lcy_poly$Last_CY == 0] <- 3333
  lcy_rast <- rasterize(lcy_poly, age_combined, field = "fid")
  
  
  # may 
  may <- leaf[[grepl("_05", names(leaf))]]
  may <- may[[1:(nlyr(may) - 2)]]
  
  # june 
  june <- leaf[[grepl("_06", names(leaf))]]
  june <- june[[1:(nlyr(june) - 2)]]
  
  # july 
  july <- leaf[[grepl("_07", names(leaf))]]
  july <- july[[1:(nlyr(july) - 2)]]
  
  # august 
  aug <- leaf[[grepl("_08", names(leaf))]]
  aug <- aug[[1:(nlyr(aug) - 2)]]
  
  # sept
  sept <- leaf[[grepl("_09", names(leaf))]]
  sept <- sept[[1:(nlyr(sept) - 2)]]
  
  # oct
  oct <- leaf[[grepl("_10", names(leaf))]]
  oct <- oct[[1:(nlyr(oct) - 2)]]
  
  oct_q <- global(oct, probs = (0.90), na.rm = TRUE)
  max_wint <- max(oct_q$mean)
  
  # april
  apr <- leaf[[grepl("_04", names(leaf))]]
  apr <- apr[[1:(nlyr(apr) - 2)]]
  apr[apr > max_wint] <- NA
  
  # for supplementary material ------------------------------------------------
  # NA assessment of winter layers 
  # valid_mask <- !is.na(values(hru_rast))
  # 
  # # Function to compute %NA for valid pixels only
  # na_percent <- function(layer_i) {
  #   vals <- values(layer_i)[valid_mask]
  #   100 * sum(is.na(vals)) / length(vals)
  # }
  # 
  # # Apply to each layer in the stack
  # na_summary <- sapply(1:nlyr(winter_lai), function(i) na_percent(winter_lai[[i]]))
  # names(na_summary) <- names(winter_lai)
  # na_summary <- na_summary %>% as.data.frame()
  # names(na_summary) <- c("% NA")
  # na_summary$date_int <- 1:nrow(na_summary)
  # 
  
  #write.csv(na_summary, 
  #          file.path(data_path, "..", "compare_LiDAR_HLS_LAI", 
  #               "LAI_corrected_may28","supp_material", "percent_NA_winter.csv"),
  #          row.names = TRUE)
  # ----------------------------------------------------------------------------
  # Compute mean per year
  annual_means <- lapply(years, function(yr) {
    # get layer names from year of interest
    layer_names <- layer_info$name[layer_info$year == yr]
    # pull out of winter stack
    yearly_stack <- winter[[layer_names]]
    
    yearly_stack[yearly_stack > max_wint] <- NA
    # calculate mean 
    yearly_mean <- terra::mean(yearly_stack, na.rm = TRUE)
    # rename 
    names(yearly_mean) <- paste0("mean_winter_fCOVER_", yr)
    return(yearly_mean)
    
  }) %>% 
    c()
  
  # compute median
  annual_median <- lapply(years, function(yr) {
    # get layer names from year of interset
    layer_names <- layer_info$name[layer_info$year == yr]
    # pull out of winter stack 
    yearly_stack <- winter[[layer_names]]
    
    yearly_stack[yearly_stack > max_wint] <- NA
    # calculate median
    yearly_med <- terra::median(yearly_stack, na.rm = TRUE)
    # rename
    names(yearly_med) <- paste0("med_winter_fCOVER_", yr)
    return(yearly_med)
  })
  
  # Combine into one SpatRaster
  mean_winter <- rast(annual_means)
  mean_winter <- mean_winter[[2:(nlyr(mean_winter)) - 2]]
  
  med_winter <- rast(annual_median)
  med_winter <- med_winter[[2:(nlyr(med_winter)) - 2]]
  
}

# ---------------------
# function to combine all the data into a single dataframe for further analysis 
join_fcover_age <- function(age_combined, zone_raster, hru_rast,
                         dem, slope_rast, aspect_rast, lcy_rast, mean_winter, 
                         med_winter, apr, may, june, july, 
                         aug, sept, oct){
  
  # function to join age, fcover, and fcover mode data together 
  #   
  # inputs:
  #   @age_combined: raster with merges NTEMS layers for each year from 2015:2020
  #   @ zone_raster : static BEC ZONE raster with id
  #   @ hru_rast : raster with frst_cl attribute from hru polygon layer 
  #   @ dem : dem 
  #   @ slope : slope  
  #   @ aspect : aspect 
  #   @ lcy_rast : raster with last change year ID : unqique ID for each polygon 
  #   @ mean_winter: mean winter fCOVER for given year
  #   @ med_winter: median winter fCOVER for given year
  #   @ [month] : raster layer with monthly fCOVER for 2015:2020
  # 
  # outputs 
  #   merged dataframe with age, x, y, and monthly fCOVER values 
  # 
  for (i in 1:nlyr(age_combined)){
    
    age_i <- age_combined[[i]]
    mean_winter_i <- mean_winter[[i]]
    med_winter_i <- med_winter[[i]]
    apr_i <- apr[[i]]
    may_i <- may[[i]]
    june_i <- june[[i]]
    jul_i <- july[[i]]
    aug_i <- aug[[i]]
    sept_i <- sept[[i]]
    oct_i <- oct[[i]]
    
    data_i <- c(age_i, zone_raster, hru_rast, dem, slope_rast, aspect_rast, lcy_rast) 
    data_ii <- c(mean_winter_i, med_winter_i, apr_i, may_i, june_i, jul_i, 
                 aug_i, sept_i, oct_i)
    
    df_1 <- as.data.frame(data_i, xy=TRUE) %>%
      na.omit() 
    df_2 <- as.data.frame(data_ii, xy=TRUE)
    
    df_i <- merge(df_1, df_2, by=c("x", "y"))
    
    names(df_i) <- c("x", "y", "age", "id", "frst_cl", "dem", "slope", "aspect", "lcy_id",
                     "mean_winter_fCOVER", "med_winter_fCOVER", "april_fCOVER", "may_fCOVER", "june_fCOVER", 
                     "july_fCOVER", "aug_fCOVER", "sept_fCOVER", "oct_fCOVER")
    
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
                        slope = slope_rast,
                        aspect = aspect_rast,
                        lcy_rast = lcy_rast,
                        mean_winter = mean_winter,
                        med_winter = med_winter,
                        apr = apr,
                        may = may, 
                        june = june, 
                        july = july, 
                        aug = aug, 
                        sept = sept,
                        oct = oct)

dist_df <- join_fcover_age(age_combined = age_dist, 
                        zone_raster = zone_raster,
                        hru_rast = hru_rast,
                        dem = dem,
                        slope = slope_rast,
                        aspect = aspect_rast,
                        lcy_rast = lcy_rast,
                        mean_winter = mean_winter,
                        med_winter = med_winter,
                        apr = apr,
                        may = may, 
                        june = june, 
                        july = july, 
                        aug = aug, 
                        sept = sept,
                        oct = oct)


undist_df <- join_fcover_age(age_combined = age_undist, 
                          zone_raster = zone_raster,
                          hru_rast = hru_rast,
                          dem = dem,
                          slope = slope_rast,
                          aspect = aspect_rast,
                          lcy_rast = lcy_rast,
                          mean_winter = mean_winter,
                          med_winter = med_winter,
                          apr = apr,
                          may = may, 
                          june = june, 
                          july = july, 
                          aug = aug, 
                          sept = sept,
                          oct = oct)

# write fully merged data 
write.csv(full_df,
          file.path(output_prod_path, "fCOVER_all_data_2014_2020_nov18.csv"),
          row.names = FALSE)

#write.csv(dist_df,
#          file.path(output_prod_path, "recovery_results", "lai_distrub_2014_2020_may29.csv"),
#          row.names = FALSE)

write.csv(undist_df,
          file.path(output_prod_path, "fCOVER_undistrub_2014_2020_nov18.csv"),
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
  group_by(age) %>%
  summarize(n = length(x))

# filtered to only forested < 100 years bc this is what we care about 
# RSE <0.07 for all and for majority < 0.05 
# originally just sampled 1000 points per age but missing a lot of data for IDF this way 
# because unequal area 
# so switching to trying to sample 1000 points, with same < 25% of one disturbance,
# in each BEC zone for each age, forested pixels only 
# sample 1000 points per year 
sampled_forested_rse <- sampled_df %>%
  group_by(age, frst_cl) %>%
  slice_sample(n = 1000) %>%
  ungroup() %>%
  filter(frst_cl %in% c("MS", "ESSF", "IDF")) %>%
  filter(age < 100) %>%
  group_by(age, frst_cl) %>% 
  summarize(rse_july = AnglerCreelSurveySimulation::calculate_rse(na.omit(july_fCOVER)))

sampled_data <- sampled_df %>% 
  group_by(age, frst_cl) %>%
  slice_sample(n = 1000) %>%
  ungroup()

#p <- vect(sampled_data, geom = c("x", "y")) %>%
#  st_as_sf()
#st_crs(p) <- crs(leaf)

sampled_n <- sampled_data %>%
  group_by(age, frst_cl) %>%
  summarize(n = length(age))

# ---------------------------------------------------
# sample undisturbed data 
# ---------------------------------------------------

# Apply the function per age group
undist_sampled_df <- undist_df %>%
  group_by(lcy_id) %>%
  group_modify(~ sample_25p(.x)) %>%
  ungroup()

# sample 1000 points per year 
undist_sampled_data <- undist_sampled_df %>%
  group_by(age) %>%
  slice_sample(n = 1000) %>%
  ungroup() 

undist_sampled_n <- undist_sampled_data %>%
  group_by(age) %>%
  summarize(n = length(x))

# dist_sample <- dist_df %>%
#   group_by(id, age) %>%
#   slice_sample(n = 1000) %>%
#   ungroup()
# 
# undist_sample <- undist_df %>%
#   group_by(id, age) %>%
#   slice_sample(n = 1000) %>%
#   ungroup()

########################################################
# write out 1000 values for each yaer 
write.csv(sampled_data,
          file.path(output_prod_path, "sampled_fCOVER_all_data_2014_2021_nov19.csv"),
          row.names = FALSE)

write.csv(sampled_forested_rse, 
          file.path(output_prod_path, "sampled_fCOVER_forested_RSE_nov19.csv"),
          row.names = FALSE)
# write.csv(dist_sample,
#           file.path(output_prod_path, "recovery_results", "sampled_lai_disturb_2015_2020_may15.csv"),
#           row.names = FALSE)
write.csv(undist_sampled_data,
          file.path(output_prod_path, "sampled_fCOVER_undisturb_2014_2021_nov19.csv"),
          row.names = FALSE)
