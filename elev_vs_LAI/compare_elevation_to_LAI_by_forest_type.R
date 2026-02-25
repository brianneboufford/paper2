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
library(strucchange)
library(changepoint)

# can try change point package

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
  group_by(dem) %>%
  summarize(lai_av = mean(LAI))

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------

lai_elev_plot <- ggplot(results_grouped, aes(x = lai_av, y = dem)) +
  geom_point(alpha = 0.1) +
  geom_vline(xintercept = lai_cp, linetype = "dashed", 
             linewidth = 0.5, colour = "red") +
  labs(x = "LAI",
       y = "Elevation (m)",
       title = ""
  ) +
  theme_bw()

ggsave(lai_elev_plot, 
       filename = file.path("data", "figs", "lai_vs_elevation", "lai_elevation_plot.png"),
       dpi = 300, 
       units = "in",
       height = 8, 
       width = 8)

# ------------------------------------------------------------------------------
# sub sample data and detect breakpoint 
# -------------------------------------------------------------------------------

# subsample 10 000 points
df_samp <- dplyr::sample_n(results_grouped, 10000)

# order sampled and full dataset 
lai_samp_ordered <- df_samp[order(df_samp$dem), ]
lai_ordered <- results_grouped[order(results_grouped$dem), ]

lai_ordered$dem <- round(lai_ordered$dem, digits = 0)
#lai_ordered$lai_av <- round(lai_ordered$lai_av, digits = 1)

lai_ordered_av <- lai_ordered %>% 
  group_by(dem) %>%
  summarize(lai_av = mean(lai_av))

#lai_ordered_av <- lai_ordered_av[lai_ordered_av$dem > 1000, ]
#lai_ordered_av <- lai_ordered_av[lai_ordered_av$dem < 1850, ]

# get single breakpoint
cp_mean <- cpt.mean(lai_ordered_av$lai_av, method = "AMOC")
# this gives us the change point index of 675

lai_cp <- lai_ordered_av$lai_av[675]
elev_cp <- lai_ordered_av$dem[675]

# lai vs elevation
el_lai_cp_plot <- ggplot(lai_ordered_av, aes(x = dem, y = lai_av)) +
  geom_point(alpha = 0.4, size = 0.8) +                 
  geom_vline(xintercept = elev_cp, linetype = "dashed", 
             linewidth = 0.5, colour = "red") +
  labs(x = "Elevation (m)",
       y = "LAI",
       title = ""
  ) +
  theme_bw()

# flip coords
el_lai_cp_plot <- ggplot(lai_ordered_av, aes(x = lai_av, y = dem)) +
  geom_point(alpha = 0.4, size = 0.8) +                 
  geom_hline(yintercept = elev_cp, linetype = "dashed", 
             linewidth = 0.5, colour = "red") +
  labs(x = "LAI",
       y = "Elevation",
       title = ""
  ) +
  theme_bw()

lai_elev_samp_plot <- ggplot(df_samp, aes(x = lai_av, y = dem)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = elev_cp, linetype = "dashed", 
             linewidth = 0.5, colour = "red") +
  labs(x = "LAI",
       y = "Elevation (m)",
       title = ""
  ) +
  theme_bw()

ggsave(el_lai_cp_plot, 
       filename = file.path("data", "figs", "lai_vs_elevation", "cp_plot_feb20.png"),
       dpi = 300, 
       units = "in",
       height = 5, 
       width = 5)
# ------------------------------------------------------------------------------
# function to get elevation and LAI for each mature forest HRU  
# ------------------------------------------------------------------------------

summarize_lai_by_elevation <- function(yr, 
                                       dem,
                                       peak_lai, 
                                       hru_files, 
                                       hru_files_base_path){
  
  # peak LAI for year of interest
  peak_lai_i <- grepl(yr, names(peak_lai))
  p_lai <- peak_lai[[peak_lai_i]]
  
  # make mask layer
  lai_mask <- p_lai 
  lai_mask[!is.na(lai_mask)] <- 1
  
  # transform dem so it matches lai 
  dem_t <- resample(dem, p_lai)
  dem_c <- crop(dem_t, p_lai)
  dem_c <- dem_c*lai_mask
  ext(dem_c) <- ext(p_lai)
  
  # stack lai and dem
  rast_data <- c(p_lai, dem_c)
  
  # hru of year of interest 
  hru_i <- hru_files[grepl(yr, hru_files)]
  hru_path <- file.path(hru_files_base_path, hru_i)
  hru <- st_read(hru_path) %>% 
    st_transform(., crs(p_lai))
  
  # subset mature forested HRUs
  hru_forested <- hru[hru$VEG_CLA == "mature", ]
  
  # mask out non forested areas from the rast data 
  rast_data_mf <- mask(rast_data, hru_forested)
  
  rast_data_df <- as.data.frame(rast_data_mf)
  names(rast_data_df) <- c("LAI", "dem")
  rast_data_df$yr <- yr
  # hru_lai <- extract(p_lai, hru_forested, fun = mean)
  # hru_elev <- extract(dem_masked, hru_forested, fun = mean)
  
  # results <- cbind(hru_elev, hru_lai, hru_forested$frst_cl) %>% 
  #   select(-c("ID")) %>% 
  #   select(-c("ID"))
  # 
  # names(results) <- c("elevation", "LAI", "forest_type")
  
  return(rast_data_df)
  
}
