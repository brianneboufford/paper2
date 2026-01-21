# ------------------------------------------------------------------------------
# adapted from 001_LAI_shift_values_sgsR.R
# 
#
# downsample data by CC class
# remove water and bare rock classes
# sample points with sgsR in classes
# take median to do shift because mroe resilient to outliers 
# 
#
# Novemember 13, 2025
# ------------------------------------------------------------------------------

# library
library(terra)
library(dplyr)
library(sf)
library(stringr)
library(ggplot2)
library(tidyr)
library(Metrics)
library(gridExtra)
library(sgsR)
library(AnglerCreelSurveySimulation)
library(purrr)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

setwd("C:/Users/blbouf/Sync/Paper2")

# output path
output_path <- file.path(".", "data", "fCOVER_hls_als_analysis")

# data paths 
lidar_fcover_path <- file.path(".", "data", "fCOVER_ALS", "percent_above_metrics.tiff")

cc_ntems_path <- file.path(".", "data", "src", "ntems", "structure.tif")

leaf_fcover_path <- file.path(".", "data", "fCOVER_clean", "fCOVER_stacked", "fcover_stacked_TC_area.tif")

vlce_path <- file.path(".", "data", "src", "ntems", "LC_Class_HMM_11S_2014_v20_v20.dat")

bec_path <- file.path(".", "data", "src", "BEC", "bec_zones_tc.shp")

# ------------------------------------------------------------------------------
# read data
#-------------------------------------------------------------------------------

# read lidar canopy cover
lidar_cc <- rast(lidar_fcover_path)$pzabove2
names(lidar_cc) <- "lidar_cc"

# read ntems canopy cover
cc_raw <- rast(cc_ntems_path)
cc_ntems_2014_raw <- cc_ntems$UTM_11S_percentage_first_returns_above_2m_2014
cc_ntems_2014 <- cc_ntems_2014_raw/100 # scaling factor provided by key 

# read land cover data
vlce_2014 <- rast(vlce_path)

# read BEC data
bec <- st_read(bec_path)
bec_r <- terra::rasterize(bec, lidar_cc, "ZONE")

# ------------------------------------------------------------------------------
# difference calculation 
# ------------------------------------------------------------------------------

# calculate monthly composite for LEAF cc 
leaf_rast <- rast(leaf_fcover_path)
leaf_comp_layers <- subset(leaf_rast, c("fCOVER_2014_08", "fCOVER_2014_09", "fCOVER_2014_10"))
                              
leaf_comp <- mean(leaf_comp_layers, na.rm = TRUE)
leaf_comp_proj <- project(leaf_comp, lidar_cc,
                          method = "average")
names(leaf_comp_proj) <- "leaf_fcover_mean"

# difference between lidar and hls 
cc_diff <- lidar_cc/100 - leaf_comp_proj
names(cc_diff) <- "diff_cc"

# project vlce, canopy cover, and overstory canopy cover to same crs 
vlce_proj <- terra::project(vlce_2014, cc_diff)
cc_proj <- terra::project(cc_ntems_2014, cc_diff)
#osc_proj <- terra::project(osc, cc_diff)

# separate cc into classes of 1-20%, 20-50%, 50-80%, >80%
cc_class_mat <- matrix(c(
  -Inf, 20, 0,
  20, 50, 1,
  50, 80, 2,
  80, Inf, 3
), ncol = 3, byrow = TRUE)

# Reclassify canopy cover raster 
cc_classified <- terra::classify(cc_proj, cc_class_mat)
cc_classified[is.na(cc_classified)] <- 0 
names(cc_classified) <- "cc_class"

# convert to df 
cc_df <- c(cc_diff, vlce_proj, cc_proj, lidar_cc/100, bec_r, leaf_comp_proj, cc_classified) %>%
  as.data.frame(xy = TRUE) #%>%
  #na.omit() # I think this drops all shrub and herbs because canopy cover is NA 

# total area for percent area calculation later 
total_area <- length(cc_df$x)

# remove water and rock/rubble class
cc_df <- cc_df[cc_df$category !=  "Water", ]
cc_df <- cc_df[cc_df$category != "Rock/Rubble", ]
cc_df <- cc_df[cc_df$category != "Exposed/Barren Land", ]

cc_df$scaling_factor <- cc_df$lidar_cc/cc_df$leaf_fcover_mean

# subest to coniferous ntems class only # CONFIEROUS #############################
cc_conif_df <- cc_df[cc_df$category == "Coniferous", ]
cc_conif_df <- cc_conif_df %>% na.omit()

# split coniferous into categories based on high, medium, and low canopy cover 
cc_conif_df$category_cc <- cc_conif_df$cc_class

# get canopy cover class info 
categor <- cc_conif_df$category
cc_class_conif <- cc_conif_df$cc_class

##3 update lia to cc
# make new column 
cc_conif_df$category_cc <- paste0(categor, "_", cc_class_conif)

# change category value to new canopy cover class and remove extra column 
cc_conif_df$category <- cc_conif_df$category_cc
cc_conif_df <- cc_conif_df %>%
  select(-c(category_cc))

# subset to non coniferous classes only # NOT CONIFEROUS #########################
cc_df_no_conif <- cc_df %>%
  subset(category != "Coniferous")

# rejoin all data so now categories are LC classes and CC classes for coniferous 
cc_df_recon <- rbind(cc_df_no_conif, cc_conif_df)

# ------------------------------------------------------------------------------
# sub sample data
# ------------------------------------------------------------------------------

set.seed(111)

cats <- unique(cc_df_recon$category)
samp_size_list <- c()

for (i in 1:length(cats)){
  categ <- cats[i]
  
  cc_cat <-cc_df_recon[cc_df_recon$category == categ, ]
  
  cc_cat_rast <- rast(cc_cat, type="xyz")
  crs(cc_cat_rast) <- crs(cc_diff)
  
  samp_size <- calculate_sampsize(cc_cat_rast$diff_cc)
  
  samp_size$category <- categ
  
  samp_size_n <- samp_size[c(1, nrow(samp_size)), ]
  samp_size_list <- rbind(samp_size_list, samp_size_n)
}

cc_recon_sample <- c()
rse_results <- c()
# sample size 
for (i in 1:length(cats)){
  
  # get category
  categ <- cats[i]
  
  # subset lai data to just this landcover category 
  cc_cat <- cc_df_recon[cc_df_recon$category == categ, ]
  
  # convert to raster 
  cc_cat_rast <- rast(cc_cat, type="xyz")
  crs(cc_cat_rast) <- crs(cc_diff)
  
  # sample 2000 points
  samp_df <- sample_srs(cc_cat_rast$diff_cc, 
                        nSamp = 792, # ideally would be 1200 bc of RSE < 0.05 for broadleaf but only 792 points so using this, all other classes are < this for RSE < 0.05 
                        mindist = 30)
  
  # sample from raster 
  samp_points_v <- terra::vect(samp_df)
  points_keep <- terra::cellFromXY(cc_cat_rast, terra::crds(samp_points_v))
  masked_cc <- cc_cat_rast * NA
  masked_cc[points_keep] <- cc_cat_rast[points_keep]
  
  vals_vect <- values(masked_cc$diff_cc)[!is.na(values(masked_cc$diff_cc))]
  
  rse <- calculate_rse(vals_vect)
  n <- length(vals_vect)
  
  rse_vals <- data.frame(category =  as.character(categ), RSE = rse, 
                         N = n, 
                         stringsAsFactors = FALSE)
  
  rse_results <- rbind(rse_results, rse_vals)
  
  # convert to dataframe and remove categories that turn NULL
  cc_df <- as.data.frame(masked_cc, xy = TRUE) %>% 
    select(-c("category", "ZONE"))
  
  # get categories back 
  cc_cat <- cc_cat %>% 
    select(c("x", "y", "category", "ZONE"))
  
  # join sampled data
  cc_df_updated <- left_join(cc_df, cc_cat, by=c("x", "y"))
  
  # merge onto final df 
  cc_recon_sample <- rbind(cc_recon_sample, cc_df_updated)
  
}

split_by_group <- function(df_group) {
  n <- nrow(df_group)
  train_indices <- sample(seq_len(n), size = 0.7 * n)
  df_group$split <- "test"
  df_group$split[train_indices] <- "train"
  return(df_group)
}

# Apply split by category
df_split <- cc_recon_sample %>%
  group_by(category) %>%
  group_split() %>%
  map_df(split_by_group) %>%
  as.data.frame()

cc_recon_sample <- df_split[df_split$split == "train", ]
cc_recon_hold <- df_split[df_split$split == "test", ]

write.csv(cc_recon_sample,
          file.path(output_path, "fCOVER_shift_train_sample.csv"),
          row.names = FALSE)

write.csv(cc_recon_hold,
          file.path(output_path, "fCOVER_shift_hold.csv"),
          row.names = FALSE)

cc_recon_sample <- read.csv(file.path(output_path, "fCOVER_shift_train_sample.csv"))
cc_recon_hold <- read.csv(file.path(output_path, "fCOVER_shift_hold.csv"))

# ------------------------------------------------------------------------------
# prep data for bar plot 
# ------------------------------------------------------------------------------

# summarize mean of leaf CC and lidar CC
cc_val_summary <- cc_recon_sample %>%
  group_by(category) %>%
  summarize(median_hls = median(leaf_fcover_mean, na.rm = TRUE),
            median_lidar = median(lidar_cc, na.rm = TRUE),
            mean_hls = mean(leaf_fcover_mean, na.rm = TRUE),
            mean_lidar = mean(lidar_cc, na.rm = TRUE)) %>%
  as.data.frame()

# pivot to format that is easier to plot as bar plot 
cc_val_summary_long <- pivot_longer(data = cc_val_summary, 
                                     cols = c("mean_hls", "mean_lidar", "median_lidar", "median_hls"),
                                     names_to = "data", 
                                     values_to = "cc_stat")

# --------------------------------
# landcover regression plot --- extra 
# --------------------------------

regression_plot <- ggplot(data=cc_recon_sample, aes(x = leaf_fcover_mean, y = lidar_cc)) + 
  geom_point(alpha=0.2) +
  theme_minimal() + 
  facet_wrap(~category)



# --------------------------------
# landcover summary statistics 
# - n 
# - mean 
# - sd 
# - rsq
# - pbias 
# --------------------------------

# summary with rounded values for mean and sd 
cc_summary_export <- cc_recon_sample %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean_diff_cc = round(mean(diff_cc),2), 
            sd_diff_cc = round(sd(diff_cc), 2), 
            median_diff_cc = round(median(diff_cc), 2)) %>%
  as.data.frame()

# calculate percent area 
percent_area <- cc_df_recon %>%
  group_by(category) %>%
  summarize(n = length(x))
percent_area$PA <- round(percent_area$n / total_area* 100, 2)

# merge together 
cc_summary_export <- merge(cc_summary_export, percent_area, by='category')
cc_summary_export <- cc_summary_export %>%
  select(c("category", "mean_diff_cc", "sd_diff_cc", "median_diff_cc", "PA"))

custom_order <- c("Shrubs", "Wetland",
                  "Wetland-Treed", "Herbs",
                  "Broadleaf", "Mixed Wood", 
                  "Coniferous_0",
                  "Coniferous_1", "Coniferous_2",
                  "Coniferous_3")

# Reorder dataframe
cc_summary_export <- cc_summary_export %>%
  mutate(category = factor(category, levels = custom_order)) %>%
  arrange(category)

# summary without rounded values 
# use this to actually shift values 
cc_summary <- cc_recon_sample %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean_diff_cc = mean(diff_cc), 
            sd_diff_cc = sd(diff_cc), 
            median_diff_cc= median(diff_cc)) %>%
  as.data.frame()

# join to percent area data for fun 
cc_summary <- merge(cc_summary, percent_area, by="category")

# write to output path
write.csv(cc_summary,
          file.path(output_path, "cc_residual_summary.csv"),
          row.names = FALSE)

# get R2 and pbias estimates for LAI sample 
cc_results <- cc_recon_sample %>%
  group_by(category) %>%
  summarise(
    R_squared = rsq(lidar_cc, leaf_fcover_mean),
    PBIAS = pbias(lidar_cc, leaf_fcover_mean),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# apply shift to held out HLS values 
# ------------------------------------------------------------------------------

# shift LAI values by adding mean LAI difference 
shifted_cc <- merge(cc_recon_hold, cc_summary, by = 'category')
shifted_cc$hls_cc_new <- shifted_cc$leaf_fcover_mean + shifted_cc$median_diff_cc

shifted_cc$hls_cc_new[shifted_cc$hls_cc_new > 1] <- 1

# take new difference 
shifted_cc$diff_shift <- shifted_cc$lidar_cc - shifted_cc$hls_cc_new

# summarize mean and sd of difference after shift
cc_shift_summary <- shifted_cc %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean = mean(diff_shift), 
            median = median(diff_shift),
            sd = sd(diff_shift)) %>%
  as.data.frame()

# Example: df with columns obs, pred, class_column
cc_shift_results <- shifted_cc %>%
  group_by(category) %>%
  summarise(
    R_squared_sh = rsq(lidar_cc, hls_cc_new),
    PBIAS_og = round(pbias(lidar_cc, leaf_fcover_mean),2),
    PBIAS_sh = round(pbias(lidar_cc, hls_cc_new),2),
    mean_lidar = round(mean(lidar_cc, na.rm = TRUE),2),
    mean_hls = round(mean(leaf_fcover_mean, ra.rm = TRUE),2),
    mean_adj = round(mean(hls_cc_new, na.rm = TRUE),2),
    .groups = "drop"
  )

# write to output path
write.csv(cc_shift_results,
          file.path(output_path, "cc_shift_mean_pbias_results.csv"),
          row.names = FALSE)

# convert back to map 
cc_shift <- shifted_cc %>% 
  select("x", "y", "hls_cc_new", "diff_shift") %>%
  rast(type="xyz", crs= crs(lidar_cc))

# ------------------------------------------------------------------------------

overall_results <- shifted_cc %>%
  group_by(category) %>%
  summarize(rmse0 = round(rmse(lidar_cc, leaf_fcover_mean),2),
            rmse1 = round(rmse(lidar_cc, hls_cc_new),2),
            pbias0 = round(pbias(lidar_cc, leaf_fcover_mean),2),
            pbias1 = round(pbias(lidar_cc, hls_cc_new),2))
# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------
# rsq function 
rsq <- function(obs, pred) {
  cor(obs, pred, use = "complete.obs")^2
}

# pbias function 
pbias <- function(obs, pred) {
  100 * sum(pred - obs, na.rm = TRUE) / sum(obs, na.rm = TRUE)
}
