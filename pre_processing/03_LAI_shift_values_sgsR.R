# ------------------------------------------------------------------------------
# tidied residual analysis to set up for shifting rasters  
# encoportating canopy cover classes into coniferous class 
###############################################
# based on LAI_boxplot_analysis.R exploratory analysis and LAI_residual_shift_values.R
# based on 01_LAI_residual_shift_values.R
#
# downsample data by class
# remove water and bare rock classes
# sample points with sgsR in classes
# take median to do shift because mroe resilient to outliers 
# 
###############################################
#
# JUne 3th, 2025
# updated march 22 for paper 2
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
# set working directory to the location of the github folder
setwd(file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", "scripts", "LAI_recovery"))

# output path 
output_path <- file.path("..", "..", "data", "compare_LiDAR_HLS_LAI", "LAI_corrected_march22")

# main project path 
project_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek")

# lidar lai clean 
#lidar_lai_clean_path <- file.path(output_path, "..", "lidar_lai.tif") # previous version used in paper 1
lidar_lai_clean_path <- file.path(output_path, "..", "lai_clean_mar20.tif")

# ntems path 
ntems_path <- file.path(output_path, "..", "from_bud", "catchments (1)", "11S")
cc_path <- file.path(ntems_path, "structure", "percentage_first_returns_above_2m")
osc_path <- file.path(ntems_path, "structure", "percentage_first_returns_above_mean")

# hls lai clean (three_month_av_clean)
hls_lai_clean_path <- file.path(output_path, "..", "three_month_av_2014.tif")

# dem 
dem_path <- file.path("..", "..", "..", "data", "lidar", "output", "dem", "dem.tif")

# slope 
slope_path <- file.path("..", "..", "..", "data", "lidar", "output", "dem", "slope.tif")

# aspect 
aspect_path <- file.path("..", "..", "..", "data", "lidar", "output", "dem", "aspect.tif")

# path to VLCE 2.0 data 
vlce_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek", "data", 
                       "raw", "11S", "VLCE2.0")
vlce_files <- list.files(vlce_path, full.names = TRUE)

bec_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek", "data", "clean", 
                      "bec_TC", "bec_tc_clipped.shp")
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

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# read in ntems structure data 
cc_raw <- list.files(cc_path, pattern = "2014.dat$", full.names = TRUE) %>%
  rast()
cc <- cc_raw/100 # scaling factor provided by key 
osc_raw <- list.files(osc_path, pattern = "2014.dat$", full.names = TRUE) %>%
  rast()
osc <- osc_raw/100 # scaling factor provided by key 

# read lidar lai data
lidar_lai <- rast(lidar_lai_clean_path)
names(lidar_lai) <- "lidar_lai"

# read hls lai data and project to the same crs  
hls_lai <- rast(hls_lai_clean_path)
hls_lai_proj <- terra::project(hls_lai, lidar_lai)
names(hls_lai_proj) <- "hls_lai"

# read land cover data
vlce <- vlce_files[grepl("_2014_v20_v20.dat$", vlce_files)] %>% rast()
bec <- st_read(bec_path)
bec_r <- terra::rasterize(bec, lidar_lai, "ZONE")

# difference between lidar and hls 
lai_diff <- lidar_lai - hls_lai_proj
names(lai_diff) <- "diff_lai"

# project vlce, canopy cover, and overstory canopy cover to same crs 
vlce_proj <- terra::project(vlce, lai_diff)
cc_proj <- terra::project(cc, lai_diff)
osc_proj <- terra::project(osc, lai_diff)

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
lai_df <- c(lai_diff, vlce_proj, hls_lai_proj, lidar_lai, bec_r, cc_classified) %>%
  as.data.frame(xy = TRUE) %>%
  na.omit()

# total area for percent area calculation later 
total_area <- length(lai_df$x)

# remove water and rock/rubble class
lai_df <- lai_df[lai_df$category !=  "Water", ]
lai_df <- lai_df[lai_df$category != "Rock/Rubble", ]
lai_df <- lai_df[lai_df$category != "Exposed/Barren Land", ]

lai_df$scaling_factor <- lai_df$lidar_lai/lai_df$hls_lai

# subest to coniferous ntems class only # CONFIEROUS #############################
cc_lai_df <- lai_df[lai_df$category == "Coniferous", ]

# split coniferous and mixed wood into categories based on high, medium, and low canopy cover 
cc_lai_df$category_cc <- cc_lai_df$cc_class

# get canopy cover class info 
categor <- cc_lai_df$category
cc_class_conif <- cc_lai_df$cc_class

# make new column 
cc_lai_df$category_cc <- paste0(categor, "_", cc_class_conif)

# change category value to new canopy cover class and remove extra column 
cc_lai_df$category <- cc_lai_df$category_cc
cc_lai_df <- cc_lai_df %>%
  dplyr::select(-c(category_cc))

# subset to non coniferous classes only # NOT CONIFEROUS #########################
lai_df_no_conif <- lai_df %>%
  subset(category != "Coniferous")

# rejoin all data so now categories are LC classes and CC classes for coniferous 
lai_df_recon <- rbind(lai_df_no_conif, cc_lai_df)

# ------------------------------------------------------------------------------
# sub sample data
# ------------------------------------------------------------------------------

set.seed(111)

cats <- unique(lai_df_recon$category)
samp_size_list <- c()
  
for (i in 1:length(cats)){
  
  categ <- as.character(cats[i])
  
  lai_cat <-lai_df_recon[lai_df_recon$category == categ, ]
  
  lai_cat_rast <- rast(lai_cat, type="xyz")
  
  crs(lai_cat_rast) <- crs(lai_diff)
  
  samp_size <- calculate_sampsize(lai_cat_rast$diff_lai)$nSamp %>%
    unlist() %>%
    as.numeric()
  
  samp_df <- as.data.frame(t(samp_size))
  samp_df$category <- categ
  samp_size_list <- rbind(samp_size_list, samp_df)

}

# convert to dataframe and write 
# samp_size_df <- samp_size_list %>% as.data.frame()
# write.csv(samp_size_df,
#           file.path("..","..", "data", "compare_LiDAR_HLS_LAI", "LAI_corrected_may28",
#                     "supp_material", "samp_size_sgsR.csv"),
#           row.names = FALSE)

lai_recon_sample <- c()
rse_results <- c()
# sample size 
for (i in 1:length(cats)){
  
  # get category
  categ <- cats[i]
  
  # subset lai data to just this landcover category 
  lai_cat <-lai_df_recon[lai_df_recon$category == categ, ]
  
  # convert to raster 
  lai_cat_rast <- rast(lai_cat, type="xyz")
  crs(lai_cat_rast) <- crs(lai_diff)
  
  # sample 2000 points
  samp_df <- sample_srs(lai_cat_rast$diff_lai, 
             nSamp = 2000, 
             mindist = 20)

  # sample from raster 
  samp_points_v <- terra::vect(samp_df)
  points_keep <- terra::cellFromXY(lai_cat_rast, terra::crds(samp_points_v))
  masked_lai <- lai_cat_rast * NA
  masked_lai[points_keep] <- lai_cat_rast[points_keep]
  
  vals_vect <- values(masked_lai$diff_lai)[!is.na(values(masked_lai$diff_lai))]

  rse <- calculate_rse(vals_vect)
  n <- length(vals_vect)
  
  rse_vals <- data.frame(category =  as.character(categ), RSE = rse, 
                         N = n, 
                         stringsAsFactors = FALSE)
  
  rse_results <- rbind(rse_results, rse_vals)
  
  # convert to dataframe and remove categories that turn NULL
  lai_df <- as.data.frame(masked_lai, xy = TRUE) %>% 
    dplyr::select(-c("category", "ZONE"))
  
  # get categories back 
  lai_cat <- lai_cat %>% 
    dplyr::select(c("x", "y", "category", "ZONE"))
  
  # join sampled data
  lai_df_updated <- left_join(lai_df, lai_cat, by=c("x", "y"))
  
  # merge onto final df 
  lai_recon_sample <- rbind(lai_recon_sample, lai_df_updated)
  
}

# write.csv(rse_results,
#           file.path(output_path, "supp_material", "RSE_results.csv"),
#           row.names = FALSE)

# split into 70 30 test train
# Grouped split function

split_by_group <- function(df_group) {
  n <- nrow(df_group)
  train_indices <- sample(seq_len(n), size = 0.7 * n)
  df_group$split <- "test"
  df_group$split[train_indices] <- "train"
  return(df_group)
}

# Apply split by category
df_split <- lai_recon_sample %>%
  group_by(category) %>%
  group_split() %>%
  map_df(split_by_group) %>%
  as.data.frame()

lai_recon_sample <- df_split[df_split$split == "train", ]
lai_recon_hold <- df_split[df_split$split == "test", ]
# 
write.csv(lai_recon_sample,
          file.path(output_path, "LAI_shift_train_sample_mar22.csv"),
          row.names = FALSE)

write.csv(lai_recon_hold,
          file.path(output_path, "LAI_shift_hold_mar22.csv"),
          row.names = FALSE)

lai_recon_sample <- read.csv(file.path(output_path, "LAI_shift_train_sample_mar22.csv"))
lai_recon_hold <- read.csv(file.path(output_path, "LAI_shift_hold_mar22.csv"))

# ------------------------------------------------------------------------------
# prep data for bar plot 
# ------------------------------------------------------------------------------

# summarize mean of HLS LAI and lidar LAI
lai_val_summary <- lai_recon_sample %>%
  group_by(category) %>%
  summarize(median_hls = median(hls_lai, na.rm = TRUE),
            median_lidar = median(lidar_lai, na.rm = TRUE),
            mean_hls = mean(hls_lai, na.rm = TRUE),
            mean_lidar = mean(lidar_lai, na.rm = TRUE)) %>%
  as.data.frame()

# pivot to format that is easier to plot as bar plot 
lai_val_summary_long <- pivot_longer(data = lai_val_summary, 
                                     cols = c("mean_hls", "mean_lidar", "median_lidar", "median_hls"),
                                     names_to = "data", 
                                     values_to = "lai_stat")

# --------------------------------
# landcover regression plot --- extra 
# --------------------------------

regression_plot <- ggplot(data=lai_recon_sample, aes(x = hls_lai, y = lidar_lai)) + 
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
lai_summary_export <- lai_recon_sample %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean_diff_lai = round(mean(diff_lai),2), 
            sd_diff_lai = round(sd(diff_lai), 2), 
            median_diff_lai = round(median(scaling_factor), 2)) %>%
  as.data.frame()

# caluclate percent area 
percent_area <- lai_df_recon %>%
  group_by(category) %>%
  summarize(n = length(x))
percent_area$PA <- round(percent_area$n / total_area* 100, 2)

# merge together 
lai_summary_export <- merge(lai_summary_export, percent_area, by='category')
lai_summary_export <- lai_summary_export %>%
  dplyr::select(c("category", "mean_diff_lai", "sd_diff_lai", "median_diff_lai", "PA"))

custom_order <- c("Shrubs", "Wetland",
                  "Wetland-Treed", "Herbs",
                  "Broadleaf", "Mixed Wood", 
                  "Coniferous_0",
                  "Coniferous_1", "Coniferous_2",
                  "Coniferous_3")

# Reorder dataframe
lai_summary_export <- lai_summary_export %>%
  mutate(category = factor(category, levels = custom_order)) %>%
  arrange(category)

# summary without rounded values 
# use this to actually shift values 
lai_summary <- lai_recon_sample %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean_diff_lai = mean(diff_lai), 
            sd_diff_lai = sd(diff_lai), 
            median_diff_lai= median(diff_lai)) %>%
  as.data.frame()

# join to percent area data for fun 
lai_summary <- merge(lai_summary, percent_area, by="category")

# write to output path
write.csv(lai_summary,
          file.path(output_path, "lai_residual_summary_mar22.csv"),
          row.names = FALSE)

# get R2 and pbias estimates for LAI sample 
lai_results <- lai_recon_sample %>%
  group_by(category) %>%
  summarise(
    R_squared = rsq(lidar_lai, hls_lai),
    PBIAS = pbias(lidar_lai, hls_lai),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# apply shift to held out HLS values 
# ------------------------------------------------------------------------------

# shift LAI values by adding mean LAI difference 
shifted_lai <- merge(lai_recon_hold, lai_summary, by = 'category')
shifted_lai$hls_lai_new <- shifted_lai$hls_lai + shifted_lai$median_diff_lai

# take new difference 
shifted_lai$diff_shift <- shifted_lai$lidar_lai - shifted_lai$hls_lai_new

# summarize mean and sd of difference after shift
lai_shift_summary <- shifted_lai %>%
  group_by(category) %>%
  summarize(n = length(x), 
            mean = mean(diff_shift), 
            median = median(diff_shift),
            sd = sd(diff_shift)) %>%
  as.data.frame()

# Example: df with columns obs, pred, class_column
lai_shift_results <- shifted_lai %>%
  group_by(category) %>%
  summarise(
    R_squared_sh = rsq(lidar_lai, hls_lai_new),
    PBIAS_og = round(pbias(lidar_lai, hls_lai),2),
    PBIAS_sh = round(pbias(lidar_lai, hls_lai_new),2),
    mean_lidar = round(mean(lidar_lai, na.rm = TRUE),2),
    mean_hls = round(mean(hls_lai, ra.rm = TRUE),2),
    mean_adj = round(mean(hls_lai_new, na.rm = TRUE),2),
    .groups = "drop"
  )

# write to output path
write.csv(lai_shift_results,
          file.path(output_path, "lai_shift_mean_pbias_results_mar22.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# NOTES: 
# adding a scale factor doesn't change the R2 and the PBIAS is worse so it is better
# to just add the mean difference to the LAI values 
# ------------------------------------------------------------------------------

# convert back to map 
lai_shift <- shifted_lai %>% 
  dplyr::select("x", "y", "hls_lai_new", "diff_shift") %>%
  rast(type="xyz", crs= crs(lidar_lai))

# ------------------------------------------------------------------------------

overall_results <- shifted_lai %>%
  group_by(category) %>%
  summarize(rmse0 = round(rmse(lidar_lai, hls_lai),2),
            rmse1 = round(rmse(lidar_lai, hls_lai_new),2),
            pbias0 = round(pbias(lidar_lai, hls_lai),2),
            pbias1 = round(pbias(lidar_lai, hls_lai_new),2))

# ----------------------------------------------------------------------------
# went up to here for march 22 update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ----------------------------------------------------------------------------

# take out EXPOSED LAND --> made it worse 

library(knitr)
library(kableExtra)

overall_results$rmse_diff <- overall_results$rmse1 - overall_results$rmse0
overall_results$pbias_diff <- overall_results$pbias1 - overall_results$pbias0
kable(overall_results, format = "html") %>%
  kable_styling()

kable(lai_summary_export, format = "html") %>%
  kable_styling()


# ------------------------------------------------------------------------------
# 4 panel grid plot 
# ------------------------------------------------------------------------------
y_tick_size <- 12
x_tick_size <- 12
y_title_size <- 14
x_title_size <- 14
legend_text_size <- 12

lai_no_conif <- lai_recon_sample[!lai_recon_sample$category %in% c("Coniferous_0",
                                                                   "Coniferous_1",
                                                                   "Coniferous_2",
                                                                   "Coniferous_3"), ]

lai_conif <- lai_recon_sample[lai_recon_sample$category %in% c("Coniferous_0",
                                                               "Coniferous_1",
                                                               "Coniferous_2",
                                                               "Coniferous_3"), ]

lai_val_summary_long_noconif <- lai_val_summary_long[!lai_val_summary_long$category %in% c("Coniferous_0",
                                                                                           "Coniferous_1",
                                                                                           "Coniferous_2",
                                                                                           "Coniferous_3"), ]

lai_val_summary_long_conif <- lai_val_summary_long[lai_val_summary_long$category %in% c("Coniferous_0",
                                                                                        "Coniferous_1",
                                                                                        "Coniferous_2",
                                                                                        "Coniferous_3"), ]


#### SUBPLOT
# box plot of difference between LAI values based on NTEMS landcover class
lc_boxplot <- ggplot(data = lai_no_conif, aes(x = category, y = diff_lai, fill = category)) + 
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = c("#fee090" , "#01665e", 
                               "#80cdc1", "#b8e186", "#276419", "#7fbc41"))+
  xlab("Landcover class") +
  ylab("ALS LAI - LEAF LAI") +
  theme_classic(14) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=y_title_size),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=y_tick_size), 
        legend.position = "none") + 
  stat_summary(fun = mean, geom = "point", shape = 8, size = 3, color = "black") +
  scale_y_continuous(limits=c(-2, 4), breaks=c(seq(-2, 4, 2))) 

#### SUBPLOT 
lc_barplot <- ggplot(data = lai_val_summary_long_noconif, aes(x = as.factor(category), y = mean_lai, fill = data)) + 
  geom_bar(stat="identity", position = "dodge", width = 0.5) +
  scale_fill_manual(labels = c("LEAF", "ALS"), values = c("darkgrey", "grey4")) +
  scale_x_discrete(labels = c("Shrubs" = "Shrub",
                              "Wetland" = "WL",
                              "Wetland-Treed" = "WL-Treed",
                              "Herbs" = "Herb",
                              "Broadlead" = "BL", 
                              "Mixed Wood" = "MW")) +
  xlab("Landcover Class") +
  ylab("Mean LAI") +
  theme_classic(14) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.text.y = element_text(size = y_tick_size),
        axis.title.y = element_text(size = y_title_size),
        axis.title.x = element_text(size = x_title_size),
        axis.text.x = element_text(size = x_tick_size, angle = 0, vjust=0),
        legend.position = c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text_size)) + 
  scale_y_continuous(limits=c(0, 4), breaks=c(0:4)) 


#### SUBPLOT 
conif_order <- c("Coniferous_0", "Coniferous_1", "Coniferous_2", "Coniferous_3")

lai_conif <- lai_conif %>% 
  mutate(category = factor(category, levels = conif_order)) %>% 
  arrange(category)

lai_val_summary_long_conif <- lai_val_summary_long_conif %>% 
  mutate(category = factor(category, levels = conif_order)) %>% 
  arrange(category)

# box plot of difference between LAI values based on NTEMS landcover class
cc_boxplot <- ggplot(data = lai_conif, aes(x = category, y = diff_lai, fill = category)) + 
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  scale_fill_manual(values = c("#d9f0d3", "#5aae61", "#1b7837", "#00441b")) + 
  xlab("Landcover class") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") +
  theme_classic(14) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.text.x = element_blank(),
        axis.text.y =  element_blank(), 
        legend.position = "none") + 
  stat_summary(fun = mean, geom = "point", shape = 8, size = 3, color = "black") + 
  scale_y_continuous(limits=c(-2, 4), breaks=c(seq(-2, 4, 2))) 

#### SUBPLOT
cc_barplot <- ggplot(data = lai_val_summary_long_conif, aes(x = as.factor(category), y = mean_lai, fill = data)) + 
  geom_bar(stat="identity", position = "dodge", width = 0.5) +
  scale_fill_manual(labels = c("LEAF", "ALS"), values = c("darkgrey", "grey4")) +
  xlab("Coniferous Canopy Cover Class") +
  ylab("") +
  theme_classic(14) +
  scale_x_discrete(labels = c("Coniferous_0" = "<20%", "Coniferous_1" ="20-50%", "Coniferous_2" ="50-80%", "Coniferous_3" = "80%")) + 
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.title.x = element_text(size = x_title_size),
        axis.text.x = element_text(size = x_tick_size, angle = 0, vjust=0.5),
        axis.text.y =  element_blank(), 
        legend.position = "none",
        legend.title = element_blank()) + 
  scale_y_continuous(limits=c(0, 4), breaks=c(0:4)) 

res_grid <- cowplot::plot_grid(lc_boxplot, NULL, cc_boxplot,    
                               lc_barplot, NULL, cc_barplot,
                               ncol = 3, rel_widths = c(3, -0.1,  2, 
                                                        3, -0.1, 2),
                               align="hv")

# grid plot with frou panels
ggsave(file.path(output_path, "results", "four_panel_lai_resid_june3.png"),
       res_grid, 
       units="in", 
       width=10, 
       height=7, 
       dpi=500)

# ------------------------------------------------------------------------------
lai_all <- shifted_lai
lai_all$group <- ifelse(shifted_lai$diff_lai >= 0, "Positive", "Negative")

conif.labs <- c("Shrubs", "Wetland", "Wetland-Treed", "Herbs", "Broadleaf", "Mixed Wood", 
                "< 20% CC Conif", "20-50% CC Conif", "50-80% CC Conif", "> 80% CC Conif")
names(conif.labs) <- c("Shrubs", "Wetland", "Wetland-Treed", "Herbs", "Broadleaf", "Mixed Wood", 
                       "Coniferous_0", "Coniferous_1", "Coniferous_2", "Coniferous_3")

lai_all <- lai_all %>% 
  mutate(category = factor(category, levels = names(conif.labs))) %>% 
  arrange(category)

p4 <- ggplot(data=lai_all, aes(x=lidar_lai, diff_lai, color = group)) + 
  geom_point(alpha=0.5, size=0.3) +
  facet_wrap(~category, labeller = labeller(category = conif.labs), nrow=4) + 
  geom_smooth(method = "lm", se= TRUE, linewidth = 0.8, color="black") + 
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  theme_classic(14) +
  theme(axis.title = element_text(size = 14), 
        title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        legend.position = "none") + 
  labs(x = 'ALS LAI', y = "(ALS LAI - LEAF LAI)") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10), labels = c(0,2,4,6,8,10)) + 
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), color = "black", geom = "label")


# grid plot with frou panels
ggsave(file.path(output_path, "results", "resid_vs_als_lm_jun3.png"),
       p4, 
       units="in", 
       width=8, 
       height=6, 
       dpi=500)

# tomorrow -- use shift values to shift all LAI values 