# ------------------------------------------------------------------------------
# adapted from Silvilaser boxplot script to generate new box plots for paper 
# october 7th, 2025
# updated october 27th to use mean value instead of residual for plot 
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

# ------------------------------------------------------------------------------
# read data 
# ------------------------------------------------------------------------------

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
lai_val_summary_long <- lai_val_summary %>% 
  dplyr::select(-c("median_lidar", "median_hls")) %>% 
  pivot_longer(data = ., 
               cols = c("mean_hls", "mean_lidar"),
               names_to = "data", 
               values_to = "lai_stat") 

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

# pivot to format that is easier to plot as bar plot 
shifted_val_summary_long <- lai_shift_results %>% 
  pivot_longer(data = ., 
               cols = c("mean_adj", "mean_lidar"),
               names_to = "data", 
               values_to = "lai_stat")

clean_box_data1 <- lai_recon_sample %>% dplyr::select(c("category", "diff_lai"))
clean_box_data1$type <- "LEAF"

clean_box_data2 <- shifted_lai %>% dplyr::select(c("category", "diff_shift")) 
names(clean_box_data2) <- c("category", "diff_lai")
clean_box_data2$type <- "Adjusted LEAF"

clean_box_data <- rbind(clean_box_data1, clean_box_data2)
clean_box_data$type <- factor(clean_box_data$type, levels = c("LEAF", "Adjusted LEAF"))
clean_box_data$category <- factor(clean_box_data$category, 
                                  levels = c("Shrubs", 
                                             "Wetland",
                                             "Wetland-Treed",
                                             "Herbs",
                                             "Broadleaf",
                                             "Mixed Wood",
                                             "Coniferous_0",
                                             "Coniferous_1",
                                             "Coniferous_2",
                                             "Coniferous_3"))

x_labels <- c(
  "Shrubs" = "Shrub",
  "Wetland" = "WL",
  "Wetland-Treed" = "WL-Treed",
  "Herbs" = "Herb",
  "Broadleaf" = "BL", 
  "Mixed Wood" = "MW",
  "Coniferous_0" = "Conif \n<20%",
  "Coniferous_1" = "Conif \n20???50%",
  "Coniferous_2" = "Conif \n50???80%",
  "Coniferous_3" = "Conif \n>80%"
)


# ----------------------------------------------------------------------------
# reorganize data for mean value plot 
# ----------------------------------------------------------------------------
hls_data <- shifted_lai %>% dplyr::select(c("category", "hls_lai"))
names(hls_data) <- c("category", "lai")
hls_data$type <- "LEAF"

shifted_hls_data <- shifted_lai %>% dplyr::select(c("category", "hls_lai_new"))
names(shifted_hls_data) <- c("category", "lai")
shifted_hls_data$type <- "Adjusted LEAF"

lidar_data <- shifted_lai %>% dplyr::select(c("category", "lidar_lai"))
names(lidar_data) <- c("category", "lai")
lidar_data$type <- "ALS"

clean_shift_data <- rbind(hls_data, shifted_hls_data) %>% 
  rbind(., lidar_data)
clean_shift_data$lai[clean_shift_data$lai < 0] <- 0
clean_shift_data$type <- factor(clean_shift_data$type, levels = c("LEAF", "Adjusted LEAF", "ALS"))
clean_shift_data$category <- factor(clean_shift_data$category, 
                                    levels = c("Shrubs", 
                                               "Wetland",
                                               "Wetland-Treed",
                                               "Herbs",
                                               "Broadleaf",
                                               "Mixed Wood",
                                               "Coniferous_0",
                                               "Coniferous_1",
                                               "Coniferous_2",
                                               "Coniferous_3"))




# box plot of difference between LAI values based on NTEMS landcover class
boxplot <- ggplot(data = clean_box_data, aes(x = category, y = diff_lai, fill = type)) + 
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05, width=0.6, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("LEAF" = "grey83", "Adjusted LEAF" = "grey26")) + 
  xlab("Landcover class") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab(expression(LAI[ALS]*" - "*LAI[LEAF])) +
  xlab("Landcover Class") + 
  theme_classic(16) +
  scale_x_discrete(labels = x_labels) + 
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "bottom",
        legend.title = element_blank()) + 
  stat_summary(fun = mean, geom = "point", shape = 8, size = 2, color = "black", position = position_dodge(width = 0.6),
               show.legend = FALSE) + 
  scale_y_continuous(limits=c(-2, 4), breaks=c(seq(-2, 4, 2))) 

boxplot_lai <- ggplot(data = clean_shift_data, aes(x = category, y = lai, fill = type)) + 
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05, width=0.6, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("LEAF" = "grey26", "Adjusted LEAF" = "grey63", "ALS" = "white")) + 
  xlab("Landcover class") +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("LAI") +
  xlab("Landcover Class") + 
  theme_classic(16) +
  scale_x_discrete(labels = x_labels) + 
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "bottom",
        legend.title = element_blank()) + 
  stat_summary(fun = mean, geom = "point", shape = 8, size = 2, color = "black", position = position_dodge(width = 0.6),
               show.legend = FALSE) 
#scale_y_continuous(limits=c(-2, 4), breaks=c(seq(-2, 4, 2))) 


# grid plot with frou panels
ggsave(boxplot, 
       filename = file.path(output_path, "results", "LEAF_residual_boxplot_march22.png"),
       units="in", 
       width=10, 
       height=5, 
       dpi=600)

ggsave(boxplot_lai, 
       filename = file.path(output_path, "results", "LEAF_boxplot_march22.png"),
       units="in", 
       width=10, 
       height=5, 
       dpi=900)
