# ------------------------------------------------------------------------------
# script to make boxplots for fCOVER before and after lidar shift
#
# adapted from 003_generate_box_plots.R
# adapted from Silvilaser boxplot script to generate new box plots for paper 

# version created November 18, 2025
# (adapted from october 7th, 2025
# updated october 27th to use mean value instead of residual for plot)
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

# ------------------------------------------------------------------------------
# read data 
# ------------------------------------------------------------------------------

cc_recon_sample <- read.csv(file.path(output_path, "fCOVER_shift_train_sample.csv"))
cc_recon_hold <- read.csv(file.path(output_path, "fCOVER_shift_hold.csv"))

# ------------------------------------------------------------------------------
# prep data for bar plot 
# ------------------------------------------------------------------------------

# summarize mean of HLS LAI and lidar LAI
cc_val_summary <- cc_recon_sample %>%
  group_by(category) %>%
  summarize(median_hls = median(leaf_fcover_mean, na.rm = TRUE),
            median_lidar = median(lidar_cc, na.rm = TRUE),
            mean_hls = mean(leaf_fcover_mean, na.rm = TRUE),
            mean_lidar = mean(lidar_cc, na.rm = TRUE)) %>%
  as.data.frame()

# pivot to format that is easier to plot as bar plot 
cc_val_summary_long <- cc_val_summary %>% 
  select(-c("median_lidar", "median_hls")) %>% 
  pivot_longer(data = ., 
               cols = c("mean_hls", "mean_lidar"),
               names_to = "data", 
               values_to = "cc_stat") 

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

# caluclate percent area 
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
shifted_cc$hls_cc_new <- shifted_cc$leaf_fcover_mean + shifted_cc$mean_diff_cc

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

# ------------------------------------------------------------------------------
# NOTES: 
# adding a scale factor doesn't change the R2 and the PBIAS is worse so it is better
# to just add the mean difference to the LAI values 
# ------------------------------------------------------------------------------

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

# take out EXPOSED LAND --> made it worse 

library(knitr)
library(kableExtra)

overall_results$rmse_diff <- overall_results$rmse1 - overall_results$rmse0
overall_results$pbias_diff <- overall_results$pbias1 - overall_results$pbias0
kable(overall_results, format = "html") %>%
  kable_styling()

kable(cc_summary_export, format = "html") %>%
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
shifted_val_summary_long <- cc_shift_results %>% 
  pivot_longer(data = ., 
               cols = c("mean_adj", "mean_lidar"),
               names_to = "data", 
               values_to = "lai_stat")

clean_box_data1 <- cc_recon_sample %>% select(c("category", "diff_cc"))
clean_box_data1$type <- "LEAF"

clean_box_data2 <- shifted_cc %>% select(c("category", "diff_shift")) 
names(clean_box_data2) <- c("category", "diff_cc")
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
  "Coniferous_1" = "Conif \n20–50%",
  "Coniferous_2" = "Conif \n50–80%",
  "Coniferous_3" = "Conif \n>80%"
)


# ----------------------------------------------------------------------------
# reorganize data for mean value plot 
# ----------------------------------------------------------------------------
hls_data <- shifted_cc %>% select(c("category", "leaf_fcover_mean"))
names(hls_data) <- c("category", "cc")
hls_data$type <- "LEAF"

shifted_hls_data <- shifted_cc %>% select(c("category", "hls_cc_new"))
names(shifted_hls_data) <- c("category", "cc")
shifted_hls_data$type <- "Adjusted LEAF"

lidar_data <- shifted_cc %>% select(c("category", "lidar_cc"))
names(lidar_data) <- c("category", "cc")
lidar_data$type <- "ALS"

clean_shift_data <- rbind(hls_data, shifted_hls_data) %>% 
  rbind(., lidar_data)
clean_shift_data$cc[clean_shift_data$cc < 0] <- 0
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
boxplot <- ggplot(data = clean_box_data, aes(x = category, y = diff_cc, fill = type)) + 
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
               show.legend = FALSE)# + 
 # scale_y_continuous(limits=c(-2, 4), breaks=c(seq(-2, 4, 2))) 

boxplot_cc <- ggplot(data = clean_shift_data, aes(x = category, y = cc, fill = type)) + 
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

figs_path <- file.path(".", "data", "figs", "fCOVER_shift")
# grid plot with frou panels
ggsave(boxplot, 
       filename = file.path(figs_path, "fCOVER_residual_boxplot_nov18.png"),
       units="in", 
       width=10, 
       height=5, 
       dpi=600)

ggsave(boxplot_cc, 
       filename = file.path(figs_path, "fCOVER_boxplot_nov18.png"),
       units="in", 
       width=10, 
       height=5, 
       dpi=900)
