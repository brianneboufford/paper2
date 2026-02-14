# -----------------------------------------------------------------------------
# plot median monthyl LAI by group  
# plot median fcover and height (NTEMS + ALS) by group 
#
# February 10th, 2026
# -----------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

setwd("C:/Users/blbouf/Sync/Paper2")

for_params_path <- file.path(".", "data", "median_forest_params_LAI_by_age")

# NTEMS forest cover and height 
h_fc_ntems_path <- file.path(for_params_path, "median_fparams_ELEV_recovery_feb10.csv") # was jun5

# LEAF LAI 
lai_path <- file.path(for_params_path, "median_lai_ELEV_recovery_feb10.csv")

# ALS forest cover and height
h_fc_als_path <- file.path(for_params_path, "median_fparams_ALS_ELEV_recovery_feb10.csv")

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

lai <- read.csv(lai_path)
lai$data <- "leaf"

h_fc_ntems <- read.csv(h_fc_ntems_path)
h_fc_ntems$data <- "ntems"
names(h_fc_ntems) <- c("age_class", 
                       "med_h", 
                       "med_fc",
                       "data")

h_fc_als <- read.csv(h_fc_als_path)
h_fc_als$data <- "als"
names(h_fc_als) <- c("age_class", 
                       "med_h", 
                       "med_fc",
                       "data")

fpar_df <- rbind(h_fc_ntems, 
                 h_fc_als) 

# clean labels
classes <- fpar_df$age_class

fpar_df$lai_grp <- NA
fpar_df$lai_grp[grepl("g0", fpar_df$age_class)] <- "g0"
fpar_df$lai_grp[grepl("g1", fpar_df$age_class)] <- "g1"

# clean group labels 
fpar_df$plot_num <- fpar_df$age_class
fpar_df$plot_num[grepl("M", fpar_df$plot_num)] <- 11
fpar_df$plot_num[grepl("D", fpar_df$plot_num)] <- 1
fpar_df$plot_num[grepl("R1", fpar_df$plot_num)] <- 2
fpar_df$plot_num[grepl("R2", fpar_df$plot_num)] <- 3
fpar_df$plot_num[grepl("R3", fpar_df$plot_num)] <- 4
fpar_df$plot_num[grepl("R4", fpar_df$plot_num)] <- 5
fpar_df$plot_num[grepl("R5", fpar_df$plot_num)] <- 6
fpar_df$plot_num[grepl("R6", fpar_df$plot_num)] <- 7
fpar_df$plot_num[grepl("R7", fpar_df$plot_num)] <- 8
fpar_df$plot_num[grepl("R8", fpar_df$plot_num)] <- 9
fpar_df$plot_num[grepl("R9", fpar_df$plot_num)] <- 10
fpar_df$plot_num <- as.numeric(fpar_df$plot_num)

# plot median height
median_height_plot <- ggplot(fpar_df, aes(x = plot_num, y = med_h, colour = data)) +
  geom_point() +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  facet_wrap(~lai_grp, ncol = 2) +
  labs(x = "recovery class", 
       y = "Median height") +   
  theme_minimal(12) +
  theme(legend.position = "right")

# plot median forest cover 
median_fcover_plot <- ggplot(fpar_df, aes(x = plot_num, y = med_fc, colour = data)) +
  geom_point() +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  facet_wrap(~lai_grp, ncol = 2) +
  labs(x = "recovery class", 
       y = "Median fcover") +   
  theme_minimal(12) +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# LAI 
# ------------------------------------------------------------------------------

# clean labels
classes <- lai$age_class

lai$lai_grp <- NA
lai$lai_grp[grepl("g0", lai$age_class)] <- "g0"
lai$lai_grp[grepl("g1", lai$age_class)] <- "g1"

lai_df_long <- lai %>%
  pivot_longer(
    cols = ends_with("lai"),
    names_to = "month",
    values_to = "lai"
  )

# Clean the month names
lai_df_long <- lai_df_long %>%
  mutate(month_clean = stringr::str_replace(month, "_lai", "")) %>%
  mutate(month_clean = stringr::str_replace(month_clean, "med_", "")) %>%
  filter(month_clean %in% c("apr", "may", "june", "july", "aug", "sept", "oct", "winter"))

# Ensure correct order of months
ordered_months <- c("apr", "may", "june", "july", "aug", "sept", "oct", "winter")
lai_df_long$month_clean <- factor(lai_df_long$month_clean, levels = ordered_months, 
                                  labels = c("04", "05", "06", "07", "08", "09", "10", "winter"))


mean_lai_from_recovery_classes <- ggplot(lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                                          color = age_class)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  geom_line(data = lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                    color = age_class), linewidth = 1) +  
  facet_wrap(~lai_grp, ncol = 2) +
  labs(x = "Month", 
       y = "Median LAI", 
       color = "Recovery Class", 
       fill = "Recovery Class") +   
  theme_minimal(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12),
        legend.position = "right")
