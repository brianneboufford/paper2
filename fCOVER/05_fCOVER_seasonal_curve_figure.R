# ------------------------------------------------------------------------------
# plots median monthly fCOVER for the ecosystem scale classes
# 
# adapted from 005_LAI_seasonal_curves_figure
# Nov 19, 2025
#-------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

output_prod_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "average_monthly_fCOVER")
figs_path <- file.path(".", "data", "figs", "fCOVER_analysis")
av_data_path <- file.path(output_prod_path, "median_seasonal_curve_fcover_nov19.csv")

av_fcover <- read.csv(av_data_path)

av_fcover$id[grepl("essf", av_fcover$age_class)] <- " ESSF"
av_fcover$id[grepl("idf", av_fcover$age_class)] <- " IDF"
av_fcover$id[grepl("ms", av_fcover$age_class)] <- " MS"

fcov_df <- filter(av_fcover, !(age_class %in% c("ALPINE", "WET_LAND", "SHRUB")))

# clean labels
classes <- fcov_df$age_class
classes <- sapply(strsplit(classes, "_"), "[", 1)
fcov_df$age_class_group <- classes
fcov_df$age_class_group[fcov_df$age_class_group == "old"] <- "Mature Forest"
fcov_df$age_class_group[fcov_df$age_class_group == "mature"] <- "Mature Forest"
fcov_df$age_class_group[fcov_df$age_class_group == "young"] <- "Self-thinning"
fcov_df$age_class_group[fcov_df$age_class_group == "early"] <- "Stand Regeneration"
fcov_df$age_class_group[fcov_df$age_class_group == "late"] <- "Rapid Canopy Development"
fcov_df$age_class_group[fcov_df$age_class_group == "recovery"] <- "Rapid Canopy Development"
fcov_df$age_class_group[fcov_df$age_class_group == "disturbed"] <- "Disturbed"

# clean group labels 
fcov_df$plotting_groups <- fcov_df$age_class_group
fcov_df$plotting_groups[fcov_df$age_class == "mature_essf"] <- "Mature ESSF"
fcov_df$plotting_groups[fcov_df$age_class == "mature_idf"] <- "Mature IDF"
fcov_df$plotting_groups[fcov_df$age_class == "mature_ms"] <- "Mature MS"

fcov_df_long <- fcov_df %>%
  pivot_longer(
    cols = matches(paste0("med_", "apr|may|june|july|aug|sept|oct|winter", "$")),
    names_to = "month",
    values_to = "fCOVER"
  )

# Clean the month names
fcov_df_long <- fcov_df_long %>%
  mutate(month_clean = stringr::str_replace(month, "med_", "")) %>%
  filter(month_clean %in% c("apr", "may", "june", "july", "aug", "sept", "oct", "winter"))

# Ensure correct order of months
ordered_months <- c("apr", "may", "june", "july", "aug", "sept", "oct", "winter")
fcov_df_long$month_clean <- factor(fcov_df_long$month_clean, levels = ordered_months, 
                                  labels = c("04", "05", "06", "07", "08", "09", "10", "winter"))

# ------------------------------------------------------------------------------
# PLOT 
# ------------------------------------------------------------------------------
mean_fcover_plot <- ggplot(fcov_df_long, aes(x = month_clean, y = fCOVER, group = plotting_groups, 
                                                          color = age_class_group)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  facet_wrap(~ id, ncol = 5) +
  labs(x = "Month", 
       y = "Median fCOVER", 
       color = "Recovery Class", 
       fill = "Recovery Class") + 
  theme_minimal(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12),
        legend.position = "right") +
  scale_color_manual(
    values = c("Disturbed" = "#8c510a", "Stand Regeneration" = "#4575b4", 
               "Rapid Canopy Development" = "#35978f", "Self-thinning" = "#5aae61", 
               "Mature Forest" = "#00441b"),
    breaks = c("Disturbed", "Stand Regeneration",
               "Rapid Canopy Development", "Self-thinning", "Mature Forest")
  ) +
  guides(
    color    = guide_legend(order = 2)    
  )
