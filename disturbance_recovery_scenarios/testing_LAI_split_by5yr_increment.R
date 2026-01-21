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

# set wd 
setwd(file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", 
                "scripts", "LAI_recovery"))

# paths 
output_path <- file.path("..","..","data", "SOM_outputs", "recovery_results")
figs_path <- file.path("..","..","data", "SOM_outputs", "recovery_results", "figs")

result_path <- file.path(file.path("..", "..", "..", "raven-runs", "_Trapping_LAI", "data", "monthly_LAI"))

basic_path <- file.path(result_path, "median_seasonal_curve_basic_aug7.csv")

simple_path <- file.path(result_path, "median_seasonal_curve_simplerecovery_jul2.csv") # was jun5

adv_path <- file.path(result_path, "median_seasonal_curve_lai_recovery_advanced_jul2.csv") # was jun5
adv_path <- file.path(result_path, "median_seasonal_curve_lai_recovery_nov26.csv") 
adv_path <- file.path(result_path, "unsampled_median_seasonal_curve_lai_recovery_nov27.csv") 
# read data and add key per dataset 
# basic_lai <- read.csv(basic_path)
# basic_lai$id <- "Simple LAI"
# simple_lai <- read.csv(simple_path)
# simple_lai$id <- "Full Region Recovery-Based LAI"
# 
adv_lai <- read.csv(adv_path)
#adv_lai <- result
adv_lai$id <- "BEC Recovery-Based LAI"

# clean adv data 
adv_lai$id[grepl("essf", adv_lai$age_class)] <- paste0(adv_lai$id[grepl("essf", adv_lai$age_class)], " ESSF")
adv_lai$id[grepl("idf", adv_lai$age_class)] <- paste0(adv_lai$id[grepl("idf", adv_lai$age_class)], " IDF")
adv_lai$id[grepl("ms", adv_lai$age_class)] <- paste0(adv_lai$id[grepl("ms", adv_lai$age_class)], " MS")

# # join data 
# lai_df <- rbind(basic_lai, 
#                 simple_lai) %>%
#   rbind(., adv_lai)

lai_df <- adv_lai

lai_df <- filter(lai_df, !(age_class %in% c("ALPINE", "WETLAND", "SHRUB")))

# clean labels
classes <- lai_df$age_class
#classes <- sapply(strsplit(classes, "_"), "[", 1)
lai_df$age_class_group <- classes
# lai_df$age_class_group[lai_df$age_class_group == "old"] <- "Mature Forest"
# lai_df$age_class_group[lai_df$age_class_group == "mature"] <- "Mature Forest"
# lai_df$age_class_group[lai_df$age_class_group == "young"] <- "Self-thinning"
# lai_df$age_class_group[lai_df$age_class_group == "early"] <- "Stand Regeneration"
# lai_df$age_class_group[lai_df$age_class_group == "late"] <- "Rapid Canopy Development"
# lai_df$age_class_group[lai_df$age_class_group == "recovery"] <- "Rapid Canopy Development"
# lai_df$age_class_group[lai_df$age_class_group == "disturbed"] <- "Disturbed"

# clean group labels 
lai_df$plotting_groups <- lai_df$age_class_group
# lai_df$plotting_groups[lai_df$age_class == "mature_essf"] <- "Mature ESSF"
# lai_df$plotting_groups[lai_df$age_class == "mature_idf"] <- "Mature IDF"
# lai_df$plotting_groups[lai_df$age_class == "mature_ms"] <- "Mature MS"


lai_df_long <- lai_df %>%
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

lai_df_long$id <- factor(lai_df_long$id, levels = c(
  "Simple LAI",
  "Full Region Recovery-Based LAI",
  "BEC Recovery-Based LAI ESSF",
  "BEC Recovery-Based LAI IDF",
  "BEC Recovery-Based LAI MS"), 
  labels = c("Baseline",
             "Catchment \nScale",
             "Ecosystem \nScale ESSF",
             "Ecosystem \nScale IDF",
             "Ecosystem \nScale MS"
  ))


lai_df_long$TYPE <- "All"
lai_df_long$TYPE[grepl("ESSF", lai_df_long$plotting_groups)] <- "ESSF"
lai_df_long$TYPE[grepl("IDF", lai_df_long$plotting_groups)] <- "IDF"
lai_df_long$TYPE[grepl("MS", lai_df_long$plotting_groups)] <- "MS"
lai_df_long$TYPE[grepl("ESSF", lai_df_long$id)] <- "ESSF"
lai_df_long$TYPE[grepl("IDF", lai_df_long$id)] <- "IDF"
lai_df_long$TYPE[grepl("MS", lai_df_long$id)] <- "MS"


mean_lai_bec_facet <- ggplot(lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                              color = age_class)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  geom_line(data = lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                    color = age_class), linewidth = 1) + 
  facet_wrap(~ id, ncol = 3) +
  theme_minimal(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12),
        legend.position = "bottom") +
  guides(color = guide_legend(order = 1)) 

# ------------------------------------------------------------------------------
# SILVILASER PLOT 
# ------------------------------------------------------------------------------
lai_df_long <- lai_df %>%
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

lai_df_svlsr <- lai_df_long[!lai_df_long$id %in% c("BEC Recovery-Based LAI ESSF",
                                                   "BEC Recovery-Based LAI IDF",
                                                   "BEC Recovery-Based LAI MS"), ]

lai_df_svlsr$id <- factor(lai_df_svlsr$id, levels = c(
  "Simple LAI",
  "Full Region Recovery-Based LAI"),
  labels = c("Simple LAI",
             "Catchment Scale \nLAI"
  ))

mean_lai_svlsr <- ggplot(lai_df_svlsr, aes(x = month_clean, y = lai, group = plotting_groups, 
                                           color = age_class_group)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  geom_line(data = lai_df_svlsr, aes(x = month_clean, y = lai, group = plotting_groups, 
                                     color = age_class_group), linewidth = 1) +  
  facet_wrap(~ id, ncol = 2, labeller = labeller(id = 
                                                   c("Simple LAI" = "Baseline LAI", 
                                                     "Catchment Scale \nLAI" = "Data-driven LAI"))) +
  labs(x = "Month", 
       y = "Median LAI", 
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
  ) 

ggsave(mean_lai_svlsr,
       filename = file.path(figs_path, "med_lai_seasonan_trends_svlsr_sept26.png"),
       units="in", 
       height= 5, width=10,
       dpi=600)
