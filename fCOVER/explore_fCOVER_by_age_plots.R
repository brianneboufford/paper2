# -----------------------------------------------------------------------------_
# script to plot and analyze fCOVER trends by age  
#
# adapted from 01_growth_curve_fitting_LAI_recovery 
# 
# november 19, 2025
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
sampled_data_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "fCOVER_age",
                               "sampled_fCOVER_all_data_2014_2021_nov19.csv")

# ------------------------------------------------------------------------------
# read data 
# ------------------------------------------------------------------------------
sampled_data <- read.csv(sampled_data_path)

# zon ekey 
zone_key <- read.csv(file.path(".", "data", "src", "BEC", "zone_key.csv"))
zone_key <- zone_key[c(1,2,5), ]

# get mean jul LAI for each age class and BEC zone -- for plotting later 
sampled_data <- sampled_data[sampled_data$age < 100, ]
sampled_data <- merge(sampled_data, zone_key, by="id")

av_jul_age <- sampled_data %>%
  group_by(age, id) %>%
  summarize(av_jul = median(july_fCOVER, na.rm = TRUE),
            n = length(july_fCOVER))

# factor id column 
av_jul_age <- av_jul_age %>%
  mutate(id = factor(id))

# ------------------------------------------------------------------------------
# exploratory plots of july/oct forest cover and median july/oct forest cover by age 
# ------------------------------------------------------------------------------
jul_plot_data <- sampled_data %>% 
  select(c("age", "ZONE", "july_fCOVER"))

oct_plot_data <- sampled_data %>% 
  select(c("age", "ZONE", "oct_fCOVER"))

Qup <- 0.9
Qlow <- 0.1

july_plot_data_grouped <- jul_plot_data %>%
  group_by(ZONE, age) %>%
  summarise(LAIhigh = quantile(july_fCOVER, Qup, na.rm = TRUE),
            LAIlow = quantile(july_fCOVER, Qlow, na.rm = TRUE),
            med_fCOVER = mean(july_fCOVER, na.rm = TRUE)) %>%
  ungroup() 

oct_plot_data_grouped <- oct_plot_data %>%
  group_by(ZONE, age) %>%
  summarise(LAIhigh = quantile(oct_fCOVER, Qup, na.rm = TRUE),
            LAIlow = quantile(oct_fCOVER, Qlow, na.rm = TRUE),
            med_fCOVER = mean(oct_fCOVER, na.rm = TRUE)) %>%
  ungroup() 


oct_av_plot <- ggplot(data = oct_plot_data_grouped, aes(x = age, y = med_fCOVER)) +
  geom_ribbon(aes(ymin = LAIlow, ymax = LAIhigh), alpha = 0.25) +
  geom_point() +
  facet_wrap(~ZONE, ncol = 1) +
  # geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
  # geom_line(data = rich_df_long, aes(x = age, y = PeakLAI), color = "black", inherit.aes = FALSE, linewidth = 1) +
  theme_minimal() +
  # scale_x_continuous(breaks = age_breaks) +
  labs(x = "Age", y = "Median october fCOVER")# + 
# theme(legend.position = 'bottom') +
#scale_fill_manual(
#  values = c("Disturbed" = "#8c510a", "Early Recovery" = "#4575b4", 
#             "Late Recovery" = "#35978f", "Young Forest" = "#5aae61", "Mature Forest" = "#00441b"),
#  breaks = c("Disturbed", "Early Recovery",
#             "Late Recovery", "Young Forest", "Mature Forest")  # ensures legend order
)

jul_av_plot <- ggplot(data = july_plot_data_grouped, aes(x = age, y = med_fCOVER)) +
  geom_ribbon(aes(ymin = LAIlow, ymax = LAIhigh), alpha = 0.25) +
  geom_point() +
  facet_wrap(~ZONE, ncol = 1) +
 # geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
 # geom_line(data = rich_df_long, aes(x = age, y = PeakLAI), color = "black", inherit.aes = FALSE, linewidth = 1) +
  theme_minimal() +
 # scale_x_continuous(breaks = age_breaks) +
  labs(x = "Age", y = "Median Peak Season fCOVER")# + 
 # theme(legend.position = 'bottom') +
  #scale_fill_manual(
  #  values = c("Disturbed" = "#8c510a", "Early Recovery" = "#4575b4", 
  #             "Late Recovery" = "#35978f", "Young Forest" = "#5aae61", "Mature Forest" = "#00441b"),
  #  breaks = c("Disturbed", "Early Recovery",
  #             "Late Recovery", "Young Forest", "Mature Forest")  # ensures legend order
  )

# ------------------------------------------------------------------------------
# plot monthly fCOVER 
# ------------------------------------------------------------------------------
df_long <- sampled_data %>%
  # pivot all columns ending in "fCOVER"
  pivot_longer(
    ends_with("fCOVER"),
    names_to = "month",
    values_to = "fCOVER"
  ) %>%
  filter(!is.na(fCOVER)) %>%
  mutate(
    # extract month
    month = str_remove(month, "_fCOVER"),
    month = factor(
      month,
      levels = c("winter", "april", "may", "june", "july",
                 "aug", "sept", "oct")
    ),
    
    # ---- BIN AGES INTO 5-YEAR CLASSES ----
    age_bin = cut(
      age,
      breaks = seq(0, max(age, na.rm = TRUE) + 5, by = 5),
      include.lowest = TRUE,
      right = FALSE
    )
  )

# summarize
df_summary <- df_long %>%
  group_by(frst_cl, age_bin, month) %>%
  summarize(mean_fCOVER = mean(fCOVER), .groups = "drop")

# plot
ggplot(df_summary, aes(x = month, y = mean_fCOVER,
                       color = age_bin, group = age_bin)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~ frst_cl) +
  labs(
    x = "Month",
    y = "Mean fCOVER",
    color = "Age class (years)",
    title = "Average Monthly fCOVER by Forest Class and Age Group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

# function to get mode - built in R function doesn't compute the statistical mode 
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
