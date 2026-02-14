# -----------------------------------------------------------------------------_
#  make plots of recovery curves for LAI, veg height, and fcover frst params 
#
# ------------------------------------------------------------------------------

# library 
library(terra)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(ggh4x)

setwd("C:/Users/blbouf/Sync/Paper2")


# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

veg_path <- file.path(".", "data", "forest_params_by_age", "sampled_veg_params_byLAI_GRP_2015_2021_feb10.csv")

# ------------------------------------------------------------------------------
# read data 
# ------------------------------------------------------------------------------

veg_data <- read.csv(veg_path) 

# age breaks  
age_breaks_df <- data.frame(
  class   = c("D_g0", "R1_g0", "R2_g0", "R3_g0", "R4_g0", "R5_g0", 
              "R6_g0", "R7_g0", "R8_g0", "R9_g0", "M_g0"),
  age_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  )

# quartiles for plotting
Qup <- 0.9
Qlow <- 0.1

# ------------------------------------------------------------------------------
# organize data 
# ------------------------------------------------------------------------------

# get med jul LAI for each age class and BEC zone -- for plotting later 
veg_data <- veg_data[veg_data$age < 100, ]

med_veg_byage  <- veg_data %>%
  group_by(age, lai_grp) %>%
  summarize(av_jul_lai = median(july_lai, na.rm = TRUE),
            med_ntems_fcover = median(ntems_fcover, na.rm = TRUE), 
            med_ntems_height = median(ntems_height, na.rm = TRUE),
            med_peak_fcover = median(peak_fcover_leaf, na.rm = TRUE),
            n = length(july_lai))

# factor id column 
med_veg_byage <- med_veg_byage %>%
  mutate(id = factor(lai_grp))

elev_zoned_classes <- list(
  g0 = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_g0", "R1_g0", "R2_g0", "R3_g0", "R4_g0", "R5_g0", 
                "R6_g0", "R7_g0", "R8_g0", "R9_g0", "M_g0")
  ),
  g1 = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_g1", "R1_g1", "R2_g1", "R3_g1", "R4_g1", "R5_g1", 
                "R6_g1", "R7_g1", "R8_g1", "R9_g1","M_g1")
  )
)

zoned_breaks <- bind_rows(elev_zoned_classes, .id = "lai_grp") %>%
  dplyr::select(lai_grp, min_age) %>%
  distinct()

plot_data_grouped <- veg_data %>%
  group_by(lai_grp, age) %>%
  summarise(LAIhigh = quantile(med_winter_lai, Qup, na.rm = TRUE),
            LAIlow = quantile(med_winter_lai, Qlow, na.rm = TRUE),
            mean_lai = mean(med_winter_lai, na.rm = TRUE),
            
            ntems_fcover_high = quantile(ntems_fcover, Qup, na.rm = TRUE),
            ntems_fcover_low = quantile(ntems_fcover, Qlow, na.rm = TRUE),
            mean_ntems_fcover = mean(ntems_fcover, na.rm = TRUE),
            
            ntems_height_high = quantile(ntems_height, Qup, na.rm = TRUE),
            ntems_height_low = quantile(ntems_height, Qlow, na.rm = TRUE),
            mean_ntems_height = mean(ntems_height, na.rm = TRUE),
            
            peak_fcover_high = quantile(peak_fcover_leaf, Qup, na.rm = TRUE),
            peak_fcover_low = quantile(peak_fcover_leaf, Qlow, na.rm = TRUE),
            mean_peak_fcover = mean(peak_fcover_leaf, na.rm = TRUE)
            ) %>%
  ungroup() 

# make the lai_group column match the class matrix 
plot_data_grouped$lai_grp[plot_data_grouped$lai_grp == 0] <- "g0"
plot_data_grouped$lai_grp[plot_data_grouped$lai_grp == 1] <- "g1"

# Function to assign age class based on zone-specific breakpoints
assign_class <- function(age, lai_grp) {
  class_df <- elev_zoned_classes[[lai_grp]]
  if (is.null(class_df)) return(NA)
  idx <- which(age >= class_df$min_age & age < class_df$max_age)
  if (length(idx) == 0) return(NA)
  class_df$class[idx[1]]
}

# Apply the function rowwise to assign age class
age_class_plot_data <- plot_data_grouped %>%
  rowwise() %>%
  mutate(age_class = assign_class(age, lai_grp)) %>%
  ungroup()

# -----------------------------------------------------------------------------
age_range <- range(plot_data_grouped$age, na.rm = TRUE)
age_breaks <- seq(from = floor(age_range[1] / 5) * 5, 
                  to = ceiling(age_range[2] / 5) * 5, 
                  by = 5)

# zoned breaks for plots
zoned_breaks <- zoned_breaks[zoned_breaks$min_age != 0, ]

# LEAF LAI PLOT 
leaf_lai_plot <- ggplot(data = age_class_plot_data, aes(x = age, y = mean_lai, fill = age_class)) + 
  geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
  geom_ribbon(aes(ymin = LAIlow, ymax = LAIhigh), alpha = 0.25, colour = NA) + 
  geom_point() + 
  facet_wrap(~lai_grp, ncol = 1) +
  theme_minimal(12) + 
  scale_x_continuous(breaks = age_breaks, expand = c(0.01, 0.01)) +
  labs(x = "Forest Age", y = "Mean Peak Season LAI", fill = "Recovery Class") + 
  theme(legend.position = 'none',
        panel.spacing =unit(1, "cm"))
 

# LEAF FCOVER PLOT 
leaf_fcover_plot <- ggplot(data = age_class_plot_data, aes(x = age, y = mean_peak_fcover, fill = age_class)) + 
  geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
  geom_ribbon(aes(ymin = peak_fcover_low, ymax = peak_fcover_high), alpha = 0.25, colour = NA) + 
  geom_point() + 
  facet_wrap(~lai_grp, ncol = 2) +
  theme_minimal(12) + 
  scale_x_continuous(breaks = age_breaks, expand = c(0.01, 0.01)) +
  labs(x = "Forest Age", y = "Mean Peak Season fcover", fill = "Recovery Class") + 
  theme(legend.position = 'none',
        panel.spacing =unit(1, "cm"))


# NTEMS FCOVER PLOT 
ntems_fcover_plot <- ggplot(data = age_class_plot_data, aes(x = age, y = mean_ntems_fcover, fill = age_class)) + 
  geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
  geom_ribbon(aes(ymin = ntems_fcover_low, ymax = ntems_fcover_high), alpha = 0.25, colour = NA) + 
  geom_point() + 
  facet_wrap(~lai_grp, ncol = 2) +
  theme_minimal(12) + 
  scale_x_continuous(breaks = age_breaks, expand = c(0.01, 0.01)) +
  labs(x = "Forest Age", y = "Mean NTEMS fcover", fill = "Recovery Class") + 
  theme(legend.position = 'none',
        panel.spacing =unit(1, "cm"))


# NTEMS HEIGHT PLOT 
ntems_height_plot <- ggplot(data = age_class_plot_data, aes(x = age, y = mean_ntems_height, fill = age_class)) + 
  geom_vline(data = zoned_breaks, aes(xintercept = min_age), linetype = "dashed", colour = "red") +
  geom_ribbon(aes(ymin = ntems_height_low, ymax = ntems_height_high), alpha = 0.25, colour = NA) + 
  geom_point() + 
  facet_wrap(~lai_grp, ncol = 2) +
  theme_minimal(12) + 
  scale_x_continuous(breaks = age_breaks, expand = c(0.01, 0.01)) +
  labs(x = "Forest Age", y = "Mean NTEMS height", fill = "Recovery Class") + 
  theme(legend.position = 'none',
        panel.spacing =unit(1, "cm"))



ggsave(av_plot,
       filename = file.path(figs_path, "recovery_curves_dec4.png"),
       units="in", 
       height= 12, width=15,
       dpi=600)

# calculate stats to go into paper here !

# subset data
all_dat <- jul_plot_data_simple 
idf_dat <- jul_plot_data_adv[jul_plot_data_adv$ZONE == "IDF", ]
ms_dat <- jul_plot_data_adv[jul_plot_data_adv$ZONE == "MS", ]
essf_dat <- jul_plot_data_adv[jul_plot_data_adv$ZONE == "ESSF", ]

# mature ages
# MS 44 
# ESSF 50
# IDF 33 
# Catchment scale 46

# print out recovery stats 
get_recovery_stats(all_dat, 46, "Catchment")
get_recovery_stats(essf_dat, 50, "ESSF")
get_recovery_stats(idf_dat, 33, "IDF")
get_recovery_stats(ms_dat, 44, "MS")



# function to get recovery stats 
get_recovery_stats <- function(df, mature_age, df_name){
  
  av_mature_LAI <-  mean(df$july_lai[df$age > mature_age], na.rm = TRUE)
  av_5_LAI <-  mean(df$july_lai[df$age == 5], na.rm = TRUE)
  av_10_LAI <-  mean(df$july_lai[df$age == 10], na.rm = TRUE)
  
  recovery_after_5 <- round(av_5_LAI/av_mature_LAI*100, 2)
  recovery_after_10 <- round(av_10_LAI/av_mature_LAI*100, 2)
  
  result <- data.frame(c(recovery_after_5, recovery_after_10, df_name))
  
  return(result)
}


get_recovery_stats(jul_plot_data_simple, 46, "Catchment_based")


