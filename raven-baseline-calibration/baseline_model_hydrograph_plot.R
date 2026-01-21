# ------------------------------------------------------------------------------
# baseline model hydrograph plot
#
# 
# 
# December 4, 2025
# ------------------------------------------------------------------------------

# packages 
library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(lubridate)
library(tidyr)
library(cowplot)
library(HyMETT)
library(Metrics)
library(hydrostats)
library(cowplot)
library(ggpattern)
library(fasstr)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

hru_run_path <- file.path(".", "raven-runs", "Baseline2", "Runs", "Trapping_HRU_baseline_Dec10")
fig_path <- file.path(".", "data", "figs", "baseline_model_calibration")

# list of hydrographs
hydrograph_list <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                              recursive = TRUE, 
                              full.name = TRUE)

hydro_data <- read.csv(hydrograph_list)
hydro_data$date <- as.Date(hydro_data$date)
hydro_data <- na.omit(hydro_data)
hydro_data <- hydro_data[hydro_data$time > 365, ]

hydroplot_baseline <- ggplot() + 
  geom_line(data = hydro_data, aes(x=date, y=TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., color = "Observed")) + 
  geom_line(data = hydro_data, aes(x=date, y=TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s., color = "Modelled")) + 
  theme_minimal(16) + 
  scale_color_manual(name = "", values = c("Observed" = "#4575b4", "Modelled" = "#d73027")) +
  labs(title="", x="Date", y=paste0("Streamflow (m\u00B3/s)")) + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), 
        legend.position = c(0.10, 0.90),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA),  # white background, no border
        legend.key = element_rect(fill = "white", colour = NA))

ggsave(hydroplot_baseline,
       filename = file.path(fig_path, "daily_hydrograph_baseline_paper2_dec10.png"),
       units = "in",
       width = 10,
       height = 4)

# average discharge calculation 
av_discharge <- mean(hydro_data$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)
