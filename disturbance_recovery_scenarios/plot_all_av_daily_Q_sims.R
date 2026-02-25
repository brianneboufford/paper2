# ------------------------------------------------------------------------------
# plot average daily streamflow for all simulations 
#
# 
# script based on 01_check_out_averag_daily_Q.R
# Feb 23rd, 2026
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
library(gridExtra)

setwd("C:/Users/blbouf/Sync/Paper2")
# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Feb20")
hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios_Feb19", "Runs")

hydrograph_list_all <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                                  recursive = TRUE, 
                                  full.name = TRUE)

# ------------------------------------------------------------------------------
# read and join data 
# ------------------------------------------------------------------------------

hyd_data <- lapply(hydrograph_list_all, prep_sim_data) %>% 
  do.call(rbind, .)

hyd_data <- hyd_data[!grepl("_40_", hyd_data$sim), ]
hyd_data <- hyd_data[!grepl("Observed", hyd_data$Type), ]

af_data <- hyd_data[hyd_data$sim == "all_forest", ]
sim_data <- hyd_data[!hyd_data$sim == "all_forest", ]

Qup = 0.9
Qlow = 0.1

sim_data <- sim_data %>% 
  mutate("elev_zone" = str_split(sim, pattern = "_", simplify = TRUE)[, 2], 
         "percent_dist" = str_split(sim, pattern = "_", simplify = TRUE)[, 3],
         "yrs_recov" = str_split(sim, pattern = "_", simplify = TRUE)[, 4])

sim_data$yrs_recov <- str_remove(sim_data$yrs_recov, "yrs")
sim_data <- sim_data %>%
  subset(yday(date) %in% 105:166)

sim_data_grouped <- sim_data %>%
  group_by(sim, date = yday(date)) %>%
  summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
            Qlow = quantile(Value, Qlow, na.rm = TRUE),
            Value = mean(Value, na.rm = TRUE), 
            elev_zone = first(elev_zone),
            percent_dist = first(percent_dist),
            yrs_recov = first(yrs_recov)) %>%
  ungroup() 

af_data_grouped <- af_data %>%
  group_by(sim, date = yday(date)) %>%
  summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
            Qlow = quantile(Value, Qlow, na.rm = TRUE),
            Value = mean(Value, na.rm = TRUE)) %>%
  ungroup() 

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------

av_daily_facet_plot <- ggplot(data = sim_data_grouped, aes(x = date, y = Value, colour = percent_dist)) +
 # geom_ribbon(aes(ymin = Qlow, ymax = Qhigh), alpha = 0.25, colour = NA) +
  geom_line() +
  facet_grid(elev_zone ~ yrs_recov) + 
  # scale_x_date(date_labels = '%b') +
 # scale_colour_manual('', values = c('turquoise4', 'chocolate2')) +
#  scale_fill_manual('', values = c( 'turquoise4', 'chocolate2')) +
  labs(y = "aveage daily Q") +
  theme_bw(14) +
  theme(legend.position = 'bottom', 
        axis.title.x = element_blank())

av_daily_facet_plot

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

prep_sim_data <- function(df_path){
  
  df <- read.csv(df_path)
  sim_name <- dirname(df_path) %>% basename() %>% str_remove("Trapping_")
  
  df$date <- as.Date(df$date)
  df <- df[!(is.na(df$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)), ]
  
  names(df) <- c("time", "date", "hour", "precip", "Simulated", "Observed")
  df <- df %>% 
    na.omit() %>%
    pivot_longer(
      cols = c("Observed", "Simulated"),  # selects the last two columns
      names_to = "Type",
      values_to = "Value"
    )
  
  df$sim <- sim_name
  
  return(df)
} 
