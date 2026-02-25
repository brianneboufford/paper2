# script based on 07_lai_baseline_report_metrics_jun3.R
# November 26th, 2025
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

hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios", "Runs")
fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Feb20")
hru_run_path <- file.path(".", "raven-runs", "Baseline3", "Runs", "Trapping_HRU_baseline_Feb20") # was dec 10
hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios_Feb19", "Runs")
# list of hydrographs
hydrograph_list_all <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                                  recursive = TRUE, 
                                  full.name = TRUE)

hydrograph_list_af <- hydrograph_list_all[1] %>% read.csv()

hydrograph_0yrs_list <- hydrograph_list_all[grepl(hydrograph_list_all, pattern = "_0yrs")]
hydrograph_0yrs_list <- hydrograph_0yrs_list[!grepl(hydrograph_0yrs_list, pattern = "_40_")]

hydroaf <- prep_model_data_af(hydrograph_list_af, "af") %>% rbind()

hydro_data <- lapply(hydrograph_0yrs_list, prep_model_data_sims) %>% 
  do.call(rbind, .)

# HERE facet plot by Site for hydro data and to each plot add the AF hydro data !!

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

prep_model_data_af <- function(df, data_name){
  names(df) <- c("time", "date", "hour", "precip", "Simulated", "Observed")
  df <- df %>% 
    na.omit() %>%
    pivot_longer(
      cols = c("Observed", "Simulated"),  # selects the last two columns
      names_to = "Type",
      values_to = "Value"
    )
  df$Site <- data_name
  return(df)
} 

prep_model_data_sims <- function(df_path){
  
  sim_name <- dirname(df_path) %>% basename() %>% str_remove("Trappingp50_") %>%
    str_remove("_0yrs") %>% str_replace("_", " ") %>% paste0(., "%")
  
  df <- read.csv(df_path)
  
  
  names(df) <- c("time", "date", "hour", "precip", "Simulated", "Observed")
  df <- df %>% 
    na.omit() %>%
    pivot_longer(
      cols = c("Observed", "Simulated"),  # selects the last two columns
      names_to = "Type",
      values_to = "Value"
    )
  df$Site <- sim_name
  return(df)
} 

plot_average_daily <- function(model_data, af_data, data_type = 'Streamflow', Qup = 0.9, Qlow = 0.1, wrap_gen = Inf, drop_nas = T){
  
  model_data_grouped <- model_data %>%
    #{if(drop_nas) spread(., Type, Value) %>% drop_na() %>% gather(., Type, Value, -Site, -date) else . } %>%
    group_by(Site, date = yday(date), Type) %>%
    #mutate(date = yday(date)) %>%
    group_by(Site, date, Type) %>%
    summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
              Qlow = quantile(Value, Qlow, na.rm = TRUE),
              Value = mean(Value, na.rm = TRUE)) %>%
    ungroup() 
  
  af_data_grouped <- af_data %>%
    #{if(drop_nas) spread(., Type, Value) %>% drop_na() %>% gather(., Type, Value, -Site, -date) else . } %>%
    group_by(Site, date = yday(date), Type) %>%
    #mutate(date = yday(date)) %>%
    group_by(Site, date, Type) %>%
    summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
              Qlow = quantile(Value, Qlow, na.rm = TRUE),
              Value = mean(Value, na.rm = TRUE)) %>%
    ungroup() 
  
  af_data_grouped <- af_data_grouped[af_data_grouped$Type == "Simulated", ]
  model_data_grouped <- model_data_grouped[model_data_grouped$Type == "Simulated", ]
  
  all_data_grouped <- rbind(af_data_grouped, model_data_grouped)
  
  all_data_grouped$Site[all_data_grouped$Site == "af"] <- "All Forested"
  
  av_plot <- ggplot(data = all_data_grouped, aes(x = as.Date(strptime(date, format = '%j')), y = Value,
                                                 colour = Site, fill = Site)) +
    geom_ribbon(aes(ymin = Qlow, ymax = Qhigh), alpha = 0.25, colour = NA) +
    geom_line() +
    # facet_wrap(~forcats::fct_reorder(Site, Value, .fun = function(x){mean(x, na.rm = TRUE)}),
    #             scales = 'free_y', labeller = label_wrap_gen(width = wrap_gen)) +
    scale_x_date(date_labels = '%b') +
    scale_colour_manual('', values = c('turquoise4', 'chocolate2')) +
    scale_fill_manual('', values = c( 'turquoise4', 'chocolate2')) +
    labs(y = y_lab) +
    theme_bw(14) +
    theme(legend.position = 'bottom', axis.title.x = element_blank())
  
  av_plot 
  
  return(av_plot)
  
}
