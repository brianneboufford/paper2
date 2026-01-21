# ------------------------------------------------------------------------------
# assess modelled streamflow - baseline
#
# 
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

hru_run_path <- file.path(".", "raven-runs", "Baseline2", "Runs", "Trapping_hru_baseline_Dec2") 

# list of hydrographs
hydrograph_list <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                              recursive = TRUE, 
                              full.name = TRUE)

run_metrics <- evalutate_hydrograph_performance(hydrograph_list)

hydro <- read.csv(hydrograph_list)

evaluate_Q_performance(hydro, "baseline")
hydro <- prep_model_data(hydro, "baseline") 

model_data <- rbind(hydro)
daily_average <- plot_average_daily(model_data)

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

prep_model_data <- function(df, data_name){
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

# FUNCTION TO CLEAN HYDROGRAPH DATA
clean_hydrograph <- function(hydro){
  
  hydro$date <- as.Date(hydro$date)
  hydro$mean <- mean(hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  hydro_clean <- hydro 
  hydro_clean <- hydro_clean[!(is.na(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)), ]
  
  # difference between modeled and observed streamflow 
  hydro$difference <- hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s. - hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.
  
  hydro$pos <- TRUE
  hydro$pos[hydro$difference < 0] <- FALSE
  
  return(hydro)
  
}

# function to get RMSE, R2, NSE, KGE, and pbias from raven runs 
evalutate_hydrograph_performance <- function(hydrograph_list){
  
  # get name of Raven run 
  model_run <- stringr::str_split(hydrograph_list, pattern = "/")[[1]]
  model_run <- model_run[2]
  
  # read hydrograph csv
  hydro <- read.csv(file.path(hydrograph_list))
  
  # convert date column to date format 
  hydro$date <- as.Date(hydro$date)
  hydro$mean <- mean(hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  hydro_clean <- hydro 
  hydro_clean <- hydro_clean[!(is.na(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)), ]
  
  # caluclaute RMSE and R2 for the time series 
  rootmeanse <- rmse(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  r2 <- cor(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.) ^ 2
  
  kge <- GOF_kling_gupta_efficiency(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                    obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                    na.rm = TRUE)
  
  nse <- GOF_nash_sutcliffe_efficiency(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                       obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                       j = 2, na.rm = TRUE)
  
  pbias_hyd <- GOF_percent_bias(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                na.rm=TRUE)
  
  results <- data.frame(model_run = model_run, 
                        RMSE = rootmeanse,
                        R2 = r2,
                        KGE = kge, 
                        NSE = nse, 
                        pbias = pbias_hyd)
  
  return(results)
}

plot_average_daily <- function(model_data, data_type = 'Streamflow', Qup = 0.9, Qlow = 0.1, wrap_gen = Inf, drop_nas = T){
  
  # Specify Lab
  if(data_type == 'Streamflow'){  y_lab = expression(Streamflow~(m^3/s))}else
    if(data_type == 'Stage'){       y_lab = 'Water Level (m)'}else
      if(data_type == 'Temperature'){ y_lab = expression(Temperature~(degree~C))}else{
        y_lab = data_type }
  
  model_data_grouped <- model_data %>%
    #{if(drop_nas) spread(., Type, Value) %>% drop_na() %>% gather(., Type, Value, -Site, -date) else . } %>%
    group_by(Site, date = yday(date), Type) %>%
    #mutate(date = yday(date)) %>%
    group_by(Site, date, Type) %>%
    summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
              Qlow = quantile(Value, Qlow, na.rm = TRUE),
              Value = mean(Value, na.rm = TRUE)) %>%
    ungroup() 
  
  av_plot <- ggplot(data = model_data_grouped, aes(x = as.Date(strptime(date, format = '%j')), y = Value,
                                                   colour = Type, fill = Type)) +
    geom_ribbon(aes(ymin = Qlow, ymax = Qhigh), alpha = 0.25, colour = NA) +
    geom_line() +
    facet_wrap(~forcats::fct_reorder(Site, Value, .fun = function(x){mean(x, na.rm = TRUE)}),
               scales = 'free_y', labeller = label_wrap_gen(width = wrap_gen)) +
    scale_x_date(date_labels = '%b') +
    scale_colour_manual('', values = c( 'navy', 'firebrick')) +
    scale_fill_manual('', values = c( 'navy', 'firebrick')) +
    labs(y = y_lab) +
    theme_bw(14) +
    theme(legend.position = 'bottom', axis.title.x = element_blank())
  
  av_plot 
  
  return(av_plot)
  
}

# function to get RMSE, R2, NSE, KGE, and pbias from raven runs 
evaluate_Q_performance <- function(hydro, model_name){
  
  # -------------------------------------------------------------------
  # overall performance metrics 
  # -------------------------------------------------------------------
  
  # convert date column to date format 
  hydro$date <- as.Date(hydro$date)
  hydro$mean <- mean(hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  hydro_clean <- hydro 
  hydro_clean <- hydro_clean[!(is.na(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)), ]
  
  hydro_clean <- hydro_clean[hydro_clean$date > as.Date("1980-12-31"), ]
  # hydro_clean <- hydro_clean[hydro_clean$date < as.Date("2023-01-01"), ]
  
  # caluclaute RMSE and R2 for the time series 
  rootmeanse <- rmse(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  r2 <- cor(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.) ^ 2
  
  kge <- GOF_kling_gupta_efficiency(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                    obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                    na.rm = TRUE)
  
  nse <- GOF_nash_sutcliffe_efficiency(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                       obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                       j = 2, na.rm = TRUE)
  
  pbias_hyd <- GOF_percent_bias(mod = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.,
                                obs = hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.,
                                na.rm=TRUE)
  
  # -------------------------------------------------------------------
  # peak flow metrics 
  # -------------------------------------------------------------------
  
  # subset data to just peak value for each year
  df_peaks <- hydro_clean %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarise(
      date_obs_peak = date[which.max(TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)],
      observed_peak = max(TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s., na.rm = TRUE),
      date_mod_peak = date[which.max(TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)],
      modelled_peak = max(TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s., na.rm = TRUE),
    ) %>%
    ungroup()
  
  # plot observed peak vs modelled peak
  # ggplot() + 
  #   geom_point(data=df_peaks, aes(x=observed_peak, y = modelled_peak))
  
  # plot observed peak date vs modelled peak date 
  # ggplot() + 
  #   geom_point(data=df_peaks, aes(x=lubridate::yday(date_obs_peak), y = lubridate::yday(date_mod_peak)))
  
  df_peaks <- df_peaks[!df_peaks$year %in% c(1980, 2023), ]
  r2_peak <- cor(df_peaks$observed_peak, df_peaks$modelled_peak) ^ 2
  
  rmse_peak <- rmse(df_peaks$observed_peak, df_peaks$modelled_peak)
  
  # double check r2 calculation
  # model <- lm(observed_peak ~ modelled_peak, data = df_peaks)
  # summary(model)$r.squared
  
  pbias_peaks <- GOF_percent_bias(mod = df_peaks$modelled_peak,
                                  obs = df_peaks$observed_peak,
                                  na.rm=TRUE)
  
  
  
  # -------------------------------------------------------------------
  # MAF metrics + mean august - september + May + June 
  # -------------------------------------------------------------------
  # not doing at this moment but could include later
  
  results <- data.frame(model_run = model_name, 
                        RMSE = rootmeanse,
                        R2 = r2,
                        KGE = kge, 
                        NSE = nse, 
                        pbias = pbias_hyd,
                        RMSE_peak = rmse_peak,
                        R2_peak = r2_peak, 
                        pbias_peak = pbias_peaks)
  
  return(results)
}
