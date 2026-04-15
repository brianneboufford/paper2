# ------------------------------------------------------------------------------
# metric difference figures for paper 2
#
# 
# 
# February 25, 2026
# updated march 22, 2026
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
library(viridis)
library(stringr)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios_Feb19 - Copy", "Runs")
outpath <- file.path(".", "data", "streamflow_analysis")

fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Mar22")

# list of hydrographs
hydrograph_list <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                              recursive = TRUE, 
                              full.name = TRUE)

results_pf <- lapply(hydrograph_list, 
                     get_peak_flows_data) %>%
  do.call(rbind, .)

results_ex <- lapply(hydrograph_list, 
                     get_extreme_flows_data) %>%
  do.call(rbind, .)

# -----------------------------------------------------------------------------
# clean outputs 
# -----------------------------------------------------------------------------
results_pf <- results_pf %>% 
  subset(!Year %in% c(1980, 2023))

# box pot for Max_1_Day, Max_1_Day_Doy and Max_7Day, Max_7_Day_Doy 
pf_af <- results_pf[results_pf$run_name == "all_forest", ]
pf_no_af <- results_pf[results_pf$run_name != "all_forest", ]
pf_no_af <- pf_no_af[!grepl("_40_", pf_no_af$run_name), ]
#pf_no_af <- pf_no_af[!grepl("50yrs", pf_no_af$run), ]

pf_no_af$sim_name <- pf_no_af$run_name %>%
  str_remove("_(\\d+?)yrs") %>% str_replace("_", " ") %>% paste0(., "%")

pf_no_af$recovery <- pf_no_af$run_name %>%
  str_remove("_(\\d+?)_") %>% str_replace("high", "") %>% str_replace("low", "")

pf_no_af <- pf_no_af %>% 
  mutate(recovery =  
           factor(recovery, 
                  levels = c("50yrs", "45yrs", "40yrs", "35yrs",
                             "30yrs", "25yrs", "20yrs", "15yrs", "10yrs",
                             "5yrs", "0yrs"), 
                  labels = c("AF", "46-50", "41-45", "36-40", "31-35",
                             "26-30", "21-25", "16-20", "11-15", "6-10",
                             "0-5"
                  )))
pf_no_af <- pf_no_af %>% 
  mutate(scenario = factor(sim_name, 
                           levels = c("low 10%", "low 15%", "low 20%", "low 30%", 
                                      "high 10%", "high 15%", "high 20%", "high 30%"), 
                           labels = c(expression("low elevation":10*'%'), 
                                      expression("low elevation":15*'%'), 
                                      expression("low elevation":20*'%'), 
                                      expression("low elevation":30*'%'), 
                                      expression("high elevation":10*'%'),
                                      expression("high elevation":15*'%'), 
                                      expression("high elevation":20*'%'), 
                                      expression("high elevation":30*'%'))))
# -----------------------------------------------------------------------------
# clean ex output 
results_ex <- results_ex %>% 
  subset(!Year %in% c(1980, 2023))

# box pot for Max_1_Day, Max_1_Day_Doy and Max_7Day, Max_7_Day_Doy 
ex <- results_ex[results_ex$run_name == "all_forest", ]
x_no_af <- results_ex[results_ex$run_name != "all_forest", ]
x_no_af <- x_no_af[!grepl("_40_", x_no_af$run_name), ]
#pf_no_af <- pf_no_af[!grepl("50yrs", pf_no_af$run), ]

x_no_af$sim_name <- x_no_af$run_name %>%
  str_remove("_(\\d+?)yrs") %>% str_replace("_", " ") %>% paste0(., "%")

x_no_af$recovery <- x_no_af$run_name %>%
  str_remove("_(\\d+?)_") %>% str_replace("high", "") %>% str_replace("low", "")

x_no_af <- x_no_af %>% 
  mutate(recovery =  
           factor(recovery, 
                  levels = c("50yrs", "45yrs", "40yrs", "35yrs",
                             "30yrs", "25yrs", "20yrs", "15yrs", "10yrs",
                             "5yrs", "0yrs"), 
                  labels = c("AF", "46-50", "41-45", "36-40", "31-35",
                             "26-30", "21-25", "16-20", "11-15", "6-10",
                             "0-5"
                  )))
x_no_af <- x_no_af %>% 
  mutate(scenario = factor(sim_name, 
                           levels = c("low 10%", "low 15%", "low 20%", "low 30%", 
                                      "high 10%", "high 15%", "high 20%", "high 30%"), 
                           labels = c(expression("low elevation":10*'%'), 
                                      expression("low elevation":15*'%'), 
                                      expression("low elevation":20*'%'), 
                                      expression("low elevation":30*'%'), 
                                      expression("high elevation":10*'%'),
                                      expression("high elevation":15*'%'), 
                                      expression("high elevation":20*'%'), 
                                      expression("high elevation":30*'%'))))
# ------------------------------------------------------------------------------
# PLOT

max_day <- ggplot(data = pf_no_af, aes(x = scenario, y = Max_1_Day_DoY, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_fill_manual(values = c("0-5" = "#8e0152", 
                               "6-10" = "#c51b7d", 
                               "11-15" = "#de77ae", 
                               "16-20" = "#f1b6da", 
                               "21-25" = "#fde0ef", 
                               "26-30" = "#f7f7f7", 
                               "31-35" = "#e6f5d0", 
                               "36-40" = "#b8e186", 
                               "41-45" = "#7fbc41", 
                               "46-50" = "#4d9221", 
                               "AF" = "#969696")) +
  labs(x = "Simulation",
       y = "Peak 1-Day Streamflow DoY",
       fill = "Recovery") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "right") + 
  guides(fill = guide_legend(reverse = TRUE))

max_day_x <- ggplot(data = x_no_af, aes(x = scenario, y = Max_1_Day_DoY, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_fill_manual(values = c("0-5" = "#8e0152", 
                               "6-10" = "#c51b7d", 
                               "11-15" = "#de77ae", 
                               "16-20" = "#f1b6da", 
                               "21-25" = "#fde0ef", 
                               "26-30" = "#f7f7f7", 
                               "31-35" = "#e6f5d0", 
                               "36-40" = "#b8e186", 
                               "41-45" = "#7fbc41", 
                               "46-50" = "#4d9221", 
                               "AF" = "#969696")) +
  labs(x = "Simulation",
       y = "Peak 1-Day Streamflow DoY",
       fill = "Recovery") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "right") + 
  guides(fill = guide_legend(reverse = TRUE))

x_no_af_mean <- x_no_af %>%
  group_by(scenario, recovery) %>% 
  summarize(mean_doy = round(mean(Max_1_Day_DoY), 0))

# ------------------------------------------------------------------------------
# test 
# ------------------------------------------------------------------------------
av_peak_doy_x <- x_no_af %>% 
  group_by(sim_name, recovery) %>% 
  summarize(av_doy = median(Max_1_Day_DoY), 0)

av_peak_doy_x <- av_peak_doy_x %>% 
  mutate(scenario = factor(sim_name, 
                           levels = c("low 10%", "low 15%", "low 20%", "low 30%", 
                                      "high 10%", "high 15%", "high 20%", "high 30%"), 
                           labels = c(expression("low elevation":10*'%'), 
                                      expression("low elevation":15*'%'), 
                                      expression("low elevation":20*'%'), 
                                      expression("low elevation":30*'%'), 
                                      expression("high elevation":10*'%'),
                                      expression("high elevation":15*'%'), 
                                      expression("high elevation":20*'%'), 
                                      expression("high elevation":30*'%'))))

max_day_x <- ggplot(data = av_peak_doy_x, aes(x = scenario, y = av_doy, colour = "recovery")) +
  geom_point() + 
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_manual(values = c("0-5" = "#8e0152",
                               "6-10" = "#c51b7d",
                               "11-15" = "#de77ae",
                               "16-20" = "#f1b6da",
                               "21-25" = "#fde0ef",
                               "26-30" = "#f7f7f7",
                               "31-35" = "#e6f5d0",
                               "36-40" = "#b8e186",
                               "41-45" = "#7fbc41",
                               "46-50" = "#4d9221",
                               "AF" = "#969696")) +
  labs(x = "Simulation",
       y = "Peak 1-Day Streamflow DoY",
       colour = "recovery") +
  theme_minimal() +
  coord_flip() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "right") + 
  guides(fill = guide_legend(reverse = TRUE))

# ------------------------------------------------------------------------------
# copied from 02_get_streamflow_metrics 
# ------------------------------------------------------------------------------
get_peak_flows_data <- function(hydro_file){
  
  run_name <- sub(".*?_", "", basename(dirname(hydro_file)))
  
  hydro <- read.csv(hydro_file)
  
  # convert date column to date format 
  hydro$date <- as.Date(hydro$date)
  hydro$mean <- mean(hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  hydro_clean <- hydro 
  hydro_clean <- hydro_clean[!(is.na(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)), ]
  names(hydro_clean) <- c("time", "date", "hour", "precip", "modelled", "observed")
  
  years_hydro <- year(hydro_clean$date)
  first_year <- min(years_hydro)
  last_year <- max(years_hydro)
  start_year <- first_year + 1
  end_year <- last_year - 1 
  
  ###############################
  # peak flows
  peak_flows <- calc_annual_highflows(data = hydro_clean,
                                      dates = date, 
                                      values = modelled,
                                      start_year = start_year, 
                                      end_year = end_year, 
                                      ignore_missing = TRUE) %>% 
    as.data.frame()
  
  peak_flows$run_name <- run_name
  
  return(peak_flows)
}

# ------------------------------------------------------------------------------
get_extreme_flows_data <- function(hydro_file){
  
  run_name <- sub(".*?_", "", basename(dirname(hydro_file)))
  
  hydro <- read.csv(hydro_file)
  
  # convert date column to date format 
  hydro$date <- as.Date(hydro$date)
  hydro$mean <- mean(hydro$TRAPPING_CREEK_NEAR_THE_MOUTH..m3.s.)
  
  hydro_clean <- hydro 
  hydro_clean <- hydro_clean[!(is.na(hydro_clean$TRAPPING_CREEK_NEAR_THE_MOUTH..observed...m3.s.)), ]
  names(hydro_clean) <- c("time", "date", "hour", "precip", "modelled", "observed")
  
  years_hydro <- year(hydro_clean$date)
  first_year <- min(years_hydro)
  last_year <- max(years_hydro)
  start_year <- first_year + 1
  end_year <- last_year - 1 
  
  ###############################
  # peak flows
  ann_extremes_full_yr <- calc_annual_extremes(data = hydro_clean,
                                               dates = date,
                                               values = modelled,
                                               start_year = start_year,
                                               end_year = end_year, 
                                               ignore_missing = TRUE)
  
  ann_extremes_full_yr$run_name <- run_name
  
  return(ann_extremes_full_yr)
}