# ------------------------------------------------------------------------------
# assess modelled SWE 
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

weather_data_path <- file.path(".", "data", "weather_validation_data")
snow_summary_path <- file.path(weather_data_path, "snow_summary")
weather_file_path <- file.path(weather_data_path, "climate-daily.csv")
snow_run_path <- file.path(".", "raven-runs", "Baseline3-mar22", "Runs", "Trapping_snow_baseline_mar22") 

swe_list <- list.files(snow_run_path,
                       pattern = "SNOW",
                       recursive = TRUE,
                       full.name = TRUE)

set.seed(1113)

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

# swe files 
man_swe_files <- list.files(file.path(weather_data_path, "snow_data", "manual_weather_stations", "swe"), full.names = TRUE)
auto_swe_files <- list.files(file.path(weather_data_path, "snow_data", "automated_weather_stations", "swe"), full.names = TRUE)

# station id key
man_key <- readxl::read_excel(file.path(weather_data_path, "snow_data", "manual_weather_stations", "manual_snow_station_key.xlsx"))
auto_key <- readxl::read_excel(file.path(weather_data_path, "snow_data", "automated_weather_stations", "automated_snow_station_key.xlsx"))

# -----------------------------------------------------------------------------
# apply reformat to manual snow data 
manual_snow_data <- lapply(man_swe_files, 
                           reformat_man_swe, 
                           man_key = man_key) %>%
  do.call(rbind,.)

# apply reformat to auto snow data
auto_snow_data <- lapply(auto_swe_files,
                         reformat_auto_swe,
                         auto_key = auto_key) %>%
  do.call(rbind, .)

# sample only once a month from auto stations to avoid overfitting to the auto
# stations 
auto_snow_data$date <- as.Date(auto_snow_data$date)
auto_snow_data$month <- lubridate::month(auto_snow_data$date)
auto_snow_data$year <- lubridate::year(auto_snow_data$date)

sampled_auto_data <- auto_snow_data %>%
  group_by(year, month, LOCATION_NAME) %>%
  slice_sample(n = 1) %>%  # randomly select one row per group
  ungroup() %>%
  dplyr::select(-c(month, year))

# * important
auto_snow_data <- sampled_auto_data

#auto_snow_data <- auto_snow_data %>% 
#  select(-c(month, year))

swe_metrics <- compare_LR_swe(swe_list, 
                              manual_snow_data,
                              auto_snow_data)

# go into compare_LR_swe and snow_year_plot_functions and run each line !!! 
snow_file_path <- swe_list 

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

# function to reformat manual snow station data and join to station meta data
reformat_man_swe <- function(man_swe_file, man_key){
  
  str_1 <- stringr::str_split(man_swe_file, "@")[[1]][2]
  id <- stringr::str_split(str_1, "-")[[1]][1]
  loc_meta <- filter(man_key, LOCATION_ID == id) %>%
    as.data.frame()
  loc_name <- loc_meta$LOCATION_NAME
  
  man_swe <- read.csv(man_swe_file, skip=2)
  names(man_swe) <- c("date", "meas_mm_swe")
  man_swe$date <- as.Date(man_swe$date)
  man_swe$LOCATION_NAME <- loc_name
  
  man_swe_meta <- full_join(man_swe, loc_meta, by="LOCATION_NAME")
  man_swe_meta$data_type <- "manual"
  
  return(man_swe_meta)
  
}

# same function for automatic snow stations 
reformat_auto_swe <- function(auto_swe_file, auto_key){
  
  str_1 <- stringr::str_split(auto_swe_file, "@")[[1]][2]
  id <- stringr::str_split(str_1, "-")[[1]][1]
  loc_meta <- filter(auto_key, LOCATION_ID == id) %>%
    as.data.frame()
  loc_name <- loc_meta$LOCATION_NAME
  
  auto_swe <- read.csv(auto_swe_file, skip=2)
  names(auto_swe) <- c("date", "meas_mm_swe")
  auto_swe$date <- as.Date(auto_swe$date)
  auto_swe$LOCATION_NAME <- loc_name 
  
  auto_swe_meta <- full_join(auto_swe, loc_meta, by="LOCATION_NAME")
  auto_swe_meta$data_type <- "auto"
  
  return(auto_swe_meta)
  
}

compare_LR_swe <- function(swe_list, 
                           manual_snow_data,
                           auto_snow_data){
  
  # ---------------------------------------------------
  # inputs: 
  # @param swe_list : file path to raven snow output 
  # @param manual_snow_data : manual snow station snow measurements 
  # @param auto_snow_data : automatic snow station snow measurements
  #
  # outputs: 
  # dataframe with R2 and pbais metrics for each station 
  # and snow parameters for the associated raven run 
  # ----------------------------------------------------
  
  # Split the remaining string by "_"
  parts <- strsplit(swe_list, "/")[[1]]
  
  file_id <- parts[9]
  
  raven_swe <- read.csv(swe_list, skip = 1) %>%# [mm] frozen snow depth (mm SWE : snow water equivalent)
    dplyr::select(-c("time", "X"))
  raven_swe$day <- as.Date(raven_swe$day)
  hru_names_ordered <- c("GRANO_CREEK", "MISSION_CREEK","CARMI", "BIG_WHITE", "MCCHULLOCH_MANUAL")
  names(raven_swe)[2:length(names(raven_swe))] <- paste0(hru_names_ordered) 
  
  manual_snow_data <- manual_snow_data %>%
    dplyr::select(-c(SNOW_MSS_LOC_ID, STATUS))
  manual_snow_data$LOCATION_NAME <- manual_snow_data$LOCATION_NAME %>% 
    stringr::str_replace("Aberdeen Lake", "ABERDEEN") %>%
    stringr::str_replace("Monashee Pass", "MONASHEE") %>%
    stringr::str_replace("Carmi", "CARMI") %>%
    stringr::str_replace("Big White Mountain", "BIG_WHITE") %>%
    stringr::str_replace("McCulloch", "MCCHULLOCH_MANUAL") %>%
    stringr::str_replace("Graystoke Lake", "GRAYSTOKE") %>%
    stringr::str_replace("Oyama Lake", "OYAMA_MANUAL") %>%
    stringr::str_replace("Vaseux Creek", "VASEAUX")
  
  auto_snow_data <- auto_snow_data %>%
    dplyr::select(-c(OBJECT_ID, STATUS))
  auto_snow_data$LOCATION_NAME <- auto_snow_data$LOCATION_NAME %>% 
    stringr::str_replace("Grano Creek", "GRANO_CREEK") %>%
    stringr::str_replace("Mission Creek", "MISSION_CREEK") %>%
    stringr::str_replace("McCulloch", "MCCHULLOCH") %>%
    stringr::str_replace("Greyback Reservoir", "GREYBACK_RESERVOIR") %>%
    stringr::str_replace("Oyama Lake", "OYAMA") 
  
  all_snow_data <- rbind(auto_snow_data, manual_snow_data)
  
  long_raven <- raven_swe %>%
    pivot_longer(
      cols = names(raven_swe)[2:6],
      names_to = "LOCATION_NAME",
      values_to = "Raven_swe_mm"
    )
  
  names(long_raven)[names(long_raven) == "day"] <- "date"
  long_raven$date <- as.Date(long_raven$date)
  
  joined_data <- left_join(long_raven, all_snow_data,
                           by=c("date", "LOCATION_NAME")) %>%
    as.data.frame() %>%
    na.omit()
  
  joined_data$month <- month(joined_data$date)
  joined_data$year <- year(joined_data$date)
  
  joined_data$LOCATION_NAME[joined_data$LOCATION_NAME == "OYAMA_MANUAL"] <- "OYAMA"
  joined_data <- joined_data[joined_data$LOCATION_NAME != "VASEAUX", ]
  
  joined_data <- joined_data[!month(joined_data$date) %in% c("7", "8", "9"), ]
  
  # looking at R2 and Pbias for each station 
  snow_results <- joined_data %>% group_by(LOCATION_NAME) %>%
    summarize(swe_r2 = round(hydroGOF::R2(obs = meas_mm_swe, sim = Raven_swe_mm),2),
              swe_pbias = hydroGOF::pbias(sim = Raven_swe_mm, 
                                          obs = meas_mm_swe),
              N = length(meas_mm_swe),
              ELEVATION = first(ELEVATION)) 
  
  snow_results$model_run <- file_id
  
  return(snow_results)
  
}

snow_year_plot <- function(snow_file_path, all_snow_data, model_run){
  # Raven snow data___________________________________________________
  
  # read raven swe 
  raven_swe <- read.csv(snow_file_path, skip = 1) %>% # [mm] frozen snow depth (mm SWE : snow water equivalent)
    dplyr::select(-c("time", "X"))
  raven_swe$day <- as.Date(raven_swe$day)
  hru_names_ordered <- c("GRANO_CREEK", "MISSION_CREEK","CARMI", "BIG_WHITE", "MCCHULLOCH_MANUAL")
  names(raven_swe)[2:length(names(raven_swe))] <- paste0(hru_names_ordered) 
  raven_swe$year <- lubridate::year(raven_swe$day)
  
  # convert to long format and fix names
  raven_long <- raven_swe %>%
    pivot_longer(cols = !starts_with(c("day", "year")), 
                 names_to = "Station", 
                 values_to = "SWE")
  names(raven_long) <- c("date", "year", "LOCATION_NAME", "SWE")
  raven_long[raven_long$LOCATION_NAME == "OYAMA_MANUAL", "LOCATION_NAME"] <- "OYAMA"
  
  # fix auto snow data dates
  all_snow_data$date <- as.Date(all_snow_data$date)
  all_snow_data$year <- lubridate::year(all_snow_data$date)
  
  joined_data <- joined_data[!month(joined_data$date) %in% c("7", "8", "9"), ]
  # merge data together
  raven_all <- merge(raven_long, all_snow_data, by=c("date", "LOCATION_NAME", "year"))
  
  raven_all <- raven_all[!month(raven_all$date) %in% c("7", "8", "9"), ]
  
  # convert to snow year instead of calendar year
  snow_year_df <- raven_all %>%
    mutate(
      date = as.Date(date),  # ensure it's in Date format
      snow_year = if_else(month(date) >= 8, year(date), year(date) - 1)
    )
  
  # arrange by snow year date
  snow_year_df <- snow_year_df %>%
    arrange(snow_year, month(date), day(date))
  snow_year_df$month <- lubridate::month(snow_year_df$date)
  
  # filter for specific snow year
  #snow_year_df <- snow_year_df %>% filter(snow_year == yr)
  
  snow_year_auto <- snow_year_df %>% filter(data_type == "auto")
  snow_year_manual <- snow_year_df %>% filter(data_type == "manual")
  
  # ------------------------------------------------------------------------------
  # Figure V paper 
  # -----------------------------------------------------------------------------
  
  month_order <- c(8:12, 1:7)
  month_labels <- month.abb[month_order]
  
  swe_summary <- snow_year_df %>%
    group_by(LOCATION_NAME, month) %>%
    summarise(
      Modelled = mean(SWE, na.rm = TRUE),
      Measured = mean(meas_mm_swe, na.rm = TRUE),
      sd_Modelled = sd(SWE, na.rm = TRUE),
      sd_Measured = sd(meas_mm_swe, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(month = factor(month, levels = month_order, labels = month_labels, ordered = TRUE))
  
  swe_long <- swe_summary %>%
    pivot_longer(cols = c(Modelled, Measured), names_to = "Type", values_to = "Mean") %>%
    mutate(
      SD = ifelse(Type == "Modelled", sd_Modelled, sd_Measured),
      ymin = Mean - SD,
      ymax = Mean + SD
    )
  
  swe_long <- swe_long[swe_long$LOCATION_NAME != "VASEAUX", ]
  
  swe_long$month <- factor(swe_long$month, 
                           levels = c("Aug", "Sep", "Oct", "Nov", "Dec",
                                      "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"))
  
  station_labels <- c(
   # ABERDEEN = "Aberdeen",
    BIG_WHITE = "Big White",
    CARMI = "Carmi",
    GRANO_CREEK = "Grano Creek",
  #  GRAYSTOKE = "Graystoke", 
  #  GREYBACK_RESERVOIR = "Greyback Reservoir",
    MCCHULLOCH_MANUAL = "Mcchulloch", 
    MISSION_CREEK = "Mission Creek"#,
  #  MONASHEE = "Monashee", 
  #  OYAMA = "Oyama"
  )
  
  
  # Use this in the labeller
  snow_plot <- ggplot(swe_long, aes(x = month, y = Mean, color = Type, fill = Type, group = Type)) +
    geom_line() +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    facet_wrap(~ LOCATION_NAME, labeller = labeller(LOCATION_NAME = station_labels)) +
    labs(x = "Month (Augb