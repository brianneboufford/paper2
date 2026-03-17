# ------------------------------------------------------------------------------
# get streamflow results for paper 2 
# adapted from 02_get_streamflow_metrics.R
#
# 
# 
# March 9, 2025
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
library(patchwork)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios_Feb19", "Runs")
outpath <- file.path(".", "data", "streamflow_analysis")

fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Feb20")

# list of hydrographs
hydrograph_list <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                              recursive = TRUE, 
                              full.name = TRUE)

# ------------------------------------------------------------------------------
# get flood distribution and flood frequency data 
# ------------------------------------------------------------------------------ 

results_pf <- lapply(hydrograph_list, 
                  get_peak_flows_data) %>%
  do.call(rbind, .)

results_pf_keep <- results_pf

results_ff <- lapply(hydrograph_list, 
                     get_flood_freq_data) %>% 
  do.call(rbind,. )

results_cf <- lapply(hydrograph_list, 
                     get_cumulative_flow) %>% 
  do.call(rbind,. )

# ------------------------------------------------------------------------------
# make plots 
# ------------------------------------------------------------------------------

results_pf <- results_pf %>% 
  subset(!Year %in% c(1981, 2023))

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
                           labels = c(expression(z<z[p50]:10*'%'), 
                                      expression(z<z[p50]:15*'%'), 
                                      expression(z<z[p50]:20*'%'), 
                                      expression(z<z[p50]:30*'%'), 
                                      expression(z>=z[p50]:10*'%'),
                                      expression(z>=z[p50]:15*'%'), 
                                      expression(z>=z[p50]:20*'%'), 
                                      expression(z>=z[p50]:30*'%'))))

max_peak <- ggplot(data = pf_no_af, aes(x = scenario, y = Max_1_Day, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_fill_manual(values = c("0yrs" = "#8e0152", 
                               "5yrs" = "#c51b7d", 
                               "10yrs" = "#de77ae", 
                               "15yrs" = "#f1b6da", 
                               "20yrs" = "#fde0ef", 
                               "25yrs" = "#f7f7f7", 
                               "30yrs" = "#e6f5d0", 
                               "35yrs" = "#b8e186", 
                               "40yrs" = "#7fbc41", 
                               "45yrs" = "#4d9221", 
                               "AF" = "#969696")) +
  labs(x = NULL,
       y = expression("Peak 1-Day Streamflow ("*m^3/s*")"),
       fill = "Recovery") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none")


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

combo_plot <- plot_grid(max_peak, max_day, 
                        ncol=1, align = "v", rel_heights = c(1, 1.5))

ggsave(max_day,
       filename = file.path(fig_path, "Peak_DoY_mar16.png"),
       units = "in",
       dpi = 300,
       width = 6,
       height = 6)

summary_doy_stats <- pf_no_af %>% 
  group_by(sim_name, recovery) %>% 
  summarize(med_DoY = median(Max_1_Day_DoY))

summary_2_doy_stats <- pf_af %>% 
  group_by(run_name) %>% 
  summarize(med_DoY = median(Max_1_Day_DoY))

ggsave(combo_plot,
       filename = file.path(fig_path, "Peak_1_day_Q_DoY.png"),
       units = "in",
       dpi = 300,
       width = 9,
       height = 7)

ggplot(data = pf_no_af, aes(x = sim_name, y = Max_7_Day, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) + 
  scale_fill_manual(values = c("0yrs" = "#8e0152", 
                               "5yrs" = "#c51b7d", 
                               "10yrs" = "#de77ae", 
                               "15yrs" = "#f1b6da", 
                               "20yrs" = "#fde0ef", 
                               "25yrs" = "#f7f7f7", 
                               "30yrs" = "#e6f5d0", 
                               "35yrs" = "#b8e186", 
                               "40yrs" = "#7fbc41", 
                               "45yrs" = "#4d9221", 
                               "AF" = "#969696")) + 
  labs(x = "Simulation", 
       y = expression("Peak 7-day Streamflow ("*m^3/s*")"), 
       legend = "Recovery") + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "right")

ggplot(data = pf_no_af, aes(x = sim_name, y = Max_7_Day_DoY, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) + 
  scale_fill_manual(values = c("0yrs" = "#8e0152", 
                               "5yrs" = "#c51b7d", 
                               "10yrs" = "#de77ae", 
                               "15yrs" = "#f1b6da", 
                               "20yrs" = "#fde0ef", 
                               "25yrs" = "#f7f7f7", 
                               "30yrs" = "#e6f5d0", 
                               "35yrs" = "#b8e186", 
                               "40yrs" = "#7fbc41", 
                               "45yrs" = "#4d9221", 
                               "AF" = "#969696")) + 
  labs(x = "Simulation", 
       y = expression("Peak 7-Day Streamflow DoY"), 
       legend = "Recovery") + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "right")

# ------------------------------------------------------------------------------
# make cumulative flow plot 
# ------------------------------------------------------------------------------

# box pot for Max_1_Day, Max_1_Day_Doy and Max_7Day, Max_7_Day_Doy 
cf_af <- results_cf[results_cf$run_name == "all_forest", ]
cf_no_af <- results_cf[results_cf$run_name != "all_forest", ]
cf_no_af <- cf_no_af[!grepl("_40_", cf_no_af$run_name), ]
#pf_no_af <- pf_no_af[!grepl("50yrs", pf_no_af$run), ]

cf_no_af$sim_name <- cf_no_af$run_name %>%
  str_remove("_(\\d+?)yrs") %>% str_replace("_", " ") %>% paste0(., "%")

cf_no_af$recovery <- cf_no_af$run_name %>%
  str_remove("_(\\d+?)_") %>% str_replace("high", "") %>% str_replace("low", "")

cf_no_af <- cf_no_af %>% 
  mutate(recovery =  
           factor(recovery, 
                  levels = c("0yrs", "5yrs", "10yrs", "15yrs",
                             "20yrs", "25yrs", "30yrs", "35yrs", "40yrs",
                             "45yrs", "50yrs"), 
                  labels = c("0-5", "6-10", "11-15", "16-20",
                             "21-25", "26-30", "31-35", "36-40", "41-45",
                             "46-50", "AF")))

cf_no_af <- cf_no_af %>% 
  mutate(scenario = factor(sim_name, 
                           levels = c("high 10%", "high 15%", "high 20%", "high 30%", 
                                      "low 10%", "low 15%", "low 20%", "low 30%"), 
                           labels = c(expression(z>=z[p50]:10*'%'),
                                      expression(z>=z[p50]:15*'%'), 
                                      expression(z>=z[p50]:20*'%'), 
                                      expression(z>=z[p50]:30*'%'), 
                                      expression(z<z[p50]:10*'%'), 
                                      expression(z<z[p50]:15*'%'), 
                                      expression(z<z[p50]:20*'%'), 
                                      expression(z<z[p50]:30*'%'))))

cf_no_af_apr_jul <- cf_no_af[cf_no_af$Month %in% c("Apr", "May", "Jun", "Jul", "Aug", "Dec"), ]

cf_res_apr <- cf_no_af_apr_jul[cf_no_af_apr_jul$Month == "Apr", ]
cf_res_tot <- cf_no_af_apr_jul[cf_no_af_apr_jul$Month == "Dec", ]

cf_res_apr <- cf_res_apr %>%
  group_by(sim_name) %>%
  mutate(percent_diff = round(Mean/Mean[recovery == "AF"]*100, 2)) %>% 
  ungroup()

cf_res_tot <- cf_res_tot %>%
  group_by(sim_name) %>%
  mutate(percent_diff = round(Mean/Mean[recovery == "AF"]*100, 2)) %>% 
  ungroup()

ggplot(data = cf_no_af_apr_jul, mapping = aes(x = Month, y = Mean, colour = recovery)) + 
  geom_point(size = 2, alpha = 0.5) +
  facet_wrap(~sim_name, ncol = 4) + 
  scale_colour_manual(values = c("0yrs" = "#8e0152", 
                               "5yrs" = "#c51b7d", 
                               "10yrs" = "#de77ae", 
                               "15yrs" = "#f1b6da", 
                               "20yrs" = "#fde0ef", 
                               "25yrs" = "#f7f7f7", 
                               "30yrs" = "#b8e186", 
                               "35yrs" = "#7fbc41", 
                               "40yrs" = "#4d9221", 
                               "45yrs" = "#276419", 
                               "AF" = "#969696")) + 
  labs(x = "Simulation", 
       y = expression("Mean Cumulative Flow"), 
       legend = "Recovery") + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "right")

cumulative_plot <- ggplot(cf_no_af_apr_jul, aes(x = Month, y = Mean)) +
  
  # regular recovery points
  geom_point(
    data = subset(cf_no_af_apr_jul, recovery != "AF"),
    aes(fill = recovery),
    shape = 21,
    colour = "black",
    size = 2.5,
    alpha = 0.7,
    stroke = 0.5
  ) +
  
  # AF star
  geom_point(
    data = subset(cf_no_af_apr_jul, recovery == "AF"),
    shape = 8,
    colour = "black",
    fill = "white",
    size = 3.5,
    stroke = 0.7,
    alpha = 0.7
  ) +
  scale_x_discrete(labels = function(x) ifelse(x == "Dec", "Total", x)) + 
  
  facet_wrap(~scenario, ncol = 4, labeller = as_labeller(label_parsed)) +
  
  scale_fill_manual(
    values = c("0-5" = "#8e0152", 
        "6-10" = "#c51b7d", 
        "11-15" = "#de77ae", 
        "16-20" = "#f1b6da", 
        "21-25" = "#fde0ef", 
        "26-30" = "#f7f7f7", 
        "31-35" = "#b8e186", 
        "36-40" = "#7fbc41", 
        "41-45" = "#4d9221", 
        "46-50" = "#276419", 
        "AF" = "#969696"),
    name = "Recovery"
  ) +
  
  labs(
    x = "Month",
    y = expression("Mean Cumulative Flow ("*m^3*")")
  ) +
  
  theme_bw(14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.background = element_rect(NA),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) + 
  guides(fill = guide_legend(nrow = 2))
  
ggsave(cumulative_plot,
       filename = file.path(fig_path, "cumulative_flow_plot_mar16.png"),
       units = "in",
       dpi = 300,
       width = 9,
       height = 7)
# ------------------------------------------------------------------------------
# generate metrics and save 
# ------------------------------------------------------------------------------ 

results <- lapply(hydrograph_list, 
                  calc_hydrograph_stats) %>%
  do.call(rbind, .)

write.csv(results, 
          file.path(outpath, "flow_metrics_mar10.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# make plot 
# ------------------------------------------------------------------------------ 

all_forested <- results[results$run == "all_forest", ]
result_no_af <- results[results$run != "all_forest", ]

joined_data <- dplyr::full_join(result_no_af, all_forested, by = "metrics")
names(joined_data) <- c("metrics", "values_run", "run", "values_af", "AF")
joined_data$value_diff <- joined_data$values_run - joined_data$values_af


# ---- 1. PARSE RUN NAME INTO scenario AND yrs ----
df <- joined_data %>%
  mutate(
    # extract "high_20" — everything before the last underscore
    scenario = sub("_(\\d+)yrs$", "", run),
    
    # extract years (20 from 20yrs)
    yrs = as.numeric(sub(".*_(\\d+)yrs$", "\\1", run))
  )

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------ 

# ------------------------------------------------------------------------
# function to get peak flow distribution data
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

# -----------------------------------------------------
# function to get peak flow distribution data
get_flood_freq_data <- function(hydro_file){
  
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
  
  
  # reformat
  peak_flows_1d <- dplyr::select(peak_flows, Year, Value = Max_1_Day)
  peak_flows_1d <- dplyr::mutate(peak_flows_1d, Measure = "1-Day")
  
  # compute freq analysis ** use_max = TRUE
  peak_freq_analysis <- compute_frequency_analysis(data = peak_flows_1d,
                                                   events = Year,
                                                   values = Value,
                                                   measures = Measure, 
                                                   use_max = TRUE)
  peak_quant <- peak_freq_analysis$Freq_Fitted_Quantiles %>% 
    as.data.frame()
  
  peak_quant$run_name <- run_name
  
  return(peak_quant)
}

# -----------------------------------------------------
# function to get peak flow distribution data
get_cumulative_flow <- function(hydro_file){
  
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
  
  ####################################################
  # calculate monthly cumulative stats
  
  cumulativesats <- calc_monthly_cumulative_stats(data = hydro_clean, 
                                dates = date,
                                values = modelled,
                                start_year = start_year, 
                                end_year = end_year)
  
  cumulativesats$run_name <- run_name
  
  return(cumulativesats)
}

# ---------------------------------------------------------------
calc_hydrograph_stats <- function(hydro_file){
  
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
  
  # --------------------------------------------------------------------
  # STATS TO CALCULATE 
  # - 7Q10
  # - 2 year peak flow
  # - 5 year peak flow
  # - 20 year peak flow
  # - mean annual flow
  # - mean aug-sept flow
  # ---------------------------------------------------------------------
  
  
  #############################################
  # 7Q10 - 7 day low flow 10 year return period
  low_flows <- calc_annual_lowflows(data = hydro_clean,
                                    dates = date, 
                                    values = modelled,
                                    roll_days = 7, 
                                    start_year = start_year,
                                    end_year = end_year,
                                    ignore_missing = TRUE)
  
  # reformat
  low_flows <- dplyr::select(low_flows, Year, Value = Min_7_Day)
  low_flows <- dplyr::mutate(low_flows, Measure = "7-Day")
  
  # compute freq analysis
  freq_analysis <- compute_frequency_analysis(data = low_flows,
                                              events = Year,
                                              values = Value,
                                              measures = Measure)
  sevQ10_quant <- freq_analysis$Freq_Fitted_Quantiles
  
  # get 7Q10
  sevQ10 <- sevQ10_quant$`7-Day`[3]
  
  ###############################
  # 20 year and 2 year peak flow 
  peak_flows <- calc_annual_highflows(data = hydro_clean,
                                      dates = date, 
                                      values = modelled,
                                      start_year = start_year, 
                                      end_year = end_year, 
                                      ignore_missing = TRUE)
  # reformat
  peak_flows <- dplyr::select(peak_flows, Year, Value = Max_1_Day)
  peak_flows <- dplyr::mutate(peak_flows, Measure = "1-Day")
  
  # compute freq analysis ** use_max = TRUE
  peak_freq_analysis <- compute_frequency_analysis(data = peak_flows,
                                                   events = Year,
                                                   values = Value,
                                                   measures = Measure, 
                                                   use_max = TRUE)
  peak_quant <- peak_freq_analysis$Freq_Fitted_Quantiles
  
  peak_1Q20 <- peak_quant$`1-Day`[2]
  peak_1Q5 <- peak_quant$`1-Day`[4]
  peak_1Q2 <- peak_quant$`1-Day`[5]
  
  ###################################################
  # calculate mean annual flow
  maf <- calc_longterm_mean(data = hydro_clean,
                            dates = date,
                            values = modelled, 
                            start_year = start_year,
                            end_year = end_year)
  
  ###################################################
  # calculate mean annual aug/sept flow 
  masf <- calc_longterm_mean(data = hydro_clean,
                             dates = date,
                             values = modelled,
                             start_year = start_year,
                             end_year = end_year,
                             months = 8:9)
  
  ####################################################
  # calculate average summer min flow DOY
  ann_extremes <- calc_annual_extremes(data = hydro_clean,
                                       dates = date,
                                       values = modelled,
                                       start_year = start_year,
                                       months = c(6,7,8,9),
                                       end_year = end_year, 
                                       ignore_missing = TRUE)
  
  av_min_doy <- round(mean(ann_extremes$Min_1_Day_DoY), 0)
  
  ####################################################
  # calculate average summer min flow DOY
  ann_extremes_full_yr <- calc_annual_extremes(data = hydro_clean,
                                               dates = date,
                                               values = modelled,
                                               start_year = start_year,
                                               end_year = end_year, 
                                               ignore_missing = TRUE)
  
  
  av_peak_doy <- round(mean(ann_extremes_full_yr$Max_1_Day_DoY), 0)
  
  metrics <- c("7Q10", "1Q20", "1Q5", "1Q2", "mean_annual_flow", "mean_aug_sep_flow", "av_peak_DOY", 
               "av_min_DOY")
  values <- c(round(sevQ10,2), 
              round(peak_1Q20,2),
              round(peak_1Q5, 2),
              round(peak_1Q2,2), 
              round(maf,2),
              round(masf,2),
              round(av_peak_doy,2),
              round(av_min_doy,2))
  
  results <- data.frame(metrics, values)
  names(results) <- c("metrics", "values")
  results$run <- run_name
  
  return(results)
}
