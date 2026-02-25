# ------------------------------------------------------------------------------
# metric difference figures for paper 2
#
# 
# 
# February 25, 2026
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
# generate metrics and save 
# ------------------------------------------------------------------------------ 

results <- lapply(hydrograph_list, 
                  calc_hydrograph_stats) %>%
  do.call(rbind, .)

write.csv(results, 
          file.path(outpath, "flow_metrics_feb25.csv"),
          row.names = FALSE)

results <- read.csv(file.path(outpath, "flow_metrics_feb24.csv"))

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

df <- df[df$scenario != "high_40", ]
df <- df[df$scenario != "low_40", ]

# ------------------------------------------------------------------------------
# make difference bar plot 
# ------------------------------------------------------------------------------

ggplot(df[df$metrics %in% c("1Q2", "1Q20", "mean_annual_flow"), ], aes(x = factor(yrs), y = value_diff)) +
  geom_bar(stat = "identity", colour = "black", fill = "grey") +
  facet_grid(metrics ~ scenario, scales = "free_y") +
  labs(
    x = "Years",
    y = "Value difference"
  ) +
  theme_bw()

ggplot(df[df$metrics %in% c("1Q2", "1Q20", "mean_annual_flow"), ], aes(x = factor(yrs), y = scenario, size = value_diff, colour = value_diff)) +
  geom_point(alpha = 0.6) +
  facet_grid(~metrics, scales = "free_y") +
  labs(
    x = "Years",
    y = "Value difference"
  ) +
  scale_colour_viridis(option = "plasma") + 
  theme_bw()


df_q20 <- df[df$metrics == "1Q20", ]
df_q2 <- df[df$metrics == "1Q2", ]
maf <- df[df$metrics == "mean_annual_flow", ]
df_7q10 <- df[df$metrics == "7Q10", ]
masf <- df[df$metrics == "mean_aug_sep_flow", ]
av_min_doy <- df[df$metrics == "av_min_DOY", ]
av_max_doy <- df[df$metrics == "av_peak_DOY", ]

q20_plot <- metric_time_series_plot(df_q20, "1Q20", fig_path) 
q2_plot <- metric_time_series_plot(df_q2, "1Q2",  fig_path)
maf_plot <- metric_time_series_plot(maf, "Mean Annual Flow", fig_path)
q710_plot <- metric_time_series_plot(df_7q10, "7Q10", fig_path)
masf_plot <- metric_time_series_plot(masf, "Mean Aug Sept Flow", fig_path)
av_min_doy_plot <- metric_time_series_plot(av_min_doy, "Average Date of Minimum Flow", fig_path)
av_peak_doy_plot <- metric_time_series_plot(av_max_doy, "Average Date of Peak Flow", fig_path)

metric_time_series_plot <- function(df_q20, metric_name, fig_path){
  
  # STEP 1 — Extract the numeric value & group (high/low)
  scenarios <- df_q20 %>%
    distinct(scenario) %>%
    mutate(
      scenario_value = parse_number(scenario),
      group = ifelse(grepl("high", scenario), "high", "low")
    )
  
  # STEP 2 — Build color palettes
  high_palette <- colorRampPalette(c("#a1d99b", "#006d2c"))   # light → dark green
  low_palette  <- colorRampPalette(c("#9ecae1", "#08519c"))   # light → dark blue
  
  # STEP 3 — Assign ONE color per scenario
  scenarios <- scenarios %>%
    group_by(group) %>%
    arrange(scenario_value) %>%
    mutate(
      scenario_color = case_when(
        group == "high" ~ high_palette(n())[row_number()],
        group == "low"  ~ low_palette(n())[row_number()]
      )
    ) %>%
    ungroup()
  
  # STEP 4 — Create named vector for scale_color_manual
  scenario_colors <- setNames(scenarios$scenario_color, scenarios$scenario)
  
  # STEP 5 — Join colors back to main df
  df2 <- df_q20 %>% left_join(scenarios, by = "scenario")
  
  legend_labels <- df2 %>%
    distinct(scenario) %>%
    mutate(
      scen_num  = parse_number(scenario),
      scen_type = ifelse(grepl("^high", scenario), "High", "Low"),
      pretty_label = paste0(scen_type, " (", scen_num, "%)")
    ) %>%
    { setNames(.$pretty_label, .$scenario) }
  
  q20 <- ggplot(df2, aes(x = yrs, y = values_run, color = scenario, group = scenario)) +
    
    geom_hline(aes(yintercept = values_af),
               linetype = "dotted", color = "black", linewidth = 0.8) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    
    scale_color_manual(values = scenario_colors,
                       labels = legend_labels) +
    
    theme_minimal(base_size = 14) +
    labs(
      x = "Years after Disturbance",
      y = expression(Streamflow~(m^3/s)),
      title = metric_name,
      color = "Simulation"
    )
  
  # ggsave(q20,
  #        filename = file.path(fig_path, paste0(metric_name, ".png")),
  #        units = "in",
  #        width = 8,
  #        height = 8)
}
# ---- 2. PLOT ----

p <- ggplot(df_q2, aes(x = scenario, y = yrs, size = value_diff)) +
  geom_point(alpha = 0.7) +
  #  facet_wrap(~ metrics, scales = "free") +
  scale_size_continuous(name = "Value magnitude") +
  labs(
    x = "Scenario",
    y = "Years",
    title = "Scenario × Years × Metric",
    subtitle = "Dot size = magnitude of metric value"
  ) +
  theme_bw()


# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------ 

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
  
  metrics <- c("7Q10", "1Q20", "1Q2", "mean_annual_flow", "mean_aug_sep_flow", "av_peak_DOY", 
               "av_min_DOY")
  values <- c(round(sevQ10,2), 
              round(peak_1Q20,2),
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
