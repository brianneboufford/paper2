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

# ------------------------------------------------------------------------------
# generate metrics and save 
# ------------------------------------------------------------------------------ 

results <- lapply(hydrograph_list, 
                  calc_hydrograph_stats) %>%
  do.call(rbind, .)

write.csv(results, 
          file.path(outpath, "flow_metrics_feb25.csv"),
          row.names = FALSE)

results <- read.csv(file.path(outpath, "flow_metrics_mar24.csv"))

# --------------------------------
# results described in paper 
# --------------------------------
# done in excel sheet for ease
 
# ------------------------------------------------------------------------------
# make date of peak flow and peak flow box plot 
# ------------------------------------------------------------------------------

hydrograph_af <- read.csv(hydrograph_list[1])
hydro_af <- prep_model_data_af(hydrograph_af, "All Forested")

sims_list <- hydrograph_list[-1]
sims_list <- sims_list[!grepl("_40_", sims_list)]

sims_data <- lapply(sims_list, 
                    prep_model_data_sims) %>%
  do.call(rbind,.)

sims_data <- sims_data[sims_data$Type == "Simulated", ]

sims_data$year <- year(as.Date(sims_data$date))

sims_data <- sims_data %>% 
  subset(!year %in% c(1980, 2023))

sims_data_ann <- sims_data %>% 
  group_by(year, Site, recovery) %>%
  summarize(
    annual_peak_flow = max(Value),
    date_peak_flow = yday(as.Date(date))[which.max(Value)] 
  )

sims_data_ann <- sims_data_ann %>% 
  mutate(recovery =  
         factor(recovery, 
                levels = c("0yrs", "5yrs", "10yrs", "15yrs",
                           "20yrs", "25yrs", "30yrs", "35yrs", "40yrs",
                           "45yrs", "50yrs"), 
                labels = c("0-5", "6-10", "11-15", "16-20",
                           "21-25", "26-30", "31-35", "36-40", "41-45",
                           "46-50", "AF")))

peak_flow_box <- ggplot(data = sims_data_ann, aes(x = Site, y = annual_peak_flow, fill = recovery)) +
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
       y = expression("Peak Streamflow ("*m^3/s*")"), 
       legend = "Recovery") + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "right")

# dont love these boxplots 

date_peak_flow_box <- ggplot(data = sims_data_ann, aes(x = Site, y = date_peak_flow, fill = recovery)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) + 
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
  coord_flip() + 
  labs(x = "Simulation", 
       y = expression("DoY Peak Streamflow"), 
       legend = "Recovery") + 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "right")


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

df <- df %>% 
  mutate(scenario = factor(scenario, 
                       levels = c("high_10", "high_15", "high_20", "high_30", 
                                  "low_10", "low_15", "low_20", "low_30"), 
                       labels = c(expression("High z":10*'%'),
                                  expression("High z":15*'%'), 
                                  expression("High z":20*'%'), 
                                  expression("High z":30*'%'), 
                                  expression("Low z":10*'%'), 
                                  expression("Low z":15*'%'), 
                                  expression("Low z":20*'%'), 
                                  expression("Low z":30*'%'))))

df_flood <- df[df$metrics %in% c("1Q2", "1Q5", "1Q20", "mean_annual_flow"), ] %>% 
  mutate(metrics = factor(metrics, 
                           levels = c("1Q2","1Q5", "1Q20", "mean_annual_flow"), 
                           labels = c(expression(1*Q*2),
                                      expression(1*Q*5),
                                      expression(1*Q*20), 
                                      expression(Mean~Annual~Q))))

df_flood <- df_flood %>% 
  mutate(yrs = factor(yrs, 
                      levels = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
                      labels = c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35", 
                                 "36-40", "41-45", "46-50", "AF")))

df_flood <- df_flood[!grepl("AF", df_flood$yrs), ]
df_flood$perc_diff <- round(df_flood$value_diff/df_flood$values_af*100, 2)

# ------------------------------------------------------------------------------
# make difference bar plot 
# ------------------------------------------------------------------------------

Qbar <- ggplot(df_flood, aes(x = factor(yrs), y = value_diff/values_af*100, fill = yrs)) +
  geom_bar(stat = "identity", colour = "grey3") +
  scale_fill_manual(values = c("0-5" = "#8e0152", 
                               "6-10" = "#c51b7d", 
                               "11-15" = "#de77ae", 
                               "16-20" = "#f1b6da", 
                               "21-25" = "#fde0ef", 
                               "26-30" = "#f7f7f7", 
                               "31-35" = "#e6f5d0", 
                               "36-40" = "#b8e186", 
                               "41-45" = "#7fbc41", 
                               "46-50" = "#4d9221")) + 
  facet_grid(metrics ~ scenario, scales = "free_y", 
             labeller = as_labeller(label_parsed)) +
  labs(
    x = "",
    y = expression("Percent Difference from All Forested (%)"),
    fill = "Recovery"
  ) +
  theme_bw(14) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_rect(NA),
        axis.title = element_text(size = 12),
        legend.position = "bottom", 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size =12)) + 
  guides(fill = guide_legend(nrow = 1))

ggsave(Qbar, 
       filename = file.path(fig_path, "Qmetric_barplot_mar24.png"),
       units = "in",
       width = 10,
       height = 6, 
       dpi = 300)

# ------------------------------------------------------------------------------
# make line and point plot April 9, 2026
# ------------------------------------------------------------------------------

df_flood$elev <-str_split(df_flood$run, pattern = "_", simplify = TRUE)[, 1]
df_flood$dist_size <-str_split(df_flood$run, pattern = "_", simplify = TRUE)[, 2]
df_flood$yrs_str <-str_split(df_flood$run, pattern = "_", simplify = TRUE)[, 3]

df_flood$elev <- dplyr::recode(df_flood$elev,
                        "high" = "High",
                        "low"  = "Low")

Qline <- ggplot(df_flood, 
                aes(x = yrs, 
                    y = value_diff/values_af*100, 
                    group = interaction(scenario, dist_size, elev),
                    color = elev,
                    shape = dist_size)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~ metrics, 
             labeller = as_labeller(label_parsed)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_shape_manual(values = c(
    "10" = 1,   # open circle
    "15" = 2,   # open triangle
    "20" = 0,   # open square
    "30" = 8    # star
  )) +
  scale_color_manual(values = c(
    "High" = "#c51b7d",
    "Low"  = "#4d9221"
  )) +
  scale_y_continuous(breaks = seq(0, 25, 5)) + 
  labs(
    x = "Recovery Stage",
    y = expression("Difference from All Forested (%)"),
    color = "Elevation Region",
    shape = "Disturbance Scenario (%)"
  ) +
  theme_bw(14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = NA),
    axis.title = element_text(size = 12),
    legend.position = "bottom", 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = NA, linewidth = 1))
  )

ggsave(Qline, 
       filename = file.path(fig_path, "Qmetric_lineplot_april9.png"),
       units = "in",
       width = 8,
       height = 10, 
       dpi = 300)

# ------------------------------------------------------------------------------
# old plots 
# ------------------------------------------------------------------------------
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
      scenario_value = readr::parse_number(scenario),
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

# function to clean hydrograph data and add column with simulation name 
prep_model_data_sims <- function(df_path){
  
  sim_name <- dirname(df_path) %>% basename() %>% str_remove("Trappingp50_") %>%
    str_remove("_(\\d+?)yrs") %>% str_replace("_", " ") %>% paste0(., "%")
  
  recovery <- dirname(df_path) %>% basename() %>% str_remove("Trappingp50_") %>%
    str_remove("_(\\d+?)_") %>% str_replace("high", "") %>% str_replace("low"
    , "")
  
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
  df$recovery <- recovery
  return(df)
} 


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

