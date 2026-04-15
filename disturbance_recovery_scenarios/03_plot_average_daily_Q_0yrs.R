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
library(stringr)

setwd("C:/Users/blbouf/Sync/Paper2")
# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios", "Runs")
fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Mar22")
hru_run_path <- file.path(".", "raven-runs", "Baseline3", "Runs", "Trapping_HRU_baseline_Feb20") # was dec 10
hru_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios_Feb19", "Runs")
# list of hydrographs
hydrograph_list_all <- list.files(hru_run_path, pattern="Hydrographs.csv", 
                                  recursive = TRUE, 
                                  full.name = TRUE)

hydrograph_list_af <- hydrograph_list_all[1] %>% read.csv()

# subset to only recently disturbed runs (0yrs)
hydrograph_0yrs_list <- hydrograph_list_all[grepl(hydrograph_list_all, pattern = "_0yrs")]
hydrograph_0yrs_list <- hydrograph_0yrs_list[!grepl(hydrograph_0yrs_list, pattern = "_40_")]

# prep all forested data
hydroaf <- prep_model_data_af(hydrograph_list_af, "af") %>% rbind()

# prep all other data
hydro_data <- lapply(hydrograph_0yrs_list, prep_model_data_sims) %>% 
  do.call(rbind, .)

Qup = 0.9
Qlow = 0.1

# group data and calculate quantiles
sim_data_grouped <- hydro_data %>%
  group_by(Site, date = yday(as.Date(date)), Type) %>%
  mutate(date = yday(as.Date(date))) %>%
  group_by(Site, date, Type) %>%
  summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
            Qlow = quantile(Value, Qlow, na.rm = TRUE),
            Value = mean(Value, na.rm = TRUE)) %>%
  ungroup() 

af_data_grouped <- hydroaf %>%
  group_by(Site, date = yday(as.Date(date)), Type) %>%
  mutate(date = yday(as.Date(date))) %>%
  group_by(Site, date, Type) %>%
  summarise(Qhigh = quantile(Value, Qup, na.rm = TRUE),
            Qlow = quantile(Value, Qlow, na.rm = TRUE),
            Value = mean(Value, na.rm = TRUE)) %>%
  ungroup() 

af_data_grouped <- af_data_grouped[af_data_grouped$Type == "Simulated", ]
sim_data_grouped <- sim_data_grouped[sim_data_grouped$Type == "Simulated", ]

#all_data_grouped <- rbind(af_data_grouped, model_data_grouped)

af_data_grouped$Site[af_data_grouped$Site == "af"] <- "All Forested"
af_data_grouped <- af_data_grouped %>% select(-c("Site"))

sim_data_grouped <- sim_data_grouped %>% 
  mutate(Site = factor(Site, 
                       levels = c("high 10%", "high 15%", "high 20%", "high 30%", 
                                  "low 10%", "low 15%", "low 20%", "low 30%"), 
                       labels = c(expression("High Elevation":~10*'%'),
                                  expression("High Elevation":~15*'%'), 
                                  expression("High Elevation":~20*'%'), 
                                  expression("High Elevation":~30*'%'), 
                                  expression("Low Elevation":~10*'%'), 
                                  expression("Low Elevation":~15*'%'), 
                                  expression("Low Elevation":~20*'%'), 
                                  expression("Low Elevation":~30*'%'))))

av_plot <- ggplot(
  data = sim_data_grouped,
  aes(x = as.Date(strptime(date, format = "%j")), y = Value)
) +
  # Simulation
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh, fill = "Simulation"),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(aes(colour = "Simulation")) +
  
  facet_wrap(~Site, nrow = 2, labeller = as_labeller(label_parsed)) +
  
  # All forested
  geom_ribbon(
    data = af_data_grouped,
    aes(
      x = as.Date(strptime(date, format = "%j")),
      y = Value,
      ymin = Qlow,
      ymax = Qhigh,
      fill = "All Forested"
    ),
    alpha = 0.25,
    colour = NA,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = af_data_grouped,
    aes(
      x = as.Date(strptime(date, format = "%j")),
      y = Value,
      colour = "All Forested"
    ),
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    name = "",
    values = c(
      "Simulation"   = "#c51b7d",
      "All Forested" =  "#4d9221"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "Simulation"   = "#c51b7d",
      "All Forested" =  "#4d9221"
    )
  ) +
  scale_x_date(
    breaks = as.Date(paste0("2026-", c("01-01","04-01","07-01","10-01"))),
    labels = c("Jan", "Apr", "Jul", "Oct")
  ) + 
  labs(y = expression("Average Daily Streamflow (" * m^3/s * ")")) +
  theme_bw(14) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10), 
    strip.background = element_rect(fill = NA)
  )

av_plot

ggsave(av_plot,
       filename = file.path(fig_path, "av_daily_Q_0yrs_mar24.png"),
       units = "in",
       dpi = 300,
       width = 8, 
       height = 5)


# ---------------------------------------------------------------
# results for paper 2 Fig 7
#----------------------------------------------------------------
sdg <- sim_data_grouped %>% as.data.frame()
af <- af_data_grouped %>% as.data.frame()
af_max <- max(af$Value)

h30_max <- max(sdg$Value[sdg$Site == expression("High Elevation":~30 * "%")])

l30_max <- max(sdg$Value[sdg$Site == expression("Low Elevation":~30 * "%")])

h10_max <- max(sdg$Value[sdg$Site == expression("High Elevation":~10 * "%")])

l10_max <- max(sdg$Value[sdg$Site == expression("Low Elevation":~10 * "%")])

round(h30_max/af_max*100, 2)
round(l30_max/af_max*100, 2)
round(h10_max/af_max*100, 2)
round(l10_max/af_max*100, 2)
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
