# ------------------------------------------------------------------------------
# assess modelled temperature and precipitation 
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
weather_run_path <- file.path(".", "raven-runs", "Baseline3", "Runs", "Trapping_weather_baseline_feb23")
fig_path <- file.path("data", "figs", "baseline_model_calibration")

# list of temp 
temp_list <- list.files(weather_run_path, 
                        pattern = "TEMP",
                        recursive = TRUE, 
                        full.name = TRUE)

# list of precip
precip_list <- list.files(weather_run_path,
                          pattern = "PRECIP",
                          recursive = TRUE,
                          full.name = TRUE)

# ------------------------------------------------------------------------------
# read data 
# ------------------------------------------------------------------------------

clim <- read.csv(weather_file_path)
stations <- unique(clim$STATION_NAME)

# read in climate data 
clim_sf <- st_read(file.path(weather_data_path, "clim_station.shp"))

# subset the climate data to 1980 to now 
all_stat <- clim[clim$LOCAL_YEAR >= 1980, ]

# names of all weather stations
hru_names_ordered <- c("BIG WHITE", "BIG WHITE MTN LODGE","MCCULLOCH", "BEAVERDELL", "BEAVERDELL NORTH", "CARMI", "CHRISTIAN VALLEY")

# ------------------------------------------------------------------------------
# analyze temperature and precip
# ------------------------------------------------------------------------------

# daily average precip 
raven_precip <- read.csv(precip_list, skip = 1) %>% # from raven manual: rain/snow precipitation rate over time step /data interval [mm/d]
  select(-c("time", "X"))
raven_precip$day <- raven_precip$day %>% as.Date()

# fix column names on Raven files
names(raven_precip)[2:length(names(raven_precip))] <- paste0(hru_names_ordered) 
names(raven_precip)[1] <- "LOCAL_DATE"

# make year_month col
raven_precip$year_month <- paste0(year(raven_precip$LOCAL_DATE), "-", month(raven_precip$LOCAL_DATE))

# convert to long format
raven_precip_long <- pivot_longer(raven_precip, -c(LOCAL_DATE, year_month), values_to = "precip_raven", names_to = "STATION_NAME") %>%
  as.data.frame()

# monthyl summaries 
raven_monthly <- raven_precip_long %>%
  group_by(year_month, STATION_NAME) %>%
  dplyr::summarize(sum_precip = ifelse(all(is.na(precip_raven)), NA, sum(precip_raven, na.rm=TRUE))) %>%
  as.data.frame()
names(raven_monthly) <- c("year_month", "STATION_NAME", "sum_precip_raven")

# daily max temp average 
raven_max_temp <- read.csv(temp_list, skip = 1) %>%# maximum air temperature over day (0:00-24:00)[ ◦C]
  select(-c("time", "X"))
raven_max_temp$day <- raven_max_temp$day %>% as.Date()

# fix column names and convert to long format 
names(raven_max_temp)[2:length(names(raven_max_temp))] <- paste0(hru_names_ordered) 
names(raven_max_temp)[1] <- "LOCAL_DATE"
raven_temp_long <- pivot_longer(raven_max_temp, -c(LOCAL_DATE), values_to = "max_temp_raven", names_to = "STATION_NAME") %>%
  as.data.frame()

# --------------------------------------------------------------------------
# MEASURED WEATHER DATA 
# --------------------------------------------------------------------------

meas_precip <- all_stat[, c("LOCAL_DATE", "TOTAL_PRECIPITATION", "STATION_NAME")] # total rainfall and the SWE of total snowfall in mm during climatological day 
meas_precip$LOCAL_DATE <- meas_precip$LOCAL_DATE %>% as.Date()

# make year_month col
meas_precip$year_month <- paste0(year(meas_precip$LOCAL_DATE), "-", month(meas_precip$LOCAL_DATE))

# monthyl summaries 
meas_monthly_precip <- meas_precip %>%
  group_by(year_month, STATION_NAME) %>%
  dplyr::summarize(sum_value = ifelse(all(is.na(TOTAL_PRECIPITATION)), NA, sum(TOTAL_PRECIPITATION, na.rm = TRUE))) %>%
  as.data.frame()
names(meas_monthly_precip) <- c("year_month", "STATION_NAME", "sum_precip_meas")


meas_max_temp <- all_stat[, c("LOCAL_DATE", "MAX_TEMPERATURE", "STATION_NAME")]  # for climatological day 6:01AM to 6:00 AM but should be OK bc daily max is in the afternoon
meas_max_temp$LOCAL_DATE <- meas_max_temp$LOCAL_DATE %>% as.Date()


# measured data is already in long format 

# JOIN # -----------------------------------------------------------------------

# precip
precip_df <- merge(raven_monthly, meas_monthly_precip, by=c("year_month", "STATION_NAME")) %>%
  na.omit()

# max temp
temp_df <- merge(raven_temp_long, meas_max_temp, by=c("LOCAL_DATE", "STATION_NAME")) %>% 
  na.omit()

# add in month column 
precip_df$month_int <- lapply(precip_df$year_month, FUN = function(x){
  month_i <- stringr::str_split(x, pattern = "-")[[1]][2]
  return(month_i)
}) %>%
  do.call(rbind,.)

temp_df$month_int <- lubridate::month(temp_df$LOCAL_DATE)

r2_precip_all <- hydroGOF::R2(obs=precip_df$sum_precip_meas, sim=precip_df$sum_precip_raven)
r2_temp_all <- hydroGOF::R2(obs=temp_df$MAX_TEMPERATURE, sim=temp_df$max_temp_raven)

precip_df$season <- "winter"
precip_df$season[precip_df$month_int %in% c(4,5,6,7,8,9)] <- "summer"

temp_df$season <- "winter"
temp_df$season[temp_df$month_int %in% c(4,5,6,7,8,9)] <- "summer"

precip_station_season_summary <- precip_df %>%
  group_by(STATION_NAME, season) %>%
  summarize(N = length(month_int),
            r2 = hydroGOF::R2(sum_precip_raven, sum_precip_meas),
            pbias = hydroGOF::pbias(sum_precip_raven, sum_precip_meas), 
            data = "precip")

precip_station_summary <- precip_df %>%
  group_by(STATION_NAME) %>%
  summarize(N = length(month_int),
            r2 = hydroGOF::R2(sum_precip_raven, sum_precip_meas),
            pbias = hydroGOF::pbias(sum_precip_raven, sum_precip_meas), 
            data = "precip",
            season = "full year")

temp_station_season_summary <- temp_df %>%
  group_by(STATION_NAME, season) %>%
  summarize(N = length(month_int),
            r2 = hydroGOF::R2(max_temp_raven, MAX_TEMPERATURE),
            pbias= hydroGOF::pbias(max_temp_raven, MAX_TEMPERATURE), 
            data = "temp")

temp_station_summary <- temp_df %>%
  group_by(STATION_NAME) %>%
  summarize(N = length(month_int),
            r2 = hydroGOF::R2(max_temp_raven, MAX_TEMPERATURE),
            pbias = hydroGOF::pbias(max_temp_raven, MAX_TEMPERATURE),
            season = "full year", 
            data = "temp")

print(paste0("r2 temp all: ", round(r2_temp_all,2)))
print(paste0("r2 precip all: ", round(r2_precip_all,2)))

#temp_df <- merge(temp_df, temp_station_summary, by="STATION_NAME")
#precip_df <- merge(precip_df, precip_station_summary, by="STATION_NAME")

results <- rbind(precip_station_season_summary, precip_station_summary) %>%
  rbind(., temp_station_season_summary) %>%
  rbind(., temp_station_summary)

# plot 
temp_plot <- ggplot(data=temp_df, aes(x=MAX_TEMPERATURE, y=max_temp_raven, color=season)) + 
  geom_point(size = 0.2) + 
  geom_line(aes(x=MAX_TEMPERATURE, y=MAX_TEMPERATURE), color = "red", linewidth=0.75) + 
  theme_minimal() +
  facet_wrap(~ STATION_NAME) +
  labs(x = "Measured Max Temp (deg C)", y = "Raven Max Temp (deg C)", color = "Season") +
  scale_color_manual(values = c("winter" = "darkblue", "summer" = "orange")) + 
  geom_text(data = temp_station_summary, aes(x = -10, y = 30, label = paste0("R2: ", round(r2,2))), 
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 4)

precip_plot <- ggplot(data=precip_df, aes(x=sum_precip_meas, y=sum_precip_raven, color=season)) + 
  geom_point(size = 1) + 
  geom_line(aes(x=sum_precip_meas, y=sum_precip_meas), color = "red", linewidth=0.75) + 
  theme_minimal() +
  facet_wrap(~ STATION_NAME) +
  labs(x = "Measured Monthly Precip (mm)", y = "Raven Monthly Precip (mm)", color = "Season") +
  scale_color_manual(values = c("winter" = "darkblue", "summer" = "orange")) + 
  geom_text(data = precip_station_summary, aes(x = 75, y = 225, label = paste0("R2: ", round(r2,2))), 
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 4)

ggsave(temp_plot,
       filename = file.path(fig_path, paste0("temp_compare_feb23",".png")),
       units = "in", 
       height = 4,
       width = 6)

ggsave(precip_plot,
       filename = file.path(fig_path, paste0("precip_compare_feb23", ".png")),
       units = "in", 
       height = 4,
       width = 6)

# average annual behaviour plot 
temp_summary <- temp_df %>%
  group_by(STATION_NAME, month_int) %>%
  summarise(
    Modelled = mean(max_temp_raven, na.rm = TRUE),
    Measured = mean(MAX_TEMPERATURE, na.rm = TRUE),
    sd_Modelled = sd(max_temp_raven, na.rm = TRUE),
    sd_Measured = sd(MAX_TEMPERATURE, na.rm = TRUE),
    .groups = "drop"
  ) 

temp_long <- temp_summary %>%
  pivot_longer(cols = c(Modelled, Measured), names_to = "Type", values_to = "Mean") %>%
  mutate(
    SD = ifelse(Type == "Modelled", sd_Modelled, sd_Measured),
    ymin = Mean - SD,
    ymax = Mean + SD
  )

temp_long$month_name <- month.abb[temp_long$month_int]

av_temp_plot <- ggplot(temp_long, aes(x = month_int, y = Mean, color = Type, fill = Type, group = Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  facet_wrap(~ STATION_NAME) +
  labs(x = "Month (Jan-Dec)", y = "Temperature (deg C)", color = "", fill = "") +
  theme_minimal() + 
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  scale_color_manual(values = c("Modelled" = "darkblue", "Measured" = "orange")) +
  scale_fill_manual(values = c("Modelled" = "darkblue", "Measured" = "orange")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1)
    #legend.position = c(0.86, 0.09),        # x, y position in normalized plot coordinates
    #legend.justification = c("right", "bottom"), 
  )

ggsave(av_temp_plot,
       filename = file.path(fig_path, paste0("av_temp_compare_feb23", ".png")),
       units = "in", 
       height = 4,
       width = 6)

precip_df$month_int <- as.numeric(precip_df$month_int)

# average annual behaviour plot 
precip_summary <- precip_df %>%
  group_by(STATION_NAME, month_int) %>%
  summarise(
    Modelled = mean(sum_precip_raven, na.rm = TRUE),
    Measured = mean(sum_precip_meas, na.rm = TRUE),
    sd_Modelled = sd(sum_precip_raven, na.rm = TRUE),
    sd_Measured = sd(sum_precip_meas, na.rm = TRUE),
    .groups = "drop"
  ) 

precip_long <- precip_summary %>%
  pivot_longer(cols = c(Modelled, Measured), names_to = "Type", values_to = "Mean") %>%
  mutate(
    SD = ifelse(Type == "Modelled", sd_Modelled, sd_Measured),
    ymin = Mean - SD,
    ymax = Mean + SD
  )

precip_long$month_name <- month.abb[as.numeric(precip_long$month_int)]

av_precip_plot <- ggplot(precip_long, aes(x = factor(month_int), y = Mean, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), color = "black", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  facet_wrap(~ STATION_NAME) +
  labs(x = "Month (Jan–Dec)", y = "Monthly Precip (mm)", fill = "") +
  theme_minimal() +
  scale_x_discrete(
    breaks = c("1", "4", "7", "10"),
    labels = c("Jan", "Apr", "Jul", "Oct")
  ) +
  scale_fill_manual(values = c("Modelled" = "darkblue", "Measured" = "orange"))# +
# theme(
#   axis.text.x = element_text(angle = 0, hjust = 0.5),  # center the label under month
#   legend.position = c(0.86, 0.09),
#   legend.justification = c("right", "bottom")
# )

ggsave(av_precip_plot,
       filename = file.path(fig_path, paste0("av_precip_compare_feb23", ".png")),
       units = "in", 
       height = 4,
       width = 6)

av_together <- grid.arrange(av_temp_plot, av_precip_plot, ncol = 2, widths=c(2,3))
ggsave(av_together,
       filename = file.path(fig_path, paste0("av_precip_temp_compare_feb23", ".png")),
       units = "in", 
       height = 4,
       width = 12)

  

