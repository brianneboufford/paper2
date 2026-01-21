################################################################################
#
# compare biosim and rvt data 
# 
# November 12, 2025 
################################################################################

# libraries
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)
library(readr)
library(plotly)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# set wd 
setwd(file.path("C:","Users", "blbouf", "Sync", "Paper2"))

# paths
raven_run_path <- file.path(".", "raven-runs", "Baseline")
biosim_path <- file.path(raven_run_path, "biosim")
daymet_path <- file.path(raven_run_path, "daymet")

rvt_files <- list.files(biosim_path)

all_stat_data <- lapply(rvt_files, 
                        prep_rvt_df, 
                        daymet_path = daymet_path,
                        biosim_path = biosim_path) %>%
  do.call(rbind, .)

sub_d <- all_stat_data[all_stat_data$fileN == "File1" & all_stat_data$date < as.Date("1990-01-01"), ]
# make plot
ggplot(sub_d, aes(x = date, y = as.numeric(PRECIP), fill = data_src)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  facet_wrap(~ fileN, scales = "free_y") +
  labs(
    title = "Daily Precipitation Over Time",
    x = "Date",
    y = "Daily Precipitation (mm)",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

ggplotly()

ggplot(sub_d, 
       aes(x = date, y = as.numeric(TEMP_DAILY_MAX), color = data_src)) + 
  geom_line() + facet_wrap(~ fileN, scales = "free_y") + 
  labs( title = "Daily Maximum Temperature Over Time", x = "Date", 
        y = "Daily Max Temperature (°C)", color = "Data Source" ) +
  theme_minimal(base_size = 14) + 
  theme( panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(), 
         strip.background = element_rect(fill = "grey90", color = NA), 
         strip.text = element_text(face = "bold") )

ggplotly()

# ------------- 
# pivot data
# -------------
# Assuming your dataframe is named df
df_wide <- all_stat_data %>%
  select(fileN, date, data_src, TEMP_DAILY_MAX) %>%
  pivot_wider(
    names_from = data_src,
    values_from = TEMP_DAILY_MAX,
    names_prefix = "TEMP_DAILY_MAX_"
  )

r2_df <- df_wide %>%
  group_by(fileN) %>%
  summarise(
    r2 = cor(as.numeric(TEMP_DAILY_MAX_daymet), as.numeric(TEMP_DAILY_MAX_biosim), use = "complete.obs")^2,
    pbias = 100 * sum(as.numeric(TEMP_DAILY_MAX_daymet) - as.numeric(TEMP_DAILY_MAX_biosim), na.rm = TRUE) /
      sum(as.numeric(TEMP_DAILY_MAX_biosim), na.rm = TRUE)
  )

ggplot(df_wide, aes(x = TEMP_DAILY_MAX_daymet, y = TEMP_DAILY_MAX_biosim)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ fileN) +
  geom_text(
    data = r2_df,
    aes(x = -Inf, y = Inf, label = paste0("R² = ", round(r2, 2))),
    hjust = -0.1, vjust = 1.5,
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    title = "Comparison of Daily Max Temperature: Daymet vs Biosim",
    x = "Daymet TEMP_DAILY_MAX (°C)",
    y = "Biosim TEMP_DAILY_MAX (°C)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggplotly()

# function to prep daymet and biosim data for a single station
# provided by rvt_file
# returns joined weather station data
prep_rvt_df <- function(rvt_file, daymet_path, biosim_path){
  
  # read daymet data
  daymet <- read.csv(file.path(daymet_path, rvt_file),
                     skip = 4,
                     header = FALSE)
  names(daymet) <- c("TEMP_DAILY_MAX", "TEMP_DAILY_MIN", "PRECIP", "REL_HUMIDITY")
  
  # read biosim data
  biosim <- read.csv(file.path(biosim_path, rvt_file),
                     skip = 4, 
                     header = FALSE)
  names(biosim) <- c("TEMP_DAILY_MAX", "TEMP_DAILY_MIN", "PRECIP", "REL_HUMIDITY")
  
  # add date column to each starting from Jan 1, 1980 
  start_date <- as.Date("1980-01-01")
  daymet$date <- seq.Date(from = start_date, by = "day", length.out = nrow(daymet))
  biosim$date <- seq.Date(from = start_date, by = "day", length.out = nrow(biosim))
  
  # add file #
  daymet$fileN <- gsub(".rvt", "", rvt_file)
  biosim$fileN <- gsub(".rvt", "", rvt_file)
  
  # add data source
  daymet$data_src <- "daymet"
  biosim$data_src <- "biosim"
  
  # join
  stat_data <- rbind(daymet, biosim)
  
  return(stat_data)
}

