# ------------------------------------------------------------------------------
# script for downloading BioSim Daily weather data for Paper 2 model 
# 
# author: Brianne Boufford
# Date : September 16, 2025 
# ------------------------------------------------------------------------------
# intall BioSimClient just once else load package 
if (length(find.package("BioSimClient_R", quiet = TRUE)) > 0){
  devtools::install_github("https://github.com/RNCan/BioSimClient_R")
  library(BioSIM)
} else {
  library(BioSIM)
}
# packages 
library(dplyr)
library(devtools)
library(sf)

# Point to your Java 8 bin folder
java_bin <- "C:/Program Files/Java/jdk-1.8/bin"

# Update PATH inside R so BioSIM can find java
Sys.setenv(PATH = paste(java_bin, Sys.getenv("PATH"), sep = ";"))

# Confirm R can find java just by name
system2("java", args = "-version", stdout = TRUE, stderr = TRUE)

# Now try BioSIM again
library(BioSIM)
getModelList()

# paths 
data_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data")
weather_points_path <- file.path(data_path, "BioSim", "src", "weather_points")

# load weather data
weather_points_df <- read.csv(file.path(weather_points_path, "weather_points.csv"))
weather_points <- st_as_sf(weather_points_df, coords = c("Lon", "Lat"), crs=4326)

# write file once 
if (!file.exists(file.path(weather_points_path, "weather_points.shp"))){
  st_write(weather_points, 
           file.path(weather_points_path, "weather_points.shp"), 
           append = FALSE)
}
# ------------------------------------------------------------------------------
# get data for first point 

ids <- weather_points_df$ID
lats <- weather_points_df$Lat %>% as.numeric()
lons <- weather_points_df$Lon %>% as.numeric()
elevs <- weather_points_df$Elev %>% as.numeric()

Weather <- generateWeather("Climatic_Daily", fromYr = 1980, toYr = 2023, id = ids, 
                           latDeg = lats, longDeg = lons, elevM = elevs)
Weather <- Weather$Climatic_Daily

# function to subset weather data for single RVT
prep_rvt_data <- function(ID_int, Weather){
  
  weather_subset <- Weather[Weather$KeyID == paste0("P", ID_int), ] %>% 
    select(c("Tmax", "Tmin", "Prcp", "RelH"))
  weather_subset$RelH <- weather_subset$RelH/100
  
  write.csv(weather_subset,
            file.path(weather_points_path, "..", "..", "outputs", "CSVs", paste0("File", ID_int, ".csv")),
            row.names = FALSE)
}

ID_int_list <- c(1,2,3,4,5,6,7,8,9)
lapply(ID_int_list, 
       prep_rvt_data,
       Weather = Weather)

# need to divide RelH by 100 
# plot to compare to daymet
  
