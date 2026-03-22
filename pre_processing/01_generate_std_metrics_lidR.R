# ------------------------------------------------------------------------------
# LiDAR processing scripts 
# 
# date: February 9, 2026
# updated: March 20, 2026 to use PC that was normalized by 1m DEM 
# author: Brianne Boufford 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# libraries 
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(lubridate)
library(units)
library(smoothr)
library(lidR)
library(future)
library(stringr)
library(raster)
library(sp)

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------

# main project path 
project_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek")

output_folder <- file.path(project_path, "data", "lidar", "output_2026")

nlas_path <- file.path(output_folder, "normalized")

# read in las catolog 
nctg <- lidR::readALSLAScatalog(nlas_path, progress=TRUE,  filter = "-drop_z_below 0")

library(future)

plan(multisession, workers = 10)

metrics <- pixel_metrics(nctg, .stdmetrics, 20) # calculate standard metrics

terra::writeRaster(metrics, file.path(project_path, "data", "products", "std_metrics", "std_metrics_mar20.tif"),
                   overwrite=TRUE)
