# ------------------------------------------------------------------------------
# LiDAR processing scripts 
# LAI spatial products 
# 
# date: August 28th, 2024
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

nctg <- lidR::readALSLAScatalog(nlas_path, progress=TRUE,  filter = "-drop_z_below 0")

# from Gonzo -- LAI metrics 

metrics_lad <- function(z, zmin=NA, dz = 1, k = 0.5, z0 = 2) {
  if (!is.na(zmin)) z <- z[z>zmin]
  lad_max <- lad_mean <- lad_cv <- lad_min <- lai <- NA_real_
  if(length(z) > 2) {
    ladprofile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)
    lad_max <- with(ladprofile, max(lad, na.rm = TRUE))
    lad_mean <- with(ladprofile, mean(lad, na.rm = TRUE))
    lad_cv <- with(ladprofile, sd(lad, na.rm=TRUE)/mean(lad, na.rm = TRUE))
    lad_min <- with(ladprofile, min(lad, na.rm = TRUE))
    lai <- with(ladprofile, sum(lad, na.rm = TRUE))
  }
  lad_metrics <- list(lad_max = lad_max,
                      lad_mean = lad_mean,
                      lad_cv = lad_cv,
                      lad_min = lad_min,
                      lai = lai)
  return(lad_metrics)
}

library(future)
plan(multisession, workers = 10)
LAI <- grid_metrics(nctg , ~metrics_lad(Z), res=20)
LAI_rast <- rast(LAI)
terra::writeRaster(LAI_rast, file.path(project_path, "data", "products", "LAI", "lad_metrics_norm_mar20.tif"),
                   overwrite=TRUE)

# subset to just the LAI layer and then remove some errounous values near big white
lai <- LAI_rast$lai 
lai[lai > 10] <- NA

# write out just LAI layer
terra::writeRaster(lai, file.path(project_path, "data", "products", "LAI", "lai_clean_mar20.tif"),
                   overwrite=TRUE)
