# ------------------------------------------------------------------------------
# create percent above 2m layer from las files 
#
# date: november 13, 2025
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# libraries 
# ------------------------------------------------------------------------------

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
library(sp)
library(lidRmetrics)

#devtools::install_github("ptompalski/lidRmetrics")

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
setwd("C:/Users/blbouf/Sync/Paper2")

las_path <- file.path("F:", "Trapper_LAS", "Trapper_LAS")
norm_las_path <- file.path("F:", "Trapper_LAS", "normalized")

las_files <- list.files(las_path, pattern = "\\.las$", full.names = TRUE)

data_path <- file.path(".", "data")

output_folder <- file.path(data_path, "fCOVER_ALS")

if(!file.exists(output_folder)){
  dir.create(output_folder, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# set up parallel processing and clean temp dir
# ------------------------------------------------------------------------------
plan("multisession", workers = 5L)

# ------------------------------------------------------------------------------
# clean temp dir 
# ------------------------------------------------------------------------------

#Set Raster Options
if(!dir.exists(path.expand("~/tmp"))) {
  dir.create(path.expand("~/tmp"))
}

clean_temp <- function() {
  files <- list.files(path = "~/tmp", full.names = TRUE)
  for (fn in files) {
    if (file.exists(fn)) {
      file.remove(fn)
    }
  }
}

clean_temp()

raster::rasterOptions(tmptime = 72, tmpdir = path.expand("~/tmp"), progress = "")

# ------------------------------------------------------------------------------
# process pabove2
# ------------------------------------------------------------------------------
ctg <- readLAScatalog(norm_las_path, filter = 'drop_z_below 0 -drop_z_above 65')
#las <- readLAS(list.files(norm_las_path, full.name = TRUE)[1],filter = 'drop_z_below 0 -drop_z_above 65')

opt_stop_early(ctg) <- TRUE
opt_output_files(ctg) <-paste0(output_folder, "/{ORIGINALFILENAME}")
opt_laz_compression(ctg) <- TRUE
opt_progress(ctg) <- TRUE

metr_perab <- pixel_metrics(ctg, ~lidRmetrics::metrics_percabove(Z), res = 30)

writeRaster(metr_perab, 
            filename = file.path(output_folder, "percent_above_metrics.tiff"),
            overwrite =- TRUE)
