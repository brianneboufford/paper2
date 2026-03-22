# ------------------------------------------------------------------------------
# LiDAR processing scripts 
# - build out DEM, CHM, LAI spatial products 
# 
# date: July 30th, 2024
# added pabove2m November 15, 2025
# copied and redone March 19, 2026
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

las_path <- file.path("F:", "Trapper_LAS", "Trapper_LAS")

las_files <- list.files(las_path, pattern = "\\.las$", full.names = TRUE)

output_folder <- file.path(project_path, "data", "lidar", "output_2026")

if(!file.exists(output_folder)){
  dir.create(output_folder, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# set up parallel processing 
# ------------------------------------------------------------------------------

# run this later ** 
plan("multisession", workers = 10L)

# ------------------------------------------------------------------------------
# clean temp dir 
# ------------------------------------------------------------------------------

# Set Raster Options
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

rasterOptions(tmptime = 72, tmpdir = path.expand("~/tmp"), progress = "")

# ------------------------------------------------------------------------------
# make DTM 
# ------------------------------------------------------------------------------

# read in las catolog 
tc_ctg <- lidR::readLAScatalog(las_path, progress = TRUE, filter = "-drop_z_below 0")

# complete las check
las_check(tc_ctg)

# show estimated progress 
opt_progress(tc_ctg) <- TRUE

# output folder 
opt_output_files(tc_ctg) <- paste0(output_folder, "/{ORIGINALFILENAME}")

# default chunk size is 0
opt_chunk_size(tc_ctg) <- 0

cell_size <- 1

# develop dem 
dem <- rasterize_terrain(tc_ctg, 
                         size = cell_size, 
                         algorithm = tin())

writeRaster(dem, 
            filename = file.path(output_folder, "dem_1m.tif"), 
            overwrite = TRUE)

# normalization is not done 
# ground classification is done 
# duplicates are not dealt with yet 

# ------------------------------------------------------------------------------
# generate footprint of TC 
# ------------------------------------------------------------------------------
footprint <- as.spatial(tc_ctg) %>%
  st_as_sf()

id_col <- footprint$filename 

# function for extracting unique Tile ID for las file 
get_id <- function(id){
  # ------------------------------------------------
  # inputs:
  # @id - full path or filename for single las file 
  #
  # outputs: 
  # unqiue lidar tile ID 
  # -------------------------------------------------
  
  # split filename string and pull out unique tile ID 
  id1 <- strsplit(id, "INTOK")[[1]][2]
  id2 <- strsplit(id1, "BCAlbers")[[1]][1]
  id3 <- str_replace(id2, "_$", "") %>%
    str_replace(., "^_", "")
  
  # unique ID 
  return(id3)
}

footprint$id <- lapply(footprint$filename, get_id) %>%
  cbind() %>%
  as.character()

footprint_outfolder <- file.path(output_folder, "footprint")

if (!dir.exists(footprint_outfolder)){
  dir.create(footprint_outfolder, recursive = TRUE)
}


st_write(footprint, 
         file.path(footprint_outfolder, "tc_footprint.shp"),
         delete_dsn = TRUE)

# ------------------------------------------------------------------------------
# nornmalize las files and build chm 
# ------------------------------------------------------------------------------

dem_path <- file.path(output_folder, "dem_1m.tif")

ctg <- readLAScatalog(las_path)

dem <- raster(dem_path)

output <- file.path(output_folder, "normalized")

if(!dir.exists(output)){
  dir.create(output)
}

opt_stop_early(ctg) <- FALSE
opt_output_files(ctg) <-paste0(output, "/{ORIGINALFILENAME}")
opt_laz_compression(ctg) <- TRUE
opt_progress(ctg) <- TRUE

ctg_norm <- normalize_height(ctg, tin(), dtm = dem)

ctg <- readLAScatalog(output, filter = 'drop_z_below 0 -drop_z_above 65')

output <- file.path(output_folder, "chm_tiles")

if(!dir.exists(output)){
  dir.create(output)
}

opt_stop_early(ctg) <- FALSE

opt_output_files(ctg) <-paste0(output, "/{ORIGINALFILENAME}")

chm <- rasterize_canopy(ctg, res = cell_size, pitfree(c(0,10,20), c(0,1.5)))

writeRaster(chm, 
            filename = file.path(output, "chm_1m.tif"),  
            overwrite = TRUE)

