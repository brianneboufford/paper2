# ------------------------------------------------------------------------------
# script to clip fCOVER to study area and join into a single raster stack 
# 
# date: november 12, 2025
# author: Brianne Boufford 
# modified from scale_lai_stack_rasters.R
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# load 
library(dplyr)
library(terra)
library(sf)
library(stringr)

setwd("C:/Users/blbouf/Sync/Paper2")

data_path <- file.path(".", "data")
fcover_path <- file.path(data_path, "src", "FCOVER")
catchment_path <- file.path(data_path, "src", "TC_catchment", "catchments.shp")


# list HLS data
land_hls_files <- list.files(fcover_path, recursive = TRUE, full.names = TRUE)
fcover_files <- land_hls_files[grepl("fCOVER_30m.tif$", land_hls_files)]
qc_files <- land_hls_files[grepl("QC_30m.tif$", land_hls_files)]

#catchments
catchments_shp <- st_read(catchment_path)

# set and create outfolder if DNE
outfolder <- file.path(data_path, "fCOVER_clean")
if (!dir.exists(outfolder)){
  dir.create(outfolder)
}

# do this once - fix files names and scale data 
lapply(fcover_files, 
       tidy_fcover_files,
       outfolder)

# stack raster and clip to catchment boundary 
new_files <- list.files(outfolder,
                        full.names = TRUE)

# STACK EM 
fcover_all <- rast(new_files)
fcover_all[fcover_all == 0] <- NA

# make stacked dir if DNE
stack_outpath <- file.path(outfolder, "fCOVER_stacked")
if (!(file.exists(stack_outpath))){
  dir.create(stack_outpath)
}

# write raster
writeRaster(fcover_all,
            file.path(stack_outpath, "fcover_stacked.tif"),
            overwrite = TRUE)

# convert to same crs as TC catchment 
fcover_all_bcalbers <- project(fcover_all, 
                            crs("EPSG:3005"))

catchments_shp <- st_transform(catchments_shp, crs=crs("EPSG:3005"))

# mask by polygons 
fcover_masked <- mask(fcover_all_bcalbers,
                   catchments_shp)
# write raster 
writeRaster(fcover_masked,
            file.path(stack_outpath, "fcover_stacked_TC_area.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

tidy_fcover_files <- function(file_name, outfolder){
  # function to get year_moont_int from filepath, filter out bad values based on 
  # QC layer, and scale by factor of 200 (based on LEAF documentation)
  
  qc_filename <- gsub("fCOVER", "QC", file_name)
  
  year <- str_split(file_name, "/")[[1]][5]
  year <- str_split(year, "_")[[1]][2]
  
  month <- str_split(file_name, "/")[[1]][6]
  month <- str_split(month, "_")[[1]][3]
  
  month_df <- as.data.frame(matrix(data=NA, nrow=12, ncol=2))
  names(month_df) <- c("month", "month_int")
  month_df$month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  month_df$month_int <- c("01", "02","03","04","05","06","07","08","09","10","11","12")
  
  month_index <- which(month_df$month == month)
  month_int <- month_df$month_int[month_index]
  
  fcover_rast <- rast(file_name)
  qc_rast <- rast(qc_filename)
  
  # convert to bc_albers 
  #lai_rast_albers <- terra::project(lai_rast, 
  #                                  crs("EPSG:3005"))
  #qc_rast_albers <- terra::project(qc_rast,
  #                                 crs("EPSG:3005"))
  
  qc_rast[qc_rast != 0] <- NA
  
  qc_rast[qc_rast == 0] <- 1  
  
  fcover_fixed <- fcover_rast/200
  fcover_fixed <- fcover_fixed*qc_rast
  
  names(fcover_fixed) <- paste0("fCOVER_", year, "_", month_int)
  
  writeRaster(fcover_fixed, 
              file.path(outfolder, paste0("fCOVER_", year, "_", month_int, ".tif")),
              overwrite = TRUE)
}
