# ------------------------------------------------------------------------------
# shift LEAF fCOVER values by residuals for whole time series
###############################################
# based on 002_shift_lai_by_residuals.R
# which was based on 
# based on LAI_residuals_shift_values.R and 02_shift_lai_by_residuals.R
# adapted to shift all months and to use new residual values
#
###############################################
#
# june 3, 2025
# ------------------------------------------------------------------------------

# library
library(terra)
library(dplyr)
library(sf)
library(stringr)
library(ggplot2)
library(future)
library(future.apply)

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
set.seed(1111)

setwd("C:/Users/blbouf/Sync/Paper2")

# output path
output_path <- file.path(".", "data", "fCOVER_hls_als_analysis")

# data paths 
lidar_fcover_path <- file.path(".", "data", "fCOVER_ALS", "percent_above_metrics.tiff")

cc_ntems_path <- file.path(".", "data", "src", "ntems", "structure.tif")

hls_cc_path <- file.path(".", "data", "fCOVER_clean", "fCOVER_stacked", "fcover_stacked_TC_area.tif")

vlce_path <- file.path(".", "data", "src", "ntems", "vlce_1984_2022")

bec_path <- file.path(".", "data", "src", "BEC", "bec_zones_tc.shp")

# path to lai residual values by class 
resid_path <- file.path(output_path, "cc_residual_summary.csv")

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
hls_cc <- rast(hls_cc_path)

hls_layer_list <- names(hls_cc)

#plan(multisession, workers=10L)
# can't read .dat file with future anymore for some reason? - assuming it cant find the header file
# works for a single file but not as batch?
t1 <- Sys.time()
adjusted_lai <- lapply(hls_layer_list,
                       shift_cc_by_residuals, 
                       hls_cc_path = hls_cc_path,
                       cc_ntems_path = cc_ntems_path,
                       vlce_path = vlce_path,
                       resid_path = resid_path,
                       output_path = file.path(output_path, "cc_corrected_monthly")) 
t2 <- Sys.time()

# need to figure out why these arent the same extent even after projecting 
# -- the adjusted ones get all squished -- smth is not right 
shifted_files <- list.files(file.path(output_path, "cc_corrected_monthly"), 
                            full.names = TRUE)

# fix that there is no april 2017 so I copied a null raster but need to change the layer name 
new_apr_2017 <- shifted_files[grepl("fCOVER_2015_11.tif", shifted_files)]
apr_2017_rast <- rast(new_apr_2017)
names(apr_2017_rast) <- "fCOVER_2017_04"
writeRaster(apr_2017_rast, file.path(output_path, "cc_corrected_monthly", "fCOVER_2017_04.tif"), overwrite = TRUE)

# read load list of files to now include april 2017 
shifted_files <- list.files(file.path(output_path, "cc_corrected_monthly"), 
                            full.names = TRUE)

# combine 
shifted_rasters <- lapply(shifted_files, 
                          rast) %>%
  rast()


writeRaster(shifted_rasters, 
            filename = file.path(output_path, "adjusted_hls_cc.tif"),
            overwrite = TRUE)
# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

shift_cc_by_residuals <- function(hls_layer, 
                                   hls_cc_path, 
                                   cc_ntems_path, 
                                   vlce_path,
                                   resid_path, 
                                   output_path){
  
  print(hls_layer)
  
  meta_i <- stringr::str_split(hls_layer, "_")[[1]]
  y <- meta_i[2]
  m <- meta_i[3]
  
  # read in hls fCOVER for year + month of interest
  hls_cc <- rast(hls_cc_path)
  hls_cc_i <- hls_lai[hls_layer]
  
  # fix 0 data 
  # hls_lai_i[hls_lai_i == 0] <- NA
  
  cc_ext <- ext(hls_cc_i)
  
  if (y != 2023){
    
    # read in LCC from NTEMS
    vlce_files <- list.files(vlce_path, full.names = TRUE)
    vlce <- normalizePath(vlce_files[grepl(paste0(y, "_v20_v20.dat$"), vlce_files)]) %>%
      rast() %>% 
      terra::project(hls_cc)
    
    # canopy over from NTEMS
    cc_ntems <- rast(cc_ntems_path)
    cc_ntems <- cc_ntems[[grepl(paste0(y), names(cc_ntems))]] 
    cc_ntems <- cc_ntems[[grepl("above_2m", names(cc_ntems))]]
    
    cc <- cc_ntems/100
    cc <- terra::project(cc, hls_cc)
    
    # residual shift values
    cc_resid <- read.csv(resid_path)
    
    # separate cc into classes of 1-20%, 20-50%, 50-80%, >80%
    cc_class_mat <- matrix(c(
      -Inf, 20, 0,
      20, 50, 1,
      50, 80, 2,
      80, Inf, 3
    ), ncol = 3, byrow = TRUE)
    
    # Reclassify canopy cover raster
    cc_classified <- terra::classify(cc, cc_class_mat)
    cc_classified[is.na(cc_classified)] <- 0
    names(cc_classified) <- "cc_class"
    
    # convert to df
    df <- c(hls_cc_i, vlce, cc_classified) %>%
      as.data.frame(xy = TRUE) #%>%
    #na.omit()
    
    # make coniferous class key column
    conif <- df[df$category == "Coniferous", ]
    conif$category <- paste0(conif$category, "_", conif$cc_class)
    
    # split into coniferous and other LC classes
    other_lc <- df[df$category != "Coniferous", ]
    
    # join rejoin data
    joined_data <- rbind(conif, other_lc)
    joined_data <- merge(joined_data, cc_resid, by="category", all = TRUE)
    joined_data[is.na(joined_data$mean_diff_cc), "median_diff_cc"] <- 0
    
    # adjust LAI values
    joined_data$fCOVER_adj <- joined_data[, hls_layer] + joined_data[, "median_diff_cc"]
    joined_data <- joined_data %>% select(all_of(c("x", "y", hls_layer, "fCOVER_adj")))
    
    r_template <- hls_cc_i
    
    # Convert the dataframe to SpatVector
    v <- vect(joined_data, geom = c("x", "y"), crs = crs(hls_cc_i))
    
    # Rasterize the vector onto the template grid
    r_new <- rasterize(v, r_template, field = "fCOVER_adj", fun = "mean")
    
    r_new[r_new < 0] <- 0
    r_new[r_new > 1] <- 1
    
    names(r_new) <- hls_layer
    
    writeRaster(r_new,
                file.path(output_path, paste0(hls_layer, ".tif")),
                overwrite = TRUE)
    
  } else {
    
    writeRaster(hls_cc_i,
                file.path(output_path, paste0(hls_layer, ".tif")),
                overwrite = TRUE)
  }
  
}



