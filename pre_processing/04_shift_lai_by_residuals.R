# ------------------------------------------------------------------------------
# shift LEAF LAI values by residuals for whole time series
###############################################
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

# set working directory to the location of the github folder
setwd(file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", "scripts", "LiDAR_HLS_LAI_compare"))

# output path 
output_path <- file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", "data", 
                         "compare_LiDAR_HLS_LAI", "LAI_corrected_march22")

# main project path 
project_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek")

# ntems path 
ntems_path <- file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", "data", 
                        "compare_LiDAR_HLS_LAI", "from_bud", "catchments (1)", "11S")
cc_path <- file.path(ntems_path, "structure", "percentage_first_returns_above_2m")
osc_path <- file.path(ntems_path, "structure", "percentage_first_returns_above_mean")

hls_lai_path <- file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis",
                          "data", "L8_HLS_stacked", "lai_stacked_TC_area.tif")

# path to VLCE 2.0 data 
vlce_path <- file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek", "data", 
                       "raw", "11S", "VLCE2.0")

# path to lai residual values by class 
resid_path <- file.path(output_path, "lai_residual_summary_mar22.csv")

# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

shift_lai_by_residuals <- function(hls_layer, 
                                   hls_lai_path, 
                                   cc_path, 
                                   vlce_path,
                                   resid_path, 
                                   output_path){
  
  print(hls_layer)
  
  meta_i <- stringr::str_split(hls_layer, "_")[[1]]
  y <- meta_i[2]
  m <- meta_i[3]
  
  # read in hls LAI for year + month of interest
  hls_lai <- rast(hls_lai_path)
  hls_lai_i <- hls_lai[hls_layer]
  
  # fix 0 data 
  hls_lai_i[hls_lai_i == 0] <- NA

  lai_ext <- ext(hls_lai_i)

  if (y != 2023){

    # read in LCC from NTEMS
    vlce_files <- list.files(vlce_path, full.names = TRUE)
    vlce <- normalizePath(vlce_files[grepl(paste0(y, "_v20_v20.dat$"), vlce_files)]) %>%
      rast() %>% 
      terra::project(hls_lai)
    
    # canopy over from NTEMS
    cc_files <- list.files(cc_path, full.names = TRUE)
    cc <- cc_files[grepl(paste0(y, ".dat$"), cc_files)] %>%
      rast()
    cc <- cc/100
    cc <- terra::project(cc, hls_lai)

    # residual shift values
    lai_resid <- read.csv(resid_path)

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
    df <- c(hls_lai_i, vlce, cc_classified) %>%
      as.data.frame(xy = TRUE) #%>%
      #na.omit()

    # make coniferous class key column
    conif <- df[df$category == "Coniferous", ]
    conif$category <- paste0(conif$category, "_", conif$cc_class)

    # split into coniferous and other LC classes
    other_lc <- df[df$category != "Coniferous", ]

    # join rejoin data
    joined_data <- rbind(conif, other_lc)
    joined_data <- merge(joined_data, lai_resid, by="category", all = TRUE)
    joined_data[is.na(joined_data$mean_diff_lai), "median_diff_lai"] <- 0

    # adjust LAI values
    joined_data$LAI_adj <- joined_data[, hls_layer] + joined_data[, "median_diff_lai"]
    joined_data <- joined_data %>% dplyr::select(all_of(c("x", "y", hls_layer, "LAI_adj")))

    r_template <- hls_lai_i

    # Convert the dataframe to SpatVector
    v <- vect(joined_data, geom = c("x", "y"), crs = crs(hls_lai_i))

    # Rasterize the vector onto the template grid
    r_new <- rasterize(v, r_template, field = "LAI_adj", fun = "mean")

    r_new[r_new < 0] <- 0

    names(r_new) <- hls_layer

    writeRaster(r_new,
                file.path(output_path, paste0(hls_layer, ".tif")),
                overwrite = TRUE)

  } else {

    writeRaster(hls_lai_i,
                file.path(output_path, paste0(hls_layer, ".tif")),
                overwrite = TRUE)
  }
  
}

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
hls_lai <- rast(hls_lai_path)

hls_layer_list <- names(hls_lai)

#plan(multisession, workers=10L)
# can't read .dat file with future anymore for some reason? - assuming it cant find the header file
# works for a single file but not as batch?
t1 <- Sys.time()
adjusted_lai <- lapply(hls_layer_list,
                shift_lai_by_residuals, 
                hls_lai_path = hls_lai_path,
                cc_path = cc_path,
                vlce_path = vlce_path,
                resid_path = resid_path,
                output_path = file.path(output_path, "LAI_corrected_monthly_mar22")) 
t2 <- Sys.time()

# need to figure out why these arent the same extent even after projecting 
# -- the adjusted ones get all squished -- smth is not right 
shifted_files <- list.files(file.path(output_path, "LAI_corrected_monthly_mar22"), 
                            full.names = TRUE)
shifted_rasters <- lapply(shifted_files, 
                          rast) %>%
  rast()

writeRaster(shifted_rasters, 
            filename = file.path(output_path, "adjusted_hls_lai_mar22.tif"),
            overwrite = TRUE)

