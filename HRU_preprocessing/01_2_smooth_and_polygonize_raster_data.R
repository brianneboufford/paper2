# ------------------------------------------------------------------------------
#
# smooth and polygonize DEM, slope, aspect layers
# lumping aspect and slope classes together to be more coarse and more closely ressemble Matts 
# to speed up model 
#
# author: Brianne Boufford 
# updated: November 25th, 2025
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# libraries 
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(plyr)
library(smoothr)
library(rmapshaper)
library(units)

# probably do N facing 
# S facing 
# and then E/W facing 
# high and low slope?

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

# path to src data
ntems_path <- file.path(".", "data", "src", "ntems")
bec_path <- file.path(".", "data", "src", "BEC")
lidar_derived_path <- file.path(".", "data", "src", "lidar_derived")
cb_path <- file.path(".","data", "src","consolidated_cutblocks", 
                      "consolidated_cutblocks_tc.shp")

# output path
outpath <- file.path(".", "data", "HRU_delineation")

# ---------------------------------------------
# functions 
# ---------------------------------------------

# function to round up to nearest 100
round_up_100 <- function(r){
  return(ceiling(r/100)*100)
}

# function to round up to nearest 10 for slope 
round_up_10 <- function(r){
  return(ceiling(r/10)*10)
}

# convert to multipolygons
split_multipolygon <- function(mpoly, i) {
  
  # Split into individual polygons
  polys <- st_cast(mpoly, "POLYGON", warn = FALSE)
  
  # Add columns
  polys$id_m <- i # id_m is for the mutlipolygon 
  polys$id_n <- seq_len(nrow(polys)) # id_n is for the single polygon within each mpoly
  
  return(polys)
}

# ------------------------------------------------------------------------------
# DEM 
# ------------------------------------------------------------------------------
dem_file <- list.files(lidar_derived_path,
                       pattern = "dem.tif$", 
                       full.names = TRUE)
# load dem
dem_src <- rast(dem_file)

# road to nearest 100
dem_100_bands <- app(dem_src, 
                     fun=round_up_100)

if(!file.exists(file.path(outpath, "DEM_interm"))){
  dir.create(file.path(outpath, "DEM_interm"))
}

writeRaster(dem_100_bands, 
            file.path(outpath, "DEM_interm", "dem_100_bands.tif"),
            overwrite = TRUE)

# smooth with modal focal window --> 25m smooth
dem_100_smooth <- focal(dem_100_bands,
                        fun = "modal", 
                        w = 5,
                        na.rm = TRUE,
                        filename = file.path(outpath, "DEM_interm", "dem_100_smooth_25m.tif"),
                        overwrite = TRUE)

# convert to polygon 
dem_poly <- as.polygons(dem_100_smooth) %>%
  st_as_sf() 

st_write(dem_poly, 
         file.path(outpath, "DEM_interm", "dem_poly.shp"),
         append = FALSE)

# split all poylgons
for (i in 1:length(dem_poly$geometry)){
  
  # subset to single multipolygon 
  mpoly <- dem_poly[i, ]
  
  # apply split function 
  mpoly_split <- split_multipolygon(mpoly = mpoly, 
                                    i = i)
  
  if (i == 1){
    split_polys <- mpoly_split
  } else {
    split_polys <- rbind(split_polys, mpoly_split)
  }
}

# add area column to polygons 
split_polys$area_m2 <- st_area(split_polys) %>% drop_units()

st_write(split_polys, 
         file.path(outpath, "DEM_interm", "dem_poly_split.shp"),
         append = FALSE)

# drop polygons that are < 1ha
split_polys_dropped <- split_polys[split_polys$area_m2 > 10000, ]

st_write(split_polys_dropped, 
         file.path(outpath, "DEM_interm", "dem_poly_biggerthan1ha.shp"),
         append = FALSE)

dem_filled <- split_polys_dropped %>% 
  fill_holes(threshold = 10000)

st_write(dem_filled, 
         file.path(outpath, "DEM_interm", "dem_filled_1ha.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# slope 
# ------------------------------------------------------------------------------
slope_file <- list.files(lidar_derived_path,
                       pattern = "slope.tif$", 
                       full.names = TRUE)
# load dem
slope_src <- rast(slope_file)

# this was used for slope_interm_3_classes
slope_classes <- matrix(c(
  0, 5, 1,
  5, 15, 2,
  15, Inf, 3
), ncol = 3, byrow = TRUE)

slope_classes <- matrix(c(
  0, 15, 1,
  15, Inf, 2
), ncol = 3, byrow = TRUE)

slope_reclass <- classify(slope_src, slope_classes)

if(!file.exists(file.path(outpath, "slope_interm"))){
  dir.create(file.path(outpath, "slope_interm"))
}

writeRaster(slope_reclass, 
            file.path(outpath, "slope_interm", "slope_reclass.tif"),
            overwrite = TRUE)

# smooth with modal focal window --> 25m smooth
slope_smooth <- focal(slope_reclass,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE, 
                      filename = file.path(outpath, "slope_interm", "slope_smooth1.tif"),
                      overwrite = TRUE)

slope_smooth <- focal(slope_smooth,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE,
                      filename = file.path(outpath, "slope_interm", "slope_smooth2.tif"),
                      overwrite = TRUE)

slope_smooth <- focal(slope_smooth,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE,
                      filename = file.path(outpath, "slope_interm", "slope_smooth3.tif"),
                      overwrite = TRUE)

# convert to polygon 
slope_poly <- as.polygons(slope_smooth) %>%
  st_as_sf() 

st_write(slope_poly, 
         file.path(outpath, "slope_interm", "slope_poly.shp"),
         append = FALSE)

# split all poylgons
for (i in 1:length(slope_poly$geometry)){
  
  # subset to single multipolygon 
  mpoly <- slope_poly[i, ]
  
  # apply split function 
  mpoly_split <- split_multipolygon(mpoly = mpoly, 
                                    i = i)
  
  if (i == 1){
    split_polys <- mpoly_split
  } else {
    split_polys <- rbind(split_polys, mpoly_split)
  }
}

st_write(split_polys, 
         file.path(outpath, "slope_interm", "slope_poly_split.shp"),
         append = FALSE)

# add area column to polygons 
split_polys$area_m2 <- st_area(split_polys) %>% drop_units()

# drop polygons that are < 1ha
split_polys_dropped <- split_polys[split_polys$area_m2 > 100000, ]

st_write(split_polys_dropped, 
         file.path(outpath, "slope_interm", "slope_poly_biggerthan1ha.shp"),
         append = FALSE)

slope_filled <- split_polys_dropped %>% 
  fill_holes(threshold = 100000)

slope_filled$area <- st_area(slope_filled$geometry)

st_write(slope_filled, 
         file.path(outpath, "slope_interm", "slope_filled_1ha.shp"),
         append = FALSE)

# ------------
# test 
# ------------

covered <- st_covered_by(split_polys, slope_filled, sparse = FALSE)
cov_sum <- rowSums(covered)

# remaining polys
rm_polys <- split_polys[cov_sum == 0, ]
rm_polys$area <- st_area(rm_polys$geometry)

st_write(rm_polys, 
         file.path(outpath, "slope_interm", "polys_holes.shp"),
         append = FALSE)


# --- Function to merge small polygons ---
merge_small_polygons <- function(rm_polys, slope_filled) {
   
  for (i in 1:length(rm_polys$area)){
    
    small_poly_i <- rm_polys[i, ]
    
    neighbors <- st_filter(slope_filled, small_poly_i, .predicate = st_touches)
    
    if (nrow(neighbors) != 0){
    
      borders <- st_intersection(st_boundary(small_poly_i), st_boundary(neighbors))
      shared_length <- st_length(borders)
      
      # Pick neighbor with maximum shared border
      best_neighbor <- neighbors[which.max(shared_length), ]
      
      merged <- st_union(best_neighbor, small_poly_i) %>% 
        select(c("focal_modal", "id_m", "id_n", "area_m2", "area", "geometry"))
      
      # remove neighbour poly from dataset 
      slope_filled <- slope_filled[!(slope_filled$id_m == best_neighbor$id_m & slope_filled$id_n == best_neighbor$id_n), ]
      
      # Add merged polygon back
      slope_filled <- rbind(slope_filled, merged)
      
      # recalculate area
      slope_filled$area <- st_area(slope_filled)
    }
  }
  
  return(slope_filled)
}

# --- Run merge process ---
merged_polys <- merge_small_polygons(rm_polys = rm_polys,
                                     slope_filled = slope_filled)

st_write(merged_polys, 
         file.path(outpath, "slope_interm", "merged_polys.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# aspect 
# ------------------------------------------------------------------------------
asp_file <- list.files(lidar_derived_path,
                       pattern = "aspect.tif$", 
                       full.names = TRUE)
# load aspect
asp_src <- rast(asp_file)

# original classes - 4 classes, renamed to asp_interm_4_classes
asp_classes <- matrix(c(
  0, 45, 1,
  45, 135, 2,
  135, 225, 3, 
  225, 315, 4, 
  315, 361, 1
), ncol = 3, byrow = TRUE)

asp_classes <- matrix(c(
  0, 45, 1,
  45, 135, 2,
  135, 225, 3, 
  225, 315, 2, 
  315, 361, 1
), ncol = 3, byrow = TRUE)

asp_reclass <- classify(asp_src, asp_classes)

if(!file.exists(file.path(outpath, "asp_interm"))){
  dir.create(file.path(outpath, "asp_interm"))
}

writeRaster(asp_reclass, 
            file.path(outpath, "asp_interm", "aspect_reclass.tif"),
            overwrite = TRUE)

# smooth with modal focal window --> 25m smooth
asp_smooth <- focal(asp_reclass,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE, 
                      filename = file.path(outpath, "asp_interm", "aspect_smooth1.tif"),
                      overwrite = TRUE)

asp_smooth <- focal(asp_smooth,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE,
                      filename = file.path(outpath, "asp_interm", "aspect_smooth2.tif"),
                      overwrite = TRUE)

asp_smooth <- focal(asp_smooth,
                      fun = "modal", 
                      w = 21,
                      na.rm = TRUE,
                      filename = file.path(outpath, "asp_interm", "asp_smooth3.tif"),
                      overwrite = TRUE)

# convert to polygon 
asp_poly <- as.polygons(asp_smooth) %>%
  st_as_sf() 

st_write(asp_poly, 
         file.path(outpath, "asp_interm", "aspect_poly.shp"),
         append = FALSE)

# split all poylgons
for (i in 1:length(asp_poly$geometry)){
  
  # subset to single multipolygon 
  mpoly <- asp_poly[i, ]
  
  # apply split function 
  mpoly_split <- split_multipolygon(mpoly = mpoly, 
                                    i = i)
  
  if (i == 1){
    split_polys <- mpoly_split
  } else {
    split_polys <- rbind(split_polys, mpoly_split)
  }
}

st_write(split_polys, 
         file.path(outpath, "asp_interm", "asp_poly_split.shp"),
         append = FALSE)

# add area column to polygons 
split_polys$area_m2 <- st_area(split_polys) %>% drop_units()

# drop polygons that are < 1ha
split_polys_dropped <- split_polys[split_polys$area_m2 > 100000, ]

st_write(split_polys_dropped, 
         file.path(outpath, "asp_interm", "asp_poly_biggerthan1ha.shp"),
         append = FALSE)

asp_filled <- split_polys_dropped %>% 
  fill_holes(threshold = 100000)

asp_filled$area <- st_area(asp_filled$geometry)

st_write(asp_filled, 
         file.path(outpath, "asp_interm", "asp_filled_1ha.shp"),
         append = FALSE)

covered <- st_covered_by(split_polys, asp_filled, sparse = FALSE)
cov_sum <- rowSums(covered)

# remaining polys
rm_polys <- split_polys[cov_sum == 0, ]
rm_polys$area <- st_area(rm_polys$geometry)

st_write(rm_polys, 
         file.path(outpath, "asp_interm", "asp_polys_holes.shp"),
         append = FALSE)

# --- Run merge process ---
merged_polys <- merge_small_polygons(rm_polys = rm_polys,
                                     slope_filled = asp_filled)

st_write(merged_polys, 
         file.path(outpath, "asp_interm", "asp_merged_polys.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# Landcover
# ------------------------------------------------------------------------------
vlce_file <- list.files(ntems_path, 
                          pattern = "VLCE2.0.tif$", 
                          full.names = TRUE)

vlce_key <- list.files(ntems_path, 
                       pattern = "vlce.csv$", 
                       full.names = TRUE)

vlce <- rast(vlce_file)vlce <- rast(TRUEvlce_file)
vlce_key <- read.csv(vlce_key)


# ------------------------------------------------------------------------------
# BEC
# -----------------------------------------------------------------------------
bec_file <- list.files(bec_path, 
                        pattern = "bec_tc_clipped.shp$",
                        full.names = TRUE)

bec <- st_read(bec_file)
 
# ESSF is the only BEC zone with multiple polygons
essf <- bec[bec$ZONE == "ESSF", ]
essf_1 <- essf[1, ]
essf_union <- st_union(essf)
essf_1$geometry <- essf_union

idf <- bec[bec$ZONE == "IDF", ]
ms <- bec[bec$ZONE == "MS", ]

bec_new <- rbind(idf, ms) %>%
  rbind(essf_1)

st_write(bec_new, 
         file.path(bec_path, "bec_zones_tc.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# clipped consolidated cutblocks layer  
# ------------------------------------------------------------------------------
cb <- st_read(cb_path)

# split all poylgons
for (i in 1:length(cb$geometry)){
  
  # subset to single multipolygon 
  mpoly <- cb[i, ]
  
  # apply split function 
  mpoly_split <- split_multipolygon(mpoly = mpoly, 
                                    i = i)
  
  if (i == 1){
    split_polys <- mpoly_split
  } else {
    split_polys <- rbind(split_polys, mpoly_split)
  }
}

cb_split <- split_polys 

st_write(cb_split, 
         file.path(outpath, "cb_interm", "cb_split.shp"),
         append = FALSE)


# add area column to polygons 
cb_split$area_m2 <- st_area(cb_split) %>% drop_units()

# drop polygons that are < 1ha
cb_dropped <- cb_split[cb_split$area_m2 > 10000, ]

st_write(cb_dropped, 
         file.path(outpath, "cb_interm", "cb_poly_biggerthan1ha.shp"),
         append = FALSE)

cb_filled <- cb_dropped %>% 
  fill_holes(threshold = 10000)

cb_filled$area_m2 <- st_area(cb_filled$geometry)

st_write(cb_filled, 
         file.path(outpath, "cb_interm", "cb_filled_1ha.shp"), # use this for now
         append = FALSE)

# cb_covered <- st_covered_by(cb_split, cb_filled, sparse = FALSE)
# cb_cov_sum <- rowSums(cb_covered)
# 
# # remaining polys
# cb_rm_polys <- split_polys[cb_cov_sum == 0, ]
# cb_rm_polys$area <- st_area(cb_rm_polys$geometry)
# 
# st_write(cb_rm_polys, 
#          file.path(outpath, "cb_interm", "cb_polys_holes.shp"),
#          append = FALSE)
# 
# # --- Run merge process --- # this needs a column called focal modal 
# cb_merged_polys <- merge_small_polygons(rm_polys = cb_rm_polys,
#                                      slope_filled = cb_filled)
# 
# st_write(cb_merged_polys, 
#          file.path(outpath, "cb_interm", "cb_merged_polys.shp"),
#          append = FALSE)
