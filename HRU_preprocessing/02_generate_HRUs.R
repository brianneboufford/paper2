# ------------------------------------------------------------------------------
#
# generate new BEC polygons
# 
# author: Brianne Boufford 
# date: November 5th, 2025
# updated: November 26th, 2025 - to use 3 aspect classes and 2 slope classes 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# packages 
# ------------------------------------------------------------------------------

# libraries 
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(units)

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

# path to src data
ntems_path <- file.path(".", "data", "src", "ntems")
bec_path <- file.path(".", "data", "src", "BEC")
interm_path <- file.path(".", "data", "HRU_delineation")
fire_path <- file.path(".", "data", "src", "fire_perimeter")
cb_path <- file.path(".", "data", "HRU_delineation", "cb_interm")
catch_path <- file.path(".", "data", "src", "TC_catchment", "catchments.shp")
og_hru_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data", "src", "HRU_paper1", "HRU_Poly_attributes.shp")

# output path
outpath <- file.path(".", "data", "HRU_delineation")

# ---------------
# functions
# ---------------
# convert to multipolygons
split_multipolygon <- function(mpoly, i) {
  
  # Split into individual polygons
  polys <- st_cast(mpoly, "POLYGON", warn = FALSE)
  
  # Add columns
  polys$id_m <- i # id_m is for the mutlipolygon 
  polys$id_n <- seq_len(nrow(polys)) # id_n is for the single polygon within each mpoly
  
  return(polys)
}

merge_small_polygons <- function(rm_polys, filled) {
  
  for (i in 1:length(rm_polys$area_m2)){
    
    small_poly_i <- rm_polys[i, ]
    
    neighbors <- st_filter(filled, small_poly_i, .predicate = st_touches)
    
    if (nrow(neighbors) != 0){
      
      borders <- st_intersection(st_boundary(small_poly_i), st_boundary(neighbors))
      shared_length <- st_length(borders)
      
      # Pick neighbor with maximum shared border
      best_neighbor <- neighbors[which.max(shared_length), ]
      
      merged <- st_union(best_neighbor, small_poly_i$geometry)
      
      # remove neighbour poly from dataset 
      filled <- filled[!(filled$id_m == best_neighbor$id_m & filled$id_n == best_neighbor$id_n), ]
      
      # Add merged polygon back
      filled <- rbind(filled, merged)
      
      # recalculate area
      filled$area_m2 <- st_area(filled)
    }
  }
  
  return(filled)
}

my_union <- function(a,b) {
  #
  # function doing a real GIS union operation such as in QGIS or ArcGIS
  #
  # a - the first sf
  # b - the second sf
  #
  st_agr(a) = "constant"
  st_agr(b) = "constant"
  op1 <- st_difference(a, st_union(b))
  op2 <- st_difference(b, st_union(a))
  op3 <- st_intersection(b, a)
  union <- dplyr::bind_rows(op1, op2, op3) 
  
  return(union)
}

# split and merge polygons 
split_and_merge_polys <- function(path_to_union_shp, threshold = 5000, 
                                  outpath_to_union_shp){
  
  union <- st_read(path_to_union_shp)
  
  # split all poylgons
  for (i in 1:length(union$geometry)){
    
    # subset to single multipolygon 
    mpoly <- union[i, ]
    
    # apply split function 
    mpoly_split <- split_multipolygon(mpoly = mpoly, 
                                      i = i)
    
    if (i == 1){
      split_polys <- mpoly_split
    } else {
      split_polys <- rbind(split_polys, mpoly_split)
    }
  }
  
  split <- split_polys 
  
  # add area column to polygons 
  split$area_m2 <- st_area(split) %>% drop_units()
  
  # drop polygons that are < 1/2 ha
  dropped <- split[split$area_m2 > threshold, ]
  
  # fill holes
  filled <- dropped %>% 
    fill_holes(threshold = threshold)
  
  # redo area
  filled$area_m2 <- st_area(filled$geometry)
  
  covered <- st_covered_by(split, filled, sparse = FALSE)
  cov_sum <- rowSums(covered)
  
  # remaining polys
  rm_polys <- split[cov_sum == 0, ]
  rm_polys$area_m2 <- st_area(rm_polys$geometry)
  
  merged_polys <- merge_small_polygons(rm_polys = rm_polys,
                                       filled = filled)
  st_write(merged_polys,
           outpath_to_union_shp,
           append = FALSE)
}

# ---------------
# read in data, reudce columns to just attributes I want to keep, then rewrite to interm 
# folders
# ---------------

# BEC
bec <- st_read(file.path(bec_path, "bec_zones_tc.shp")) %>% 
  select(c("ZONE", "SUBZONE", "geometry"))
st_write(bec,
         file.path(interm_path,"BEC_interm", "bec_for_overlay.shp"), 
         append = FALSE)

# DEM
dem <- st_read(file.path(interm_path, "DEM_interm", "dem_filled_1ha.shp")) %>% 
  select(c("fcl_mdl", "geometry")) 
names(dem) <- c("elev", "geometry")
st_write(dem,
         file.path(interm_path,"DEM_interm", "dem_for_overlay.shp"), 
         append = FALSE)

# slope 
slp <- st_read(file.path(interm_path, "slope_interm", "merged_polys.shp")) %>% 
  select(c("fcl_mdl", "geometry"))
names(slp) <- c("slp", "geometry")
st_write(slp,
         file.path(interm_path, "slope_interm", "slope_for_overlay.shp"), 
         append = FALSE)
# aspect 
asp <- st_read(file.path(interm_path, "asp_interm", "asp_merged_polys.shp")) %>% 
  select(c("fcl_mdl", "geometry"))
names(asp) <- c("asp", "geometry")
st_write(asp,
         file.path(interm_path,"asp_interm", "asp_for_overlay.shp"), 
         append = FALSE)

# fire
fire <- st_read(file.path(fire_path, "fire_perim_tc.shp")) %>%
  select(c("FIRE_YE", "geometry"))
names(fire) <- c("fire_yr", "geometry")
st_write(fire,
         file.path(interm_path,"fire_interm", "fire_for_overlay.shp"), 
         append = FALSE)

cb <- st_read(file.path(cb_path, "cb_filled_1ha.shp")) %>% 
  select(c("harv_st", "geometry"))
names(cb) <- c("harv_yr", "geometry")
st_write(cb,
         file.path(interm_path, "cb_interm", "cb_for_overlay.shp"), 
         append = FALSE)

catch <- st_read(catch_path)

################################################################################
# 1. UNION FIRE AND CB
################################################################################

# union fire and cb 
hrus_new <- my_union(fire, cb) 
hrus_new$id <- 1:length(hrus_new$geometry)

st_write(hrus_new, 
         file.path(outpath, "interm_unions", "1_fire_cb.shp"),
         append = FALSE)

################################################################################
# 2. UNION ASP and SLP in Q
################################################################################
################################################################################
# 3. fix all small polygons 
################################################################################
path_to_union_shp <- file.path(interm_path, "interm_unions", "2_asp_slp.shp")
outpath_to_union_shp <- file.path(interm_path, "interm_unions", "2_asp_slp_filled.shp")
split_and_merge_polys(path_to_union_shp = path_to_union_shp, 
                      threshold = 10000, 
                      outpath_to_union_shp = outpath_to_union_shp)

################################################################################
# 4. UNION with DEM in Q
################################################################################
################################################################################
# 5. fix all small polygons 
################################################################################
path_to_union_shp <- file.path(interm_path, "interm_unions", "3_asp_slp_dem.shp")
outpath_to_union_shp <- file.path(interm_path, "interm_unions", "3_asp_slp_dem_filled.shp")
split_and_merge_polys(path_to_union_shp = path_to_union_shp, 
                      threshold = 5000, 
                      outpath_to_union_shp = outpath_to_union_shp)

################################################################################
# 6. clip 3_asp_slp_dem_filled catchments shapefile in QGIS
################################################################################


################################################################################
# 2. UNION with BEC in Q
################################################################################
################################################################################
# 3. fix all small polygons 
################################################################################
path_to_union_shp <- file.path(interm_path, "interm_unions", "4_asp_slp_dem_bec.shp")
outpath_to_union_shp <- file.path(interm_path, "interm_unions", "4_asp_slp_dem_bec_filled.shp")
split_and_merge_polys(path_to_union_shp = path_to_union_shp, 
                      threshold = 5000,
                      outpath_to_union_shp = outpath_to_union_shp)

################################################################################
# 9. union with 1_fire_cb_in Q
################################################################################

################################################################################
# 10. fix all small polygons in fire_cb_bec_dem_slp_asp_7
################################################################################
path_to_union_shp <- file.path(interm_path, "interm_unions", "5_asp_slp_dem_bec_fire_cb.shp")
outpath_to_union_shp <- file.path(interm_path, "interm_unions", "5_asp_slp_dem_bec_fire_cb_filled.shp")
split_and_merge_polys(path_to_union_shp = path_to_union_shp, 
                      threshold = 5000, 
                      outpath_to_union_shp = outpath_to_union_shp)

# redid this in QGIS -- keep fire boundaries better
# 1. convert multipart to singlepart 
# 2. select by attribute ($area < 10000)
# 3. eliminate selected polygons --> merge by largest common boundary

# ------------------------------------------------------------------------------
# read in HRUS 
# ------------------------------------------------------------------------------
hrus <- st_read(file.path(interm_path, "interm_unions", "5_asp_slp_dem_bec_fire_cb_filledQGIS_1ha.shp"))
hrus$area_m2 <- st_area(hrus$geometry)
hrus$id <- 1:length(hrus$geometry)

# extract average dem and aspect and slp in each HRU 
# keep class distinction 
# make col labels the same as old ones 
# add in alpine and wetland HRUs from last time

# ------------------------------------------------------------------------------
# read in og HRUs 
# ------------------------------------------------------------------------------

og_hrus <- st_read(og_hru_path)
non_forested_hrus <- og_hrus[og_hrus$VEG_CLA %in% c("WET_LAND", "SHRUB", "ALPINE"), ]
st_write(non_forested_hrus, 
         file.path(interm_path, "interm_unions", "non_forested_og_hrus.shp"),
         append = FALSE)

veg_class_dissolved_path <- file.path(interm_path, "interm_unions", "non_forested_dissolved_by_vegclass.shp")
veg_class_output_path <- file.path(interm_path, "interm_unions", "non_forested_merged.shp")

split_and_merge_polys(path_to_union_shp = veg_class_dissolved_path, 
                      threshold = 5000, 
                      outpath_to_union_shp = veg_class_output_path)

# read this into QGIS with 5_asp_slp....1ha.shp
# 0. intersected layers
# 1. convert multipart to singlepart 
# 2. select by attribut e($area < 10000)
# 3. eliminate selected polygons --> merge by largest common boundary


