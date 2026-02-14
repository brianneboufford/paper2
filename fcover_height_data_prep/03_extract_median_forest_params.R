# -----------------------------------------------------------------------------_
# script to get median veg height, LAI, and fcover by age 
#
#
# adapated from 004_mean_LAI_by_age_advanced_LAI_recovery.R
# 2 -- after 00 join fcover height age and smaple 
# date created: february 4, 2026
# date last modified: february 9, 2026
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(smoothr)
library(stringr)
library(ggplot2)
library(future)
library(future.apply)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# fcover, veg height, and LAI
veg_params_path <- file.path(".", "data", "forest_params_by_age", "sampled_veg_params_byLAI_GRP_2015_2021_feb10.csv") # fcover and height
als_veg_params_path <- file.path(".", "data", "forest_params_by_age", "sampled_fc_height_2015_feb10.csv") # fcover and height

zone_path <- file.path(".", "data", "src", "ntems", "zone_key.csv")

# output path
result_path <- file.path(file.path(".", "data", "median_forest_params_LAI_by_age")) 

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# BEC zone key 
zone_key <- read.csv(zone_path)
zone_key <- zone_key[c(1,2,5), ]

# read sampled forst param data (lai, height, fc)
veg <- read.csv(veg_params_path)
veg_data <- merge(veg, zone_key, by="id")

# read sampled als forest param data (height, fc)
als_veg <- read.csv(als_veg_params_path)
als_veg_data <- merge(als_veg, zone_key, by="id")

# ------------------------------------------------------------------------------
# change group label to string 
# ------------------------------------------------------------------------------
# forest params 
veg_data$lai_grp[veg_data$lai_grp == 0] <- "g0"
veg_data$lai_grp[veg_data$lai_grp == 1] <- "g1"

# als forst params 
als_veg_data$lai_grp[als_veg_data$lai_grp == 0] <- "g0"
als_veg_data$lai_grp[als_veg_data$lai_grp == 1] <- "g1"
als_veg_data <- als_veg_data[als_veg_data$frst_cl %in% c("ESSF", "IDF", "MS"), ]

# ------------------------------------------------------------------------------
# define age classes
# ------------------------------------------------------------------------------

# look up table for BEC zone recovery class distinction
bec_zoned_classes <- list(
  MS = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_ms", "R1_ms", "R2_ms", "R3_ms", "R4_ms", "R5_ms", 
                "R6_ms", "R7_ms", "R8_ms", "R9_ms", "M_ms")
  ),
  ESSF = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_essf", "R1_essf", "R2_essf", "R3_essf", "R4_essf", "R5_essf", 
                "R6_essf", "R7_essf", "R8_essf", "R9_essf","M_essf")
  ),
  IDF = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_idf", "R1_idf", "R2_idf", "R3_idf", "R4_idf", "R5_idf", 
                "R6_idf", "R7_idf", "R8_idf", "R9_idf","M_idf")
  )
)

# look up table basedc on LAI groups (from elevation LAI breakpoint)
elev_zoned_classes <- list(
  g0 = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_g0", "R1_g0", "R2_g0", "R3_g0", "R4_g0", "R5_g0", 
                "R6_g0", "R7_g0", "R8_g0", "R9_g0", "M_g0")
  ),
  g1 = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_g1", "R1_g1", "R2_g1", "R3_g1", "R4_g1", "R5_g1", 
                "R6_g1", "R7_g1", "R8_g1", "R9_g1","M_g1")
  )
)

# ------------------------------------------------------------------------------
# median LAI - BEC recovery classes
# ------------------------------------------------------------------------------

result_lai_bec <- med_LAI_by_ageclass(veg_data, bec_zoned_classes, "eco")

write.csv(result_lai_bec,
          file.path(result_path, "median_lai_BEC_recovery_feb10.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------
# median veg height and fcover - BEC recovery classes
# ------------------------------------------------------------------------------

result_fparam_bec <- med_fparams_by_ageclass(veg_data, bec_zoned_classes, "eco")

write.csv(result_fparam_bec, 
          file.path(result_path, "median_fparams_BEC_recovery_feb10.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------
# median LAI - Elevation recovery classes
# ------------------------------------------------------------------------------

result_lai_elev <- med_LAI_by_ageclass(veg_data, elev_zoned_classes, "elev")

write.csv(result_lai_elev,
          file.path(result_path, "median_lai_ELEV_recovery_feb10.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------
# median veg height and fcover - Elevation recovery classes
# ------------------------------------------------------------------------------

result_fparam_elev <- med_fparams_by_ageclass(veg_data,  elev_zoned_classes, "elev")

write.csv(result_fparam_elev, 
          file.path(result_path, "median_fparams_ELEV_recovery_feb10.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------
# median veg height and fcover from ALS - Elevation recovery classes
# ------------------------------------------------------------------------------

result_fparam_ALS_elev <- med_fparamsALS_by_ageclass(als_veg_data,  elev_zoned_classes, "elev")

write.csv(result_fparam_ALS_elev, 
          file.path(result_path, "median_fparams_ALS_ELEV_recovery_feb10.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

# LAI 
med_LAI_by_ageclass <- function(df, zoned_classes, class_type) {
  
  # df: data frame with columns: zone, age, winter, apr, ..., oct
  # zoned_classes: named list of data frames as described
  # class_type: either "eco" or "elev" to determine the df column to dictate the subzones for the age classification
  
  # Function to assign age class based on zone-specific breakpoints
  assign_class <- function(age, zone) {
    class_df <- zoned_classes[[as.character(zone)]]
    if (is.null(class_df)) return(zone)
    idx <- which(age >= class_df$min_age & age < class_df$max_age)
    if (length(idx) == 0) return(zone)
    class_df$class[idx[1]]
  }
  
  # Apply the function rowwise to assign age class
  if (class_type == "eco"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, frst_cl)) %>%
      ungroup()
    
  } else if(class_type == "elev"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, lai_grp)) %>%
      ungroup()
  }

  # take median LAI for each age_class 
  reclass_lai <- df %>% 
    group_by(age_class) %>%
    summarize(med_apr_lai = round(median(april_lai, na.rm = TRUE),2),
              med_may_lai = round(median(may_lai, na.rm = TRUE),2),
              med_june_lai = round(median(june_lai, na.rm = TRUE),2),
              med_july_lai = round(median(july_lai, na.rm = TRUE),2),
              med_aug_lai = round(median(aug_lai, na.rm = TRUE),2),
              med_sept_lai = round(median(sept_lai, na.rm = TRUE),2),
              med_oct_lai = round(median(oct_lai, na.rm = TRUE),2),
              med_winter_lai = round(median(med_winter_lai, na.rm = TRUE),2))
  
  # clean up some erouneously high april, and winter values or Missing values 
  # if may is missing, use June 
  # if may is > jul make may = october 
  # if april is > may then make april = may 
  # if winter is > oct then make make winter october 
  reclass_lai <- reclass_lai %>%
    mutate(med_may_lai = if_else(is.na(med_may_lai), med_june_lai, med_may_lai)) %>%
    
    mutate(med_may_lai = if_else(is.na(med_may_lai) | med_may_lai > med_july_lai, med_oct_lai, med_may_lai)) %>%
    
    mutate(med_apr_lai = if_else(is.na(med_apr_lai) | med_apr_lai > med_may_lai, med_may_lai, med_apr_lai)) %>%
    
    mutate(med_winter_lai = if_else(is.na(med_winter_lai) | med_winter_lai > med_oct_lai, med_may_lai, med_winter_lai))
  
  return(reclass_lai)
  
}

# MEDIAN FCOVER and HEIGHT (NTEMS) 
med_fparams_by_ageclass <- function(df, zoned_classes, class_type) {
  
  # df: data frame with columns: zone, age, winter, apr, ..., oct
  # zoned_classes: named list of data frames as described
  # class_type: either "eco" or "elev" to determine the df column to dictate the subzones for the age classification
  
  # Function to assign age class based on zone-specific breakpoints
  assign_class <- function(age, zone) {
    class_df <- zoned_classes[[as.character(zone)]]
    if (is.null(class_df)) return(zone)
    idx <- which(age >= class_df$min_age & age < class_df$max_age)
    if (length(idx) == 0) return(zone)
    class_df$class[idx[1]]
  }
  
  # Apply the function rowwise to assign age class
  if (class_type == "eco"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, frst_cl)) %>%
      ungroup()
    
  } else if(class_type == "elev"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, lai_grp)) %>%
      ungroup()
  }
  
  # take median LAI for each age_class 
  reclass_fparams <- df %>% 
    group_by(age_class) %>%
    summarize(med_h_ntems = round(median(ntems_height, na.rm = TRUE),2),
              med_fc_ntems = round(median(ntems_fcover, na.rm = TRUE),2))
  
  return(reclass_fparams)
  
}

# MEDIAN FCOVER and HEIGHT (ALS) 
med_fparamsALS_by_ageclass <- function(df, zoned_classes, class_type) {
  
  # df: data frame with columns: zone, age, winter, apr, ..., oct
  # zoned_classes: named list of data frames as described
  # class_type: either "eco" or "elev" to determine the df column to dictate the subzones for the age classification
  
  # Function to assign age class based on zone-specific breakpoints
  assign_class <- function(age, zone) {
    class_df <- zoned_classes[[as.character(zone)]]
    if (is.null(class_df)) return(zone)
    idx <- which(age >= class_df$min_age & age < class_df$max_age)
    if (length(idx) == 0) return(zone)
    class_df$class[idx[1]]
  }
  
  # Apply the function rowwise to assign age class
  if (class_type == "eco"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, frst_cl)) %>%
      ungroup()
    
  } else if(class_type == "elev"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, lai_grp)) %>%
      ungroup()
  }
  
  # take median LAI for each age_class 
  reclass_fparams <- df %>% 
    group_by(age_class) %>%
    summarize(med_h_als = round(median(height, na.rm = TRUE),2),
              med_fc_als = round(median(fc, na.rm = TRUE),2))
  
  return(reclass_fparams)
  
}


