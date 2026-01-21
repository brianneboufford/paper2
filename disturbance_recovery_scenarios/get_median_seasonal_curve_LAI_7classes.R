# -----------------------------------------------------------------------------_
# LAI seasonal curves (monthly values) based on 7-class system for paper 2 
# 
# 
# Nov 27th, 2025
# Brianne Boufford
#
# current wd = C:\Users\blbouf\Sync\TrappingCreek\LAI_analysis\scripts\LAI_recovery
# github repo: LAI_recovery 
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# set wd 
setwd(file.path("C:","Users", "blbouf", "Sync", "TrappingCreek", "LAI_analysis", 
                "scripts", "LAI_recovery"))

output_path <- file.path("..","..","data", "SOM_outputs", "recovery_results")
figs_path <- file.path("..","..","data", "SOM_outputs", "recovery_results", "figs")

# zon ekey 
zone_key <- read.csv(file.path("..","..","data", "SOM_outputs","zone_key.csv"))
zone_key <- zone_key[c(1,2,5), ]

# path to lai data 
sampled_lai_data_path <- file.path("..","..","data", "SOM_outputs", "recovery_results", 
                                   "sampled_lai_all_data_2014_2021_june5.csv")
#sampled_lai_data_path <- file.path("..","..","data", "SOM_outputs", "recovery_results", 
#                                   "lai_all_data_2014_2020_nov27.csv")

# read LAI data 
sampled_lai_data <- read.csv(sampled_lai_data_path)
sampled_lai_data <- merge(sampled_lai_data, zone_key, by="id")

# look up table
zoned_classes <- list(
  MS = data.frame(
    min_age = c(0, 5, 10, 20, 30, 40, 50),
    max_age = c(5, 10, 20, 30, 40, 50, Inf),
    class   = c("D_ms", "R1_ms", "R2_ms", "R3_ms", "R4_ms", "R5_ms", "M_ms")
  ),
  ESSF = data.frame(
    min_age = c(0, 5, 10, 20, 30, 40, 50),
    max_age = c(5, 10, 20, 30, 40, 50, Inf),
    class   = c("D_essf", "R1_essf", "R2_essf", "R3_essf", "R4_essf", "R5_essf", "M_essf")
  ),
  IDF = data.frame(
    min_age = c(0, 5, 10, 20, 30, 40, 50),
    max_age = c(5, 10, 20, 30, 40, 50, Inf),
    class   = c("D_idf", "R1_idf", "R2_idf", "R3_idf", "R4_idf", "R5_idf", "M_idf")
  )
)

# ------------------------------------------------------------------------------

# CHAT GPT --- TEST 
df <- sampled_lai_data

group_by_zone_ageclass_lai <- function(df, zoned_classes) {
  
  # df: data frame with columns: zone, age, winter, apr, ..., oct
  # zoned_classes: named list of data frames as described
  
  # Function to assign age class based on zone-specific breakpoints
  assign_class <- function(age, zone) {
    class_df <- zoned_classes[[as.character(zone)]]
    if (is.null(class_df)) return(zone)
    idx <- which(age >= class_df$min_age & age < class_df$max_age)
    if (length(idx) == 0) return(zone)
    class_df$class[idx[1]]
  }
  
  # Apply the function rowwise to assign age class
  df <- df %>%
    rowwise() %>%
    mutate(age_class = assign_class(age, frst_cl)) %>%
    ungroup()
  
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

result <- group_by_zone_ageclass_lai(df, zoned_classes)

# result_nowinter <- result %>% select(-c("winter_lai", "april_lai", "oct_lai"))

# write out lai data
result_path <- file.path(file.path("..", "..", "..", "raven-runs", "_Trapping_LAI", "data", "monthly_LAI"))
write.csv(result, 
          file.path(result_path, "median_seasonal_curve_lai_recovery_dec2.csv"), # was jun5
          row.names = FALSE)
