# -----------------------------------------------------------------------------_
# fCOVER seasonal curves (monthly values) based on zoned classes derived from 
# growth curve by BEC zone 
# 
#  Novmeber 19, 2025
# adapted from June 5, 2025 script 004_mean_LAI_by_age_advanced_LAI_recovery
# Brianne Boufford
#
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

output_prod_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "average_monthly_fCOVER")
figs_path <- file.path(".", "data", "figs", "fCOVER_analysis")
sampled_data_path <- file.path(".", "data", "fCOVER_hls_als_analysis", "fCOVER_age",
                               "sampled_fCOVER_all_data_2014_2021_nov19.csv")

# zon ekey 
zone_key <- read.csv(file.path(".", "data", "src", "BEC","zone_key.csv"))
zone_key <- zone_key[c(1,2,5), ]

# read fCOVER data 
sampled_data <- read.csv(sampled_data_path)
sampled_data <- merge(sampled_data, zone_key, by="id")

# look up table
zoned_classes <- list(
  MS = data.frame(
    min_age = c(0, 5, 10, 37, 44),
    max_age = c(5, 10, 37, 44, Inf),
    class   = c("disturbed_ms", "early_recovery_ms", "late_recovery_ms", "young_ms", "old_ms")
  ),
  ESSF = data.frame(
    min_age = c(0, 5, 12, 40, 50),
    max_age = c(5, 12, 40, 50, Inf),
    class   = c("disturbed_essf", "early_recovery_essf", "late_recovery_essf", "young_essf", "old_essf")
  ),
  IDF = data.frame(
    min_age = c(0, 5, 8, 26, 33),
    max_age = c(5, 8, 26, 33, Inf),
    class   = c("disturbed_idf", "early_recovery__idf","late_recovery_idf", "young_idf", "old_idf")
  )
)

# CHAT GPT --- TEST 
df <- sampled_data

result <- group_by_zone_ageclass_fcover(df, zoned_classes)


# result_nowinter <- result %>% select(-c("winter_lai", "april_lai", "oct_lai"))

# write out lai data
write.csv(result, 
          file.path(output_prod_path, "median_seasonal_curve_fcover_nov19.csv"), # was jun5
          row.names = FALSE)

# ------------------------------------------------------------------------------

group_by_zone_ageclass_fcover <- function(df, zoned_classes) {
  
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
  reclass_fcover <- df %>% 
    group_by(age_class) %>%
    summarize(med_apr = round(median(april_fCOVER, na.rm = TRUE),2),
              med_may = round(median(may_fCOVER, na.rm = TRUE),2),
              med_june = round(median(june_fCOVER, na.rm = TRUE),2),
              med_july = round(median(july_fCOVER, na.rm = TRUE),2),
              med_aug = round(median(aug_fCOVER, na.rm = TRUE),2),
              med_sept = round(median(sept_fCOVER, na.rm = TRUE),2),
              med_oct = round(median(oct_fCOVER, na.rm = TRUE),2),
              med_winter = round(median(med_winter_fCOVER, na.rm = TRUE),2))
  
  # clean up some erouneously high april, and winter values or Missing values 
  # if may is missing, use June 
  # if may is > jul make may = october 
  # if april is > may then make april = may 
  # if winter is > oct then make make winter october 
  reclass_fcover <- reclass_fcover %>%
    mutate(med_may = if_else(is.na(med_may), med_june, med_may)) %>%
    
    mutate(med_may = if_else(is.na(med_may) | med_may > med_july, med_oct, med_may)) %>%
    
    mutate(med_apr = if_else(is.na(med_apr) | med_apr > med_may, med_may, med_apr)) %>%
    
    mutate(med_winter = if_else(is.na(med_winter) | med_winter > med_oct, med_may, med_winter))
  
  return(reclass_fcover)
  
}
