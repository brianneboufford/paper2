# -----------------------------------------------------------------------------_
# make RVP and RVH file info for ALL FOREST model run -- 
# RVC and RVI stay the same for all model runs with no landcover change 
# running from jan 1 1980 to jan 1 2023
# 
# copied from 01_rvn_files_allforested.R
# 
# Dec 2nd, 2025
# modified Feb 19, 2026 
# Brianne Boufford
#
# current wd = C:\Users\blbouf\Sync\TrappingCreek\raven-runs\_Trapping_LAI\Trapping_model_runs_reprod
# github repo: Trapping_model_runs
# updated March 22 , 2026
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)
library(RavenR)

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

# set up output directory ****************
# output path
outpath <- file.path(".", "data", "rvp_rvh_data", "all_forested_Mar22")

if (!file.exists(outpath)){
  dir.create(outpath, recursive = TRUE)
}

# source functions for making rvp file data 
source(file.path(".", "github", "paper2", "raven-baseline-calibration", "Raven_rvp_functions_paper2.R"))

# ------------------------------------------------------------------------------
# RVH FILE
#   - done directly with RavenR functions
#   - update all veg and land use classes to just the forest cover class
# ------------------------------------------------------------------------------

# 2023 hru for rvh file 
hru_data <- st_read(file.path(".", "data", "LCC_HRU_files", "Feb19", "HRUs_1923_2023", "HRU_2023.shp")) # was dec 2nd

# template rvh
rvh_path <- file.path(".", "raven-runs", "Baseline2", "Trapping_HRU_baseline.rvh")

# read in template rvh and separate out HRU table 
rvh <- rvn_rvh_read(rvh_path)
hru_table <- rvh$HRUtable

# change landse and Vegetation to forest class 
# hru_table$LandUse[grepl("_idf", hru_table$LandUse)] <- "M_idf"
# hru_table$LandUse[grepl("_essf", hru_table$LandUse)] <- "M_essf"
# hru_table$LandUse[grepl("_ms", hru_table$LandUse)] <- "M_ms"
# 
# hru_table$Vegetation[grepl("_idf", hru_table$Vegetation)] <- "M_idf"
# hru_table$Vegetation[grepl("_essf", hru_table$Vegetation)] <- "M_essf"
# hru_table$Vegetation[grepl("_ms", hru_table$Vegetation)] <- "M_ms"

# instead change the forested types to eitehr M_g0 or M_g1
hru_table <- hru_table %>%
  mutate(
    LandUse = if_else(
      Elevation >= 1560 & !LandUse %in% c("WET_LAND", "SHRUB", "ALPINE"),
      "M_g1",
      LandUse
    ),
    Vegetation = if_else(
      Elevation >= 1560 & !LandUse %in% c("WET_LAND", "SHRUB", "ALPINE"),
      "M_g1",
      Vegetation
    ),
    LandUse = if_else(
      Elevation < 1560 & !LandUse %in% c("WET_LAND", "SHRUB", "ALPINE"),
      "M_g0",
      LandUse
    ),
    Vegetation = if_else(
      Elevation < 1560 & !LandUse %in% c("WET_LAND", "SHRUB", "ALPINE"),
      "M_g0",
      Vegetation
    )
  ) 

# write as new rvh file
rvn_rvh_write(file.path(outpath, "Trapping_all_forest.rvh"), 
              SBtable = rvh$SBtable,
              HRUtable = hru_table)

# ------------------------------------------------------------------------------
# RVP FILE
#   - done manually .. no RavenR function
#   - 5 sections: VegetationClasses, SeasonalCanopyLAI, VegetationParameterList, 
#                 LandUse Classes, LandUseParameterList ; others are OK are default
# ------------------------------------------------------------------------------

# read all forested LAI data 
old_all_forested_lai <- read.csv(file.path("..", "TrappingCreek", "raven-runs", "_Trapping_LAI", 
                                        "data", "monthly_LAI", "median_seasonal_curve_lai_recovery_dec2.csv")) %>%
  dplyr::select(-c("med_apr_lai"))
non_forested_lai <- old_all_forested_lai[old_all_forested_lai$age_class %in% c("WETLAND", "ALPINE", "SHRUB"), ]

all_forested_lai <- read.csv(file.path(".", "data", "median_forest_params_LAI_by_age", 
                                       "median_LAI_ALS_ELEV_recovery_mar22.csv")) %>% #was mar19
  pivot_wider(values_from = lai, names_from = month_num) %>%
  dplyr::select(-c("lai_grp"))
names(all_forested_lai) <- c("age_class", "med_may_lai", "med_june_lai", "med_july_lai",
                             "med_aug_lai", "med_sept_lai", "med_oct_lai", "med_winter_lai")
all_forested_lai <- rbind(all_forested_lai, non_forested_lai)

# read in forest cover and height data 
fc_h <- read.csv(file.path(".", "data", "med_forest_params_curve_fitted", 
                           "height_fcover_med_Scurve_mar22.csv")) #%>% # was feb17
  #dplyr::select(-c("ages"))

# make small df for non-forested classes 
fc_h_non_forested <- data.frame(
  age_class = c("ALPINE", "WET_LAND", "SHRUB"),
  med_fc = c(0.0, 0.5, 0.6),
  med_h = c(0, 0, 1) 
)

# join together 
fc_h <- rbind(fc_h, fc_h_non_forested)

# fix spelling of Wetland veg class
all_forested_lai$age_class <- stringr::str_replace(all_forested_lai$age_class, "WETLAND", "WET_LAND")

# # : VegetationClasses
# veg_classes <- c("WET_LAND", "ALPINE", "SHRUB",
#                  "D_idf", "R1_idf", "R2_idf", "R3_idf", "R4_idf", "R5_idf", "M_idf",
#                  "D_ms", "R1_ms", "R2_ms", "R3_ms", "R4_ms", "R5_ms", "M_ms",
#                  "D_essf", "R1_essf", "R2_essf", "R3_essf", "R4_essf", "R5_essf", "M_essf")

# : VegetationClasses
veg_classes <- unique(all_forested_lai$age_class)

vegclasses_df <- data.frame(
  ID = veg_classes,
  MAX_HT = rep("_DEFAULT", length(veg_classes)),
  MAX_LAI = rep("_DEFAULT", length(veg_classes)),
  MAX_LEAF_COND = rep("_DEFAULT", length(veg_classes))
)

# update all to be 0 
vegclasses_df$MAX_LEAF_COND <- 0

# grab july LAI for peak LAI 
vegclasses_df <- vegclasses_df %>%
  left_join(all_forested_lai %>% dplyr::select(age_class, med_july_lai), by = c("ID" = "age_class")) %>%
  mutate(MAX_LAI = ifelse(!is.na(med_july_lai), med_july_lai, MAX_LAI)) %>%
  dplyr::select(-med_july_lai) 

vegclasses_df <- vegclasses_df %>% 
  left_join(fc_h %>% dplyr::select(c("age_class", "med_h")), by = c("ID" = "age_class"))

vegclasses_df$MAX_LAI <- round(vegclasses_df$MAX_LAI, digits = 2)
vegclasses_df$MAX_HT <- round(as.numeric(vegclasses_df$med_h), digits = 1)
vegclasses_df <- vegclasses_df %>% dplyr::select(-c("med_h"))

# :SeasonalCanopyLAI
seasonalcanopylai <- get_seasonal_canopy_lai_all_forest(all_forested_lai)

# : VegetationParameterList
vegparam_df <- data.frame(
  ID = veg_classes,
  MAX_CAPACITY = rep("_DEFAULT", length(veg_classes)),
  MAX_SNOW_CAPACITY = rep("_DEFAULT", length(veg_classes))
)

# max capacity stays as _DEFAULT for all 
vegparam_df$MAX_SNOW_CAPACITY <- get_max_snow_c(veg_classes)

# : Landuse classes
landuse_df <- data.frame(
  ID = veg_classes,
  IMPERM = rep("_DEFAULT", length(veg_classes)),
  FOREST_COVER = rep("_DEFAULT", length(veg_classes))
)

#landuse_df[, c(2,3)] <- get_LU_classes(veg_classes)
landuse_df$IMPERM <- 0.0 

landuse_df <- landuse_df %>% 
  left_join(fc_h %>% dplyr::select(c("age_class", "med_fc")), by = c("ID" = "age_class"))

landuse_df <- landuse_df %>% dplyr::select(-c("FOREST_COVER"))
names(landuse_df) <- c("ID", "IMPERM", "FOREST_COVER")
landuse_df$FOREST_COVER <- round(landuse_df$FOREST_COVER, 0)/100

# : LanduseParameterList 
landuseparam_df <- data.frame(
  ID = veg_classes,
  MELT_FACTOR = rep("_DEFAULT", length(veg_classes)),
  MIN_MELT_FACTOR = rep("_DEFAULT", length(veg_classes)),
  HBV_MELT_FOR_CORR = rep("_DEFAULT", length(veg_classes)),
  REFREEZE_FACTOR = rep("_DEFAULT", length(veg_classes)),
  HBV_MELT_ASP_CORR = rep("_DEFAULT", length(veg_classes)),
  PRIESTLYTAYLOR_COEFF = rep("_DEFAULT", length(veg_classes))
)

landuseparam_df$HBV_MELT_FOR_CORR <- get_LU_params(veg_classes)[,1]
landuseparam_df$PRIESTLYTAYLOR_COEFF <- get_LU_params(veg_classes)[, 2]

# ------------------------------------------------------------------------------
# write all tables 
# ------------------------------------------------------------------------------
write.csv(vegclasses_df,
          file.path(outpath, "vegclasses.csv"),
          row.names = FALSE)

write.csv(seasonalcanopylai, 
          file.path(outpath, "seasonalcanopyLAI.csv"),
          row.names = FALSE)

write.csv(vegparam_df, 
          file.path(outpath, "vegparam.csv"),
          row.names = FALSE)

write.csv(landuse_df, 
          file.path(outpath, "landuse.csv"),
          row.names = FALSE)

write.csv(landuseparam_df, 
          file.path(outpath, "landuseParam.csv"),
          row.names = FALSE)
