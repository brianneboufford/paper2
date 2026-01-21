# -----------------------------------------------------------------------------_
# creating LandCover Change files
# updated to make land cover change files for updated vegetation classes with 
# data from May 9th from Matt
# 
# adapted from Dynamic_HRUs.R
# Date: November 10th, 2025
# updated December 2nd to use 7-class recovery classification
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(smoothr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# set wd 
setwd(file.path("C:","Users", "blbouf", "Sync", "Paper2"))

allHRUs_path <- file.path(".", "data", "LCC_HRU_files", "Dec2")

allHRUs_file <- list.files(allHRUs_path, full.names = TRUE)
allHRUs_file <- allHRUs_file[grepl(pattern = "HRUs_1923_2023_updated_lai_classes.csv", 
                                   allHRUs_file)]
allHRUs <- read.csv(allHRUs_file)
allHRUs_init <- allHRUs[allHRUs$Year %in% 1980:2023, ]

allHRUs_rvp <- allHRUs[allHRUs$Year == 1980, ] %>% select(-c("harv_yr", "fire_yr", "forest_class", 
                                                              "age", "Year", "ZONE"))
allHRUs_rvp$AQUIFER[is.na(allHRUs_rvp$AQUIFER)] <- "[NONE]"

allHRUs_rvp$TERRAIN[is.na(allHRUs_rvp$TERRAIN)] <- "[NONE]"

write.csv(allHRUs_rvp,
          file.path(allHRUs_path, "HRU_rvh_1980.csv"),
          row.names = FALSE)

output_dir <- allHRUs_path

#' Make Dynamic Forest Change files
#' @return a tibble of HRUs for based on assumed land cover for a specifed year
#' @param allHRUs hru data frame containing all years of sim (typical for recode_forest output)
#' @param scenario name of land cover scenario
#' @param dir directory of model (defaults to working directory)
#' @import dplyr
#' @import tidyr
#' @import readr
#' @export
#' 
#################
# for testing
################
d <- output_dir
scenario <- "Trapping_Dec2_baseline"

names(allHRUs) <-  c("ID", "AREA", "ELEVATION", "LATITUDE", "LONGITUDE", "BASIN_ID",
                      "LAND_USE_CLASS", "VEG_CLASS", "SOIL_PROFILE" , "AQUIFER_PROFILE","TERRAIN_CLASS" ,
                      "SLOPE", "ASPECT", "ZONE", "burn_year", "harv_year", "age", "Year", "forest_class")
names(allHRUs_init) <-  c("ID", "AREA", "ELEVATION", "LATITUDE", "LONGITUDE", "BASIN_ID",
                     "LAND_USE_CLASS", "VEG_CLASS", "SOIL_PROFILE" , "AQUIFER_PROFILE","TERRAIN_CLASS" ,
                     "SLOPE", "ASPECT", "ZONE", "burn_year", "harv_year", "age", "Year", "forest_class")


Dynamic_HRUs(allHRUs_init, scenario, d)

Dynamic_HRUs = function(allHRUs, scenario, d = getwd()){
  # Make Directory
  if(!dir.exists(file.path(d, 'LandCoverChange'))){ dir.create(file.path(d, 'LandCoverChange'), recursive = T) }
  folder_dir = file.path(d,'LandCoverChange', scenario)
  if(!dir.exists(folder_dir)){ dir.create(folder_dir, recursive = T) }
  
  # Create annual table of land use class
  annual_table = allHRUs %>% dplyr::select(ID, LAND_USE_CLASS, Year) %>% spread(Year, LAND_USE_CLASS)
  
  # Function to get all changes
  get_annual_changes = function(yr){
    x = which(colnames(annual_table) == yr)
    temp = annual_table[which(annual_table[,x] != annual_table[,x-1]),]
    
    # Make HRU Groups for land cover change
    make_aHRUGroup = function(LC, yr){
      c(paste0(':HRUGroup ',paste(LC, yr, sep = '_')),
        temp$ID[temp[,x] == LC] %>% paste0(., collapse = ','),
        paste0(':EndHRUGroup'),
        '#')
    }
    
    lapply(unique(temp[,x]), make_aHRUGroup, yr = yr) %>% unlist()
  }
  
  group_list = lapply(unique(allHRUs$Year)[-1], get_annual_changes) %>% unlist # don't use first year
  
  # Generate list of HRU Groups
  HGroups = gsub(':HRUGroup ', '', group_list[grepl('HRUGroup ', group_list)])
  hru_groups = c(paste(':DefineHRUGroups',paste0(HGroups, collapse = ',')))
  
  second_number <- function(x) stringr::str_extract_all(x, "\\d+")[[1]][2]
  
  # Land Use Change Calls
  make_change = function(x){
    
    if (grepl("D_", x) | grepl("M_", x)){
      c(paste(':LandUseChange', x, gsub('_[0-9]{4}', '', x), paste0(parse_number(x),'-01-01')),
        paste(':VegetationChange', x, gsub('_[0-9]{4}', '', x), paste0(parse_number(x),'-01-01'))
      )
    } else {
    
    c(paste(':LandUseChange', x, gsub('_[0-9]{4}', '', x), paste0(second_number(x),'-01-01')),
      paste(':VegetationChange', x, gsub('_[0-9]{4}', '', x), paste0(second_number(x),'-01-01'))
    )
    }
  }
  veg_change = lapply(HGroups, make_change) %>% unlist
  
  # Write files
  write_lines(group_list, file.path(folder_dir, "HRUGroups.txt"))
  write_lines(hru_groups, file.path(folder_dir, "HRUList.txt"))
  write_lines(veg_change, file.path(folder_dir, "RVPList.txt"))
  
}
