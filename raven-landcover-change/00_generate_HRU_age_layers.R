# -----------------------------------------------------------------------------_
# creating annual HRU files to represent landcover change 
#
# copied from HRU_age_layers_2.R
# November 10th, 2025
# updated December 2nd to have the 7-class ecosystem scale classfification
# updated again Feb 19, 2026 to have 11 class Elevation based ecosystem scale classification 
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
setwd(file.path("C:","Users", "blbouf", "Sync", "Paper2"))

# hru shapefile 
hru_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp")

# ------------------------------------------------------------------------------
# read  
# ------------------------------------------------------------------------------
# all forested Catchment
hru_start <- st_read(hru_path)
hru_start$lai_grp <- "g0"
hru_start$lai_grp[hru_start$ELEVATI > 1560] <- "g1"
hru_start$lai_grp[hru_start$VEG_CLA == "WET_LAND"] <- "WET_LAND"
hru_start$lai_grp[hru_start$VEG_CLA == "ALPINE"] <- "ALPINE"
hru_start$lai_grp[hru_start$VEG_CLA == "SHRUB"] <- "SHRUB"
# ------------------------------------------------------------------------------
# start at some mature age 
# ------------------------------------------------------------------------------
first_year <- 1923 
end_year <- 2023

hru_start$age <- 50 # make age column and set everything to be mature 
hru_start$Year <- first_year # set start year to be first year

# initialize non-forested HRUs 
non_forested <- hru_start[hru_start$VEG_CLA %in% c("WET_LAND", "SHRUB","ALPINE"), ]
non_forested_ids <- non_forested$ID

# add forest class column 
hru_start$forest_class <- hru_start$ZONE
hru_start$forest_class[hru_start$ID %in% non_forested_ids] <- non_forested$VEG_CLA


# look up table for new 7-class ecosystem scale
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


# ----------------------------------------------------------
# hydrographs start at 1924
output_path <- file.path(".", "data", "LCC_HRU_files", "Feb19")

# make input to Dynamic_HRUs.R
make_annual_HRUs(hru_start = hru_start, 
                 start_year = first_year,
                 end_year = end_year, 
                 zoned_classes = elev_zoned_classes,
                 output_path = output_path)

#####################################################################
# functions 
#####################################################################
# function to update vegeation class 
update_vegetation_class <- function(df, zoned_classes) {
  
  # df : hru df for current year
  # zoned_classes : zoned_classes matrix 
  #
  # return updates df with new veg_class
  assign_class <- function(age, zone, existing_veg) {
    zone_class_df <- zoned_classes[[zone]]
    
    # If zone not in zoned_classes, keep original VEG_CLA
    if (is.null(zone_class_df) || is.na(age)) {
      return(existing_veg)
    }
    
    match_row <- which(age >= zone_class_df$min_age & age < zone_class_df$max_age)
    if (length(match_row) == 0) {
      return(existing_veg)
    }
    
    return(zone_class_df$class[match_row[1]])
  }
  
  # Apply row-wise using mapply, preserving original VEG_CLA if needed
  df$VEG_CLA <- mapply(assign_class, df$age, df$lai_grp, df$VEG_CLA)
  df$LAND_US <- df$VEG_CLA
  
  return(df)
}

########################################################################
# function to make each HRU year layer
make_annual_HRUs <- function(hru_start, start_year, end_year, zoned_classes, output_path){
  # @param hru_start : sf object for starting year
  # @param start_year : first year 
  # @param end_year : last year for HRU change 
  # @param output_path : path where output csv will be written 
  
  # @return 
  # csv with HRU details for each year 
  
  # make current year start - keep geometry 
  hru_current <- hru_start 
  
  # drop geometry 
  hru_start <- st_drop_geometry(hru_start)
  
  
  for (i in 1:(end_year - start_year)){
    
    print(i)
    
    # first current year is 1981 and goes to end year (2023)
    current_year <- start_year + i 
    hru_current$Year <- current_year
    
    # add 1 to the age 
    hru_current$age <- hru_current$age + 1
    
    # turn age to 0 if harvested this year
    disturb_this_year <- which(
      mapply(function(harv_yr, Year, fire_yr) {
        harv_yr == Year || Year %in% fire_yr}, 
        hru_current$harv_yr, hru_current$Year, hru_current$fire_yr)
    )
    
    hru_current$age[disturb_this_year] <- 0
    
    # update veg class
    hru_current <- update_vegetation_class(hru_current, zoned_classes)
    
    if (!file.exists(file.path(output_path, paste0("HRUs_", start_year, "_", end_year)))){
      dir.create(file.path(output_path, paste0("HRUs_", start_year, "_", end_year)))
    }
    
    # remove burn year so it will write properly 
    hru_write_current <- hru_current %>% select(-c("fire_yr"))
    st_write(hru_write_current, 
             file.path(output_path, paste0("HRUs_", start_year, "_", end_year), paste0("HRU_", current_year, ".shp")),
             append = FALSE)
    
    hru_current_bind <- hru_current %>% st_drop_geometry()
    # bind rows 
    hru_start <- bind_rows(hru_start, hru_current_bind)
    
  }
  
  # remove lai_grp so it will write properly
  hru_start_write <- hru_start %>% select(-c("lai_grp")) 
  write.csv(hru_start_write, 
            file.path(output_path, paste0("HRUs_", start_year, "_", end_year, "_updated_lai_classes.csv")),
            row.names = FALSE)
  
  # remove other columns 
  hru_write <- hru_start_write %>%
    select(-c("harv_yr", "forest_class", 
              "age", "Year"))
  
  write.csv(hru_write, 
            file.path(output_path, paste0("HRUs_", start_year, "_", end_year, "_for_rvh.csv")),
            row.names = FALSE)
  
  return(hru_start)
}
