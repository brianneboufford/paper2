# ------------------------------------------------------------------------------
# 
# make raven HRU files for each disuturbance and recovery scenario
#
# Dec 3rd, 2025
# updated Feb 19, 2026
# Brianne Boufford
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# libraries 

library(terra)
library(sf)

# ------------------------------------------------------------------------------
# paths 

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)


# output path
outpath <- file.path(".", "data", "rvp_rvh_data", "disturbance_recovery_scenarios_Feb19")

if (!file.exists(outpath)){
  dir.create(outpath)
}

hru_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp")
af_rvh_file <- file.path(".", "data", "rvp_rvh_data", "all_forested_Feb19", "Trapping_all_forest.rvh")
disturbance_scen_path <- file.path(".", "data", "disturbance_recovery_scenarios")

# ------------------------------------------------------------------------------
# data 
hru <- st_read(hru_path)
af_rvh <- rvn_rvh_read(af_rvh_file)
dist_files <- list.files(disturbance_scen_path,
                         pattern = ".shp$"
                         )
# age classes 
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

lapply(dist_files, 
       make_rvi, 
       disturbance_scen_path,
       elev_zoned_classes,
       af_rvh, 
       outpath)


make_rvi <- function(dist_file, 
                     disturbance_scen_path, 
                     zoned_classes,
                     af_rvh, 
                     outpath){
  
  dist_i <- st_read(file.path(disturbance_scen_path, dist_file))
  
  dist_filename_new <- stringr::str_replace(dist_file, ".shp$", ".rvh")
  
  # add new lai_grp column to dist_i 
  dist_i$lai_grp <- "g0"
  dist_i$lai_grp[dist_i$ELEVATI >= 1560] <- "g1"
  
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
  dist_age_class <- mapply(assign_class, dist_i$age, dist_i$lai_grp, dist_i$VEG_CLA) %>%
    cbind(dist_i$ID, .) %>%
    as.data.frame()
  
  names(dist_age_class) <- c("ID", "Vegetation") 
  dist_age_class$ID <- as.integer(dist_age_class$ID)
  
  hru_table <- af_rvh$HRUtable
  
  df_out <- hru_table %>%
    left_join(dist_age_class, by = "ID", suffix = c("", "_new")) %>%
    mutate(Vegetation = coalesce(Vegetation_new, Vegetation)) %>%
    select(-Vegetation_new)
  
  df_out$LandUse <- df_out$Vegetation
  
  # write as new rvh file
  rvn_rvh_write(file.path(outpath, paste0("Trapping", dist_filename_new)), 
                SBtable = af_rvh$SBtable,
                HRUtable = df_out)
  
  return(paste0("done", dist_filename_new))
}
  
  

