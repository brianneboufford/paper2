# ------------------------------------------------------------------------------
# Using Trapping_paper2_baseline as a template to generate new baseline RVH file 
# 
# November 26th, 2025
# last used/ edited Feb 19
#-------------------------------------------------------------------------------

# packages 
library(RavenR)

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

baseline_raven_path <- file.path(".", "raven-runs", "Baseline3")

hru_path <- file.path(".", "data", "LCC_HRU_files", "Feb19", "HRUs_1923_2023", "HRU_1980.shp")

hrus <- st_read(hru_path)
n_hrus <- length(hrus$ID)


bl_hru <- rvn_rvh_blankHRUdf(nHRUs = n_hrus)
bl_hru$Area <- hrus$AREA
bl_hru$Elevation <- hrus$ELEVATI
bl_hru$Latitude <- hrus$LATITUD
bl_hru$Longitude <- hrus$LONGITU
bl_hru$LandUse <- hrus$LAND_US
bl_hru$Vegetation <- hrus$LAND_US
bl_hru$Vegetation <- hrus$VEG_CLA
bl_hru$SoilProfile <- hrus$SOIL_PR
bl_hru$Slope <- hrus$SLOPE
bl_hru$Aspect <- hrus$ASPECT
bl_hru$fire_yr <- hrus$fire_yr
bl_hru$harv_yr <- hrus$harv_yr
bl_hru$ZONE <- hrus$ZONE

bl_sb <- rvn_rvh_blankSBdf(nSubBasins = 1)
bl_sb$Name <- "TRAPPING_CREEK_NEAR_THE_MOUTH"
bl_sb$Profile <- "08NN019"
bl_sb$ReachLength <- 24.20766
bl_sb$Gauged <- 1

rvn_rvh_write(filename = file.path(".", "raven-runs", "Baseline3", "Trapping_HRU_baseline.rvh"),
              HRUtable = bl_hru,
              SBtable = bl_sb)

# ------------------------------------------------------------------------------
# RVI files 
# ------------------------------------------------------------------------------
rvi_test <- rvn_rvi_read(file.path(baseline_raven_path, "Trapping.rvi"))
rvn_rvi_getparams(rvi_test)
conn <- rvn_rvi_connections(rvi_test)

rvn_rvi_process_diagrammer(rvi_conn = conn)
rvn_rvi_process_ggplot(rvi_conn = conn)

rvn_rvi_commandupdate(filename = file.path(baseline_raven_path, "Trapping.rvi"),
                      command = ":StartDate",
                      value = "2022-10-01 00:00:00",
                      outputfile = file.path(baseline_raven_path, "trapping_test.rvi"))

# ------------------------------------------------------------------------------
# RVP file 
#-------------------------------------------------------------------------------
