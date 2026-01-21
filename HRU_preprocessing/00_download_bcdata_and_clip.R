# -----------------------------------------------------------------------------

# read in and clip 2020 landcover 
# download and clip cutblocks and fire perimeter layers from dataBC

# october 22, 2025

# -----------------------------------------------------------------------------

# libraries 
library(terra) 
library(dplyr) 
library(bcdata)

can_lc <- rast(file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data", "src", 
                        "2020_landcover", "landcover-2020-classification.tif"))

tc <- st_read(file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data", "src",
                        "TC_catchment", "catchments.shp")) %>% st_transform(crs(can_lc)) %>%
  vect()
  

tc_lc <- crop(can_lc, tc, mask = TRUE, 
              filename = file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                                   "data", "src","2020_landcover", 
                                   "landcover-2020-TC.tif"),
              overwrite = TRUE)

# ------------------------------------------------------------------------------
# download BC fire data 
# ------------------------------------------------------------------------------

bcdc_search("fire perimeter") # id 22c7cb44-1463-48f7-8e47-88857f207702

fire_perim <- bcdc_get_data("22c7cb44-1463-48f7-8e47-88857f207702")

st_write(fire_perim, 
         file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                   "data", "src","fire_perimeter", 
                   "fire_perim_BC.shp"))

fire_perim_read <- st_read(file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                                     "data", "src","fire_perimeter", 
                                     "fire_perim_BC.shp"))

tc_firecrs <- st_read(file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data", "src",
                        "TC_catchment", "catchments.shp")) %>% 
  st_transform(crs=crs(fire_perim_read)) 

fire_perim_crop <- st_intersection(fire_perim_read, tc_firecrs)

st_write(fire_perim_crop,
         file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                              "data", "src","fire_perimeter", 
                              "fire_perim_tc.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# download BC cutblock data 
# ------------------------------------------------------------------------------

bcdc_search("cutblocks") # id b1b647a6-f271-42e0-9cd0-89ec24bce9f7

cb <- bcdc_get_data("b1b647a6-f271-42e0-9cd0-89ec24bce9f7")

cb_new <- cb %>% select(c("HARVEST_START_YEAR_CALENDAR", "HARVEST_END_DATE",
                           "AREA_HA", "geometry"))
names(cb_new) <- c("harv_st", "harv_end", 
                   "area_ha", "geometry")

st_write(cb_new, 
         file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                   "data", "src","consolidated_cutblocks", 
                   "consolidated_cutblocks_BC_2.shp"),
         append = FALSE)

cb_read <- st_read(file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                                     "data", "src","consolidated_cutblocks", 
                             "consolidated_cutblocks_BC_2.shp"))

tc_cbcrs <- st_read(file.path("C:", "Users", "blbouf", "Sync", "Paper2", "data", "src",
                                "TC_catchment", "catchments.shp")) %>% 
  st_transform(crs=crs(cb_read)) 

cb_crop <- st_intersection(cb_read, tc_cbcrs)

st_write(cb_crop,
         file.path("C:", "Users", "blbouf", "Sync", "Paper2", 
                   "data", "src","consolidated_cutblocks", 
                   "consolidated_cutblocks_tc.shp"),
         append = FALSE)

