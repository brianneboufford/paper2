
setwd("C:/Users/blbouf/Sync/Paper2")

# get info about Trapping for climate NA
dem <- rast(file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek", "data", "lidar", "output", "dem", "dem.tif")) 

c <- st_read(file.path("C:", "Users", "blbouf", "Sync", "TrappingCreek", "data", "products", "TC_centroid.shp")) %>% 
  st_transform(r, crs = crs(dem))

shp_latlon <- st_transform(c, crs = 4326)

elev_c <- terra::extract(dem, c)

clim_path <- file.path("F:", "climate_1991_2020")
clim_files <- list.files(clim_path, full.names = TRUE)
clim_names <- list.files(clim_path) %>% 
  str_replace(., ".tif", "")

clim_stack <- rast(clim_files)
names(clim_stack) <- clim_names

catchment_path <- file.path(".", "data", "src", "TC_catchment", "catchments.shp")
catchment <- st_read(catchment_path) %>% 
  st_transform(., crs(clim_stack))

clim_tc <- crop(clim_stack, catchment$geometry) %>% 
  mask(., vect(catchment$geometry))
clim_averages <- terra::extract(clim_stack, 
                                vect(catchment$geometry), 
                                fun = mean)
