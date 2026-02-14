# ------------------------------------------------------------------------------
# script to make daymet location shapefile
#
# february 3rd, 2026
# ------------------------------------------------------------------------------
 
library(sf)
library(dplyr)

setwd(file.path("C:", "Users", "blbouf", "Sync", "Paper2"))
output_path <- file.path("data", "daymet")

gauges_df <- tibble::tibble(
  gauge = paste0("File", 1:9),
  latitude = c(
    49.5804335092447,
    49.5942996362676,
    49.6220318903135,
    49.6393645490921,
    49.6555416972855,
    49.6555416972855,
    49.6902070148428,
    49.6890515042576,
    49.7491380546903
  ),
  longitude = c(
    -119.036113714223,
    -119.017625544859,
    -118.987582269643,
    -118.967938589694,
    -119.003759417836,
    -118.920562655699,
    -119.011847991933,
    -118.912474081602,
    -118.885897338141
  ),
  elevation_m = c(
    959, 996, 1027, 1092, 1204, 1192, 1481, 1408, 1662
  )
)

# Convert to sf point object
gauges_sf <- st_as_sf(
  gauges_df,
  coords = c("longitude", "latitude"),
  crs = 4326   # WGS84
)

st_write(gauges_sf, 
         file.path(output_path, "daymet_locations_1_9.shp"),
         append = FALSE)
