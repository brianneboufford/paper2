# ------------------------------------------------------------------------------
# make plot of disturbance scenarios 
# 
# Dec 8th
# ------------------------------------------------------------------------------

setwd("C:/Users/blbouf/Sync/Paper2")
# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

dist_shp_path <- file.path(".", "data", "disturbance_recovery_scenarios") %>% 
  list.files(pattern = "_0yrs",
             full.names = TRUE)

hru_all_forest_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp")
fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios")


hru_af <- st_read(hru_all_forest_path)
 
# percent forested calculation 
hru_forest_area <- sum(st_area(hru_af[!hru_af$LAND_US %in% c("WET_LAND", "ALPINE", "SHRUB"), ]))/(1000*1000)
percent_for <- hru_forest_area/(sum(st_area(hru_af)/(1000*1000)))

dist_shp_path <- dist_shp_path[!grepl("_40_", dist_shp_path)]
dist_shp_path <- dist_shp_path[grepl(".shp$", dist_shp_path)]

dist_shp <- dist_shp_path[1]

lapply(dist_shp_path, add_dist_polygons, 
       hru_af = hru_af, 
       fig_path = fig_path)

add_dist_polygons <- function(dist_shp, hru_af, fig_path){
  
  fig_filename <- basename(dist_shp) %>% str_replace(.,
                                                     pattern = ".shp$", 
                                                     replacement = ".png")
  fig_filename_full <- file.path(fig_path, fig_filename)
  dist_obj <- st_read(dist_shp) %>% 
    select(-c("p50", "w", "age"))
  dist_obj$VEG_CLA <- "DISTURBANCE"
  dist_ids <- dist_obj$ID
  
  hru_forested <- hru_af[!hru_af$ID %in% dist_ids, ]

  hru_sim <- rbind(hru_forested, dist_obj)
  
  veg_colors <- c(
    WET_LAND     = "#7DA9A3",  # muted teal
    SHRUB        = "#C7BFA5",  # muted tan
    ALPINE       = "#D9D9D9",  # light stone gray
    FOREST_IDF   = "#5B7F4A",  # muted green
    FOREST_MS    = "#45663A",  # darker muted green
    FOREST_ESSF  = "#2F4F2F",   # dark forest green
    DISTURBANCE    = "orange"
  )
  
  sim_plot <- ggplot(hru_sim) +
    geom_sf(aes(fill = VEG_CLA, color = VEG_CLA)) +  # keep color aes but hide borders
    scale_fill_manual(values = veg_colors) +
    scale_color_manual(values = veg_colors) +
    theme_void() +
    labs(fill = "Vegetation Class") +
    guides(
      color = "none",   # remove separate color legend
      fill = guide_legend(
        override.aes = list(color = NA)  # remove borders from legend swatches
      )
    )
  
  ggsave(filename = fig_filename_full,
         plot = sim_plot, 
         units = "in",
         width = 4, height = 3)
  
}
