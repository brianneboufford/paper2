# ------------------------------------------------------------------------------
# make plot of disturbance scenarios 
# 
# Dec 8th
# ------------------------------------------------------------------------------
library(patchwork)
library(ggpubr)
setwd("C:/Users/blbouf/Sync/Paper2")
# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------ 

dist_shp_path <- file.path(".", "data", "disturbance_recovery_scenarios") %>% 
  list.files(pattern = "_0yrs",
             full.names = TRUE)

hru_all_forest_path <- file.path(".", "data", "HRU_delineation", "new_HRUs", "TrappingCreek_HRUs_no_slp_asp_1125.shp")
aspect_path <- file.path(".", "data", "src", "lidar_derived", "aspect.tif")
slope_path <- file.path(".", "data", "src", "lidar_derived", "slope.tif")

fig_path <- file.path(".", "data", "figs", "disturbance_recovery_scenarios_Feb20")

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

plots <- lapply(dist_shp_path, make_dist_scen_plot,
       aspect_path, 
       slope_path, 
       hru_af, 
       fig_path) 

leg <- get_colour_legend(dist_shp, 
                         hru_af)

ggsave(leg, 
       filename = file.path(fig_path, "legend.png"),
       units = "in", 
       width = 8, 
       height = 1)

test <- plot_grid(
  plotlist = plots,
  nrow = 2,
  align = "hv",
  labels = "auto"
)

test2 <- plot_grid(test, 
                   leg, 
                   nrow = 2, 
                   rel_heights = c(10, 1))

ggsave(test2, 
       filename = file.path(fig_path, "simulation_areas_plots_legend.png"),
       units = "in",
       width = 10, 
       height = 10)

ggsave(test, 
       filename = file.path(fig_path, "simulation_areas_plots.png"),
       units = "in",
       width = 10, 
       height = 10)

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------
get_colour_legend <- function(dist_shp,
                              hru_af){
  # get figure filename
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
    DISTURBANCE    = "#543005"
  )
  
  veg_labs <- c(
    WET_LAND     = "Wetland",  # muted teal
    SHRUB        = "Shrub",  # muted tan
    ALPINE       = "Alpine",  # light stone gray
    FOREST_IDF   = "IDF",  # muted green
    FOREST_MS    = "MS",  # darker muted green
    FOREST_ESSF  = "ESSF",   # dark forest green
    DISTURBANCE    = "Disturbance"
  )
  
  
  sim_plot <- ggplot(hru_sim) +
    geom_sf(aes(fill = VEG_CLA, color = VEG_CLA)) +  # keep color aes but hide borders
    scale_fill_manual(values = veg_colors, labels = veg_labs) +
    scale_color_manual(values = veg_colors) +
    theme_void() +
    labs(fill = "Landcover Type") +
    guides(
      color = "none",   # remove separate color legend
      fill = guide_legend(
        override.aes = list(color = NA)  # remove borders from legend swatches
      )
    ) + 
    theme(legend.position = "bottom")
  
  leg <- get_legend(sim_plot)
  as_ggplot(leg)
  
}

make_dist_scen_plot <- function(dist_shp,
                                aspect_path, 
                                slope_path, 
                                hru_af, 
                                fig_path){
  
  # read topography layers 
  asp <- rast(aspect_path)
  slp <- rast(slope_path)
  
  # get figure filename
  fig_filename <- basename(dist_shp) %>% str_replace(.,
                                                     pattern = ".shp$", 
                                                     replacement = ".png")
  fig_filename_full <- file.path(fig_path, fig_filename)
  dist_obj <- st_read(dist_shp) %>% 
    select(-c("p50", "w", "age"))
  dist_obj$VEG_CLA <- "DISTURBANCE"
  dist_ids <- dist_obj$ID
  
  # clip slope raster to disturbed area only 
  disturbed_slope <- crop(slp %>% terra::project(crs(dist_obj)), dist_obj$geometry %>% vect() ) %>% 
    mask(., dist_obj$geometry %>% vect())
  
  # clip aspect to disturbed area only 
  disturbed_aspect <- crop(asp %>% terra::project(crs(dist_obj)), dist_obj$geometry %>% vect() ) %>% 
    mask(., dist_obj$geometry %>% vect())
  
  hru_forested <- hru_af[!hru_af$ID %in% dist_ids, ]
  
  hru_sim <- rbind(hru_forested, dist_obj)
  
  veg_colors <- c(
    WET_LAND     = "#7DA9A3",  # muted teal
    SHRUB        = "#C7BFA5",  # muted tan
    ALPINE       = "#D9D9D9",  # light stone gray
    FOREST_IDF   = "#5B7F4A",  # muted green
    FOREST_MS    = "#45663A",  # darker muted green
    FOREST_ESSF  = "#2F4F2F",   # dark forest green
    DISTURBANCE    = "#543005"
  )
  
  # ----------------------------------------------------------------
  # slope dataframe 
  # ----------------------------------------------------------------
  
  slp_df <- as.data.frame(disturbed_slope, na.rm = TRUE)
  colnames(slp_df) <- "slope"
  
  slp_df_binned <- slp_df %>%
    mutate(
      slope_bin = cut(
        slope,
        breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
        labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30+"),
        right = FALSE
      )
    )
  
  slp_summary_df <- slp_df_binned %>%
    count(slope_bin) %>%
    mutate(percent = n / sum(n) * 100)
  
  # slp_summary_df$slope_bin <- factor(
  #   slp_summary_df$slope_bin,
  #   levels = c("30+", "25-30", "20-25", "15-20", "10-15", "5-10", "0-5")
  # )
  
  # -------------------------------------------------------------------
  # aspect df
  # --------------------------------------------------------------------
  
  asp_df <- as.data.frame(disturbed_aspect, na.rm = TRUE)
  colnames(asp_df) <- "aspect"
  
  asp_df_class <- asp_df %>%
    mutate(
      aspect_class = case_when(
        aspect >= 337.5 | aspect < 22.5  ~ "N",
        aspect >= 22.5  & aspect < 67.5  ~ "NE",
        aspect >= 67.5  & aspect < 112.5 ~ "E",
        aspect >= 112.5 & aspect < 157.5 ~ "SE",
        aspect >= 157.5 & aspect < 202.5 ~ "S",
        aspect >= 202.5 & aspect < 247.5 ~ "SW",
        aspect >= 247.5 & aspect < 292.5 ~ "W",
        aspect >= 292.5 & aspect < 337.5 ~ "NW"
      )
    )
  asp_summary_df <- asp_df_class %>%
    count(aspect_class) %>%
    mutate(percent = n / sum(n) * 100)
  
  asp_summary_df$aspect_class <- factor(
    asp_summary_df$aspect_class,
    levels = c("NW", "W", "SW", "S", "SE", "E", "NE", "N")
  )
  
  asp_plot <- ggplot(asp_summary_df, aes(x = aspect_class, y = percent)) +
    geom_col(fill = "white", color = "black", position = "dodge", width = 0.9) +
    labs(
      x = "Aspect",
      y = "Percent area (%)",
      title = ""
    ) +
    theme_bw() + coord_flip() +
    theme(axis.text = element_text(size = 6), 
          axis.title = element_text(size = 8),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  slp_plot <- ggplot(slp_summary_df, aes(x = slope_bin, y = percent)) +
    geom_col(fill = "grey", colour = "black", position = "dodge", width = 0.9) +
    labs(
      x = "Slope (deg)",
      y = "Percent area (%)",
      title = ""
    ) +
    theme_bw() + coord_flip() +
    theme(axis.text = element_text(size = 6), 
          axis.title = element_text(size = 8),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
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
    ) + 
    theme(legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  asp_plot  <- asp_plot  +
    theme(axis.title.x = element_blank())
  
  slp_plot <- slp_plot +
    theme(axis.title.x = element_blank())
  
  x_label <- ggdraw() +
    draw_label(
      "Percent Area (%)",
      size = 8
    )
  
  final_plot <- plot_grid(
    sim_plot,
    plot_grid(asp_plot, slp_plot, ncol = 2),
    x_label,
    ncol = 1,
    rel_heights = c(4, 1, 0.15)
  )
  
}

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
    DISTURBANCE    = "#543005"
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
