# ------------------------------------------------------------------------------
# fit S-curves to ALS-derived forest parameters 
#
# date: February 12, 2026 
# update March 22, 2026 to use new LAI, height, and fcover
# updated April 14th to fix forest param figure
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(smoothr)
library(stringr)
library(ggplot2)
library(future)
library(future.apply)
library(FlexParamCurve)
library(nlme)
library(minpack.lm)
library(cowplot)

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# fcover, veg height
als_veg_params_path <- file.path(".", "data", "forest_params_by_age", "fc_height_2015_mar22.csv") # fcover and height # was feb 11

zone_path <- file.path(".", "data", "src", "ntems", "zone_key.csv")

fig_path <- file.path(".", "data", "figs", "forest_params")

# output path
result_path <- file.path(file.path(".", "data", "med_forest_params_curve_fitted")) 

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# BEC zone key 
zone_key <- read.csv(zone_path)
zone_key <- zone_key[c(1,2,5), ]

# read sampled als forest param data (height, fc)
als_veg <- read.csv(als_veg_params_path)
als_veg_data <- merge(als_veg, zone_key, by="id")
als_veg_data$als_grp_key <- paste0(als_veg_data$age, "_", als_veg_data$lai_grp)

# ------------------------------------------------------------------------------
# filter out data that has n < 100
# ------------------------------------------------------------------------------

als_summary_df <- als_veg_data %>% 
  group_by(als_grp_key) %>% 
  summarize(n = length(height),
            h_p10 = quantile(height, 0.1, na.rm = TRUE),
            h_p90 = quantile(height, 0.9, na.rm = TRUE),
            h_med = median(height, na.rm = TRUE),
            h_mean = mean(height, na.rm = TRUE),
            fc_p10 = quantile(fc, 0.1, na.rm = TRUE),
            fc_p90 = quantile(fc, 0.9, na.rm = TRUE),
            fc_med = median(fc, na.rm = TRUE),
            fc_mean = mean(fc, na.rm = TRUE),
            als_grp_key = first(als_grp_key)) #%>% 
 # subset(n > 100)

# remove the ones where RSE cannot be calculated because n < 2
keys_to_remove <- als_summary_df$als_grp_key[als_summary_df$n < 2]
als_veg_data <- als_veg_data[!als_veg_data$als_grp_key %in% keys_to_remove, ]

# get rse for each age 
sampled_forested_rse <- als_veg_data %>%
  group_by(als_grp_key) %>%
  slice_sample(n = 10000) %>%
  ungroup() %>%
  filter(age < 100) %>%
  group_by(als_grp_key) %>% 
  summarize(rse_fcover = AnglerCreelSurveySimulation::calculate_rse(na.omit(fc)),
            rse_height = AnglerCreelSurveySimulation::calculate_rse(na.omit(height)), 
            n_fc = length(fc),
            n_h = length(height))

# data to remove 
remove_fc_rse <- sampled_forested_rse$als_grp_key[sampled_forested_rse$rse_fcover > 0.1] # was 0.04
remove_h_rse <- sampled_forested_rse$als_grp_key[sampled_forested_rse$rse_height > 0.1]

# unique IDs to remove bc RSE > 0.04 for either height or fc 
remove_rse <- c(remove_fc_rse, remove_h_rse)

remove_n <- sampled_forested_rse$als_grp_key[sampled_forested_rse$n_fc < 100]

# ------------------------------------------------------------------------------
# change group label to string 
# ------------------------------------------------------------------------------

# als forst params 
als_veg_data$lai_grp[als_veg_data$lai_grp == 0] <- "g0"
als_veg_data$lai_grp[als_veg_data$lai_grp == 1] <- "g1"
als_veg_data <- als_veg_data[als_veg_data$frst_cl %in% c("ESSF", "IDF", "MS"), ]

# ------------------------------------------------------------------------------
# define age classes
# ------------------------------------------------------------------------------

# look up table for BEC zone recovery class distinction
bec_zoned_classes <- list(
  MS = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_ms", "R1_ms", "R2_ms", "R3_ms", "R4_ms", "R5_ms", 
                "R6_ms", "R7_ms", "R8_ms", "R9_ms", "M_ms")
  ),
  ESSF = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_essf", "R1_essf", "R2_essf", "R3_essf", "R4_essf", "R5_essf", 
                "R6_essf", "R7_essf", "R8_essf", "R9_essf","M_essf")
  ),
  IDF = data.frame(
    min_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    max_age = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
    class   = c("D_idf", "R1_idf", "R2_idf", "R3_idf", "R4_idf", "R5_idf", 
                "R6_idf", "R7_idf", "R8_idf", "R9_idf","M_idf")
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

# ------------------------------------------------------------------------------
# assign age class to ALS data
# ------------------------------------------------------------------------------

# assign age class
als_veg_data <- assign_ALS_ageclass(als_veg_data, elev_zoned_classes, "eco")

# filter out years > 100 bc low accuracy for age 
als_veg_data <- als_veg_data[als_veg_data$age <= 100, ]
als_veg_data1 <- als_veg_data[!als_veg_data$als_grp_key %in% remove_rse, ]

als_veg_data1 <- als_veg_data1 %>%
  dplyr::filter(!(age < 10 & height > 10)) %>%
  dplyr::filter(!(age < 10 & fc > 50)) %>% 
  dplyr::filter(!(age > 10 & age < 17 & height < 1 & fc < 10)) %>% 
  dplyr::filter(!(age > 10 & age < 17  & height > 15 & fc > 60))

g0_data <- als_veg_data1[als_veg_data1$lai_grp == "g0", ]
g1_data <- als_veg_data1[als_veg_data1$lai_grp == "g1", ]


test_g0 <- g0_data[g0_data$age < 5, ]
test_g0 <- test_g0[test_g0$age > 12, ] 

ggplot(test_g0, aes(x = fc, y = height)) + 
  geom_point(alpha = 0.1)

# ------------------------------------------------------------------------------
# g0 Forest Cover 
# ------------------------------------------------------------------------------

df1 <- g0_data

avg_fc_g0 <- df1 %>%
  group_by(age) %>%
  summarise(
    mean_fc = median(fc, na.rm = TRUE),
    Qhigh = quantile(fc, 0.9, na.rm = TRUE),
    Qlow = quantile(fc, 0.1, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc_g0 <- nls(
  mean_fc ~ SSlogis(age, Asym, xmid, scal),
  data = avg_fc_g0
)

age_seq_fc_g0 <- data.frame(
  age = seq(min(avg_fc_g0$age), max(avg_fc_g0$age), length.out = 101)
)

age_seq_fc_g0$fc_fit <- predict(fit_fc_g0, newdata = age_seq_fc_g0)

fc_g0_plot <- ggplot(avg_fc_g0, aes(age, mean_fc)) +
  #  geom_point(size = 2) +
  #  geom_line(color = "grey60") +
  geom_line(
    data = age_seq_fc_g0,
    aes(age, fc_fit),
    color = "black",
    linewidth = 1.2
  ) +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh),
    alpha = 0.25,
    colour = NA
  ) +
  labs(
    x = "Age",
    y = "Forest Cover (%)",
    title = ""
  ) +
  ylim(c(0, 100)) +
  theme_minimal(14)

summary(fit_fc_g0)

# ------------------------------------------------------------------------------
# g1 Forest Cover 
# ------------------------------------------------------------------------------

df2 <- g1_data

avg_fc_g1 <- df2 %>%
  group_by(age) %>%
  summarise(
    mean_fc = median(fc, na.rm = TRUE),
    Qhigh = quantile(fc, 0.9, na.rm = TRUE),
    Qlow = quantile(fc, 0.1, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc_g1 <- nls(
  mean_fc ~ SSlogis(age, Asym, xmid, scal),
  data = avg_fc_g1
)

age_seq_fc_g1 <- data.frame(
  age = seq(min(avg_fc_g1$age), max(avg_fc_g1$age), length.out = 101)
)

age_seq_fc_g1$fc_fit <- predict(fit_fc_g1, newdata = age_seq_fc_g1)

fc_g1_plot <- ggplot(avg_fc_g1, aes(age, mean_fc)) +
  #geom_point(size = 2) +
  #eom_line(color = "grey60") +
  geom_line(
    data = age_seq_fc_g0,
    aes(age, fc_fit),
    color = "black",
    linewidth = 1.2
  ) +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh),
    alpha = 0.25,
    colour = NA
  ) +
  labs(
    x = "Age",
    y = "Forest Cover (%)",
    title = ""
  ) +
  ylim(c(0, 100)) +
  theme_minimal(14)

summary(fit_fc_g1)

# ------------------------------------------------------------------------------
# g1 Height 
# ------------------------------------------------------------------------------

df3 <- g1_data
# df <- df[!df$age %in% c(30, 31, 32, 33, 34, 35, 36, 37, 39, 40), ]

avg_height_g1 <- df3 %>%
  group_by(age) %>%
  summarise(
    mean_h = median(height, na.rm = TRUE),
    Qhigh = quantile(height, 0.9, na.rm = TRUE),
    Qlow = quantile(height, 0.1, na.rm = TRUE),
    .groups = "drop"
  )

fit_h_g1 <- nlsLM(
  mean_h ~ L / (1 + exp(-k * (age - x0))),
  data = avg_height_g1,
  start = list(
    L  = max(avg_height_g1$mean_h),
    k  = 0.03,
    x0 = median(avg_height_g1$age)
  ),
  lower = c(L = 10, k = 0.005, x0 = min(avg_height_g1$age)),
  upper = c(L = 80, k = 0.2,   x0 = max(avg_height_g1$age)),
  control = nls.lm.control(maxiter = 500)
)

age_seq_h_g1 <- data.frame(
  age = seq(min(avg_height_g1$age), max(avg_height_g1$age), length.out = 101)
)

age_seq_h_g1$fc_fit <- predict(fit_h_g1, newdata = age_seq_h_g1)

h_g1_plot <- ggplot(avg_height_g1, aes(age, mean_h)) +
  # geom_point(size = 2) +
  # geom_line(color = "grey60") +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(
    data = age_seq_h_g1,
    aes(age, fc_fit),
    color = "black",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Height (m)",
    title = ""
  ) +
  ylim(c(0, 20)) +
  theme_minimal(14)

summary(fit_h_g1)

# ------------------------------------------------------------------------------
# g0 Height 
# ------------------------------------------------------------------------------

df4 <- g0_data

avg_height_g0 <- df4 %>%
  group_by(age) %>%
  summarise(
    mean_h = median(height, na.rm = TRUE),
    Qhigh  = quantile(height, 0.9, na.rm = TRUE),
    Qlow = quantile(height, 0.1, na.rm = TRUE),
    .groups = "drop"
  )

fit_h_g0 <- nls(
  mean_h ~ SSlogis(age, Asym, xmid, scal),
  data = avg_height_g0
)

fit_h_g0 <- nlsLM(
  mean_h ~ L / (1 + exp(-k * (age - x0))),
  data = avg_height_g0,
  start = list(
    L  = max(avg_height_g0$mean_h),
    k  = 0.03,
    x0 = median(avg_height_g0$age)
  ),
  lower = c(L = 10, k = 0.005, x0 = min(avg_height_g0$age)),
  upper = c(L = 80, k = 0.2,   x0 = max(avg_height_g0$age)),
  control = nls.lm.control(maxiter = 500)
)

age_seq_h_g0 <- data.frame(
  age = seq(min(avg_height_g0$age), max(avg_height_g0$age), length.out = 101)
)

age_seq_h_g0$fc_fit <- predict(fit_h_g0, newdata = age_seq_h_g0)

h_g0_plot <- ggplot(avg_height_g0, aes(age, mean_h)) +
  #geom_point(size = 2) +
  #geom_line(color = "grey60") +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(
    data = age_seq_h_g0,
    aes(age, fc_fit),
    color = "black",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Height (m)",
    title = ""
  ) +
  ylim(c(0, 20)) +
  theme_minimal(14)

summary(fit_h_g0)

# ------------------------------------------------------------------------------
# make plot 
# ------------------------------------------------------------------------------

forest_plots <- plot_grid(fc_g0_plot, fc_g1_plot, h_g0_plot, h_g1_plot, 
                          nrow = 2, 
                          rel_widths = c(1,1),
                          align  = "hv", 
                          labels = "auto" 
                          )

ggsave(forest_plots, 
       filename = file.path(fig_path, "fc_h_plots_mar22.png"),
       dpi = 600, 
       units = "in", 
       width = 7, 
       height = 7)

# ------------------------------------------------------------------------------
# make new plot 
# ------------------------------------------------------------------------------
# forest cover 
avg_fc_g0$prod <- "High"
avg_fc_g1$prod <- "Low"
age_seq_fc_g0$prod <- "High"
age_seq_fc_g1$prod <- "Low"

all_fc <- rbind(avg_fc_g0, avg_fc_g1)
all_fc <- all_fc[all_fc$age < 65, ]
all_fc_fit <- rbind(age_seq_fc_g0, age_seq_fc_g1)
all_fc_fit <- all_fc_fit[all_fc_fit$age < 65, ]

# height 
avg_height_g0$prod <- "High"
avg_height_g1$prod <- "Low"
age_seq_h_g0$prod <- "High"
age_seq_h_g1$prod <- "Low"

all_h <- rbind(avg_height_g0, avg_height_g1)
all_h <- all_h[all_h$age < 65, ]
all_h_fit <- rbind(age_seq_h_g0, age_seq_h_g1)
all_h_fit <- all_h_fit[all_h_fit$age < 65, ]

fc_plot_highlow <- ggplot(all_fc, aes(age, mean_fc, colour = prod)) +
  #geom_point(size = 2) +
  #eom_line(color = "grey60") +
  geom_line(
    data = all_fc_fit,
    aes(age, fc_fit, colour = prod),
    linewidth = 1.2
  ) +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh, fill = prod),
    alpha = 0.25,
    colour = NA
  ) +
  scale_color_manual(values = c(
    "High" = "#c51b7d",
    "Low"  = "#4d9221"
  )) +
  scale_fill_manual(values = c(
    "High" = "#c51b7d",
    "Low"  = "#4d9221"
  )) + 
  geom_vline(
    xintercept = seq(5, 55, by = 5),
    linetype = "dashed",
    color = "grey50"
  ) +  
  labs(
    x = "Age",
    y = "Forest Cover (%)",
    colour = "Productivity Region",
    fill = "Productivity Region",
    title = ""
  ) +
  scale_x_continuous(breaks = seq(0, max(all_fc$age , na.rm = TRUE), by = 5)) +
  ylim(c(0, 100)) +
  theme_minimal(12) + 
  theme(legend.position = "none")


h_plot_highlow <- ggplot(all_h, aes(age, mean_h, color = prod)) +
  # geom_point(size = 2) +
  # geom_line(color = "grey60") +
  geom_ribbon(
    aes(ymin = Qlow, ymax = Qhigh, fill = prod),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(
    data = all_h_fit,
    aes(age, fc_fit, colour = prod),
    linewidth = 1.2
  ) +
  scale_color_manual(values = c(
    "High" = "#c51b7d",
    "Low"  = "#4d9221"
  )) +
  scale_fill_manual(values = c(
    "High" = "#c51b7d",
    "Low"  = "#4d9221"
  )) + 
  geom_vline(
    xintercept = seq(5 , 55, by = 5),
    linetype = "dashed",
    color = "grey50"
  ) + 
  labs(
    x = "Age",
    y = "Height (m)",
    colour = "Productivity Region",
    fill = "Productivity Region",
    title = ""
  ) +
  scale_x_continuous(breaks = seq(0, max(all_h$age, na.rm = TRUE), by = 5)) + 
  ylim(c(0, 20)) +
  theme_minimal(12) + 
  theme(legend.position = "bottom")

h_plot_highlow
fc_plot_highlow

hilow_plots <- plot_grid(fc_plot_highlow, h_plot_highlow,
                          nrow = 2, 
                          align  = "hv"
)

ggsave(hilow_plots, 
       filename = file.path(fig_path, "fc_h_plots_april14.png"),
       dpi = 600, 
       units = "in", 
       width = 5, 
       height = 8)


# ------------------------------------------------------------------------------
# get summary stats for results section 
# ------------------------------------------------------------------------------
fit_fc_g0 
fit_fc_g1
fit_h_g0
fit_h_g1

age_seq_h_g1$compare_asym <- age_seq_h_g1/10.00
age_seq_h_g0$compare_asym <- age_seq_h_g0/11.7670

age_seq_fc_g1$compare_asym <- age_seq_fc_g1/69.433
age_seq_fc_g0$compare_asym <- age_seq_fc_g0/75.614

age_seq_fc_g0
age_seq_fc_g1

age_seq_h_g0
age_seq_h_g1
# ------------------------------------------------------------------------------
# join all fit data and turn into useful CSV format for raven input file 
# development
# ------------------------------------------------------------------------------

# rename columns and add lai_grp before join 
names(age_seq_fc_g0) <- c("age", "forest_cover")
names(age_seq_h_g0) <- c("age", "height")
g0_fitted_data <- merge(age_seq_fc_g0, age_seq_h_g0, by = "age")
g0_fitted_data$lai_grp <- "g0"

names(age_seq_fc_g1) <- c("age", "forest_cover")
names(age_seq_h_g1) <- c("age", "height")
g1_fitted_data <- merge(age_seq_fc_g1, age_seq_h_g1, by = "age")
g1_fitted_data$lai_grp <- "g1"

# merge all data into one 
h_fc_data <- rbind(g0_fitted_data, g1_fitted_data)

# pull function to assign class outside of wrapper function 
assign_age_lai_grp_class <- function(age, zone, zoned_classes) {
  class_df <- zoned_classes[[as.character(zone)]]
  if (is.null(class_df)) return(zone)
  idx <- which(age >= class_df$min_age & age < class_df$max_age)
  if (length(idx) == 0) return(zone)
  class_df$class[idx[1]]
}

# apply 
h_fc_data <- h_fc_data %>%
  rowwise() %>%
  mutate(age_class = assign_age_lai_grp_class(age, lai_grp, elev_zoned_classes)) %>%
  ungroup()

h_fc_grouped <- h_fc_data %>% 
  group_by(age_class) %>% 
  summarize(med_fc = median(forest_cover),
            med_h = median(height))

# -------
# write output 
# -------

write.csv(h_fc_grouped, 
          file.path(result_path, "height_fcover_med_Scurve_mar22.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# functions 
# ------------------------------------------------------------------------------

# function to get mode - built in R function doesn't compute the statistical mode 
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

assign_ALS_ageclass <- function(df, zoned_classes, class_type) {
  
  # df: data frame with columns: zone, age, winter, apr, ..., oct
  # zoned_classes: named list of data frames as described
  # class_type: either "eco" or "elev" to determine the df column to dictate the subzones for the age classification
  
  # Function to assign age class based on zone-specific breakpoints
  assign_class <- function(age, zone) {
    class_df <- zoned_classes[[as.character(zone)]]
    if (is.null(class_df)) return(zone)
    idx <- which(age >= class_df$min_age & age < class_df$max_age)
    if (length(idx) == 0) return(zone)
    class_df$class[idx[1]]
  }
  
  # Apply the function rowwise to assign age class
  if (class_type == "eco"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, frst_cl)) %>%
      ungroup()
    
  } else if(class_type == "elev"){
    df <- df %>%
      rowwise() %>%
      mutate(age_class = assign_class(age, lai_grp)) %>%
      ungroup()
  }
  
  return(df)
  
}

