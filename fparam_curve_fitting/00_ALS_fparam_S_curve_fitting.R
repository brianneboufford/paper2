# ------------------------------------------------------------------------------
# fit S-curves to ALS-derived forest parameters 
#
# date: February 12, 2026 
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

setwd("C:/Users/blbouf/Sync/Paper2")

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# fcover, veg height
als_veg_params_path <- file.path(".", "data", "forest_params_by_age", "fc_height_2015_feb11.csv") # fcover and height

zone_path <- file.path(".", "data", "src", "ntems", "zone_key.csv")

# output path
result_path <- file.path(file.path(".", "data", "median_forest_params_LAI_by_age")) 

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
remove_fc_rse <- sampled_forested_rse$als_grp_key[sampled_forested_rse$rse_fcover > 0.04]
remove_h_rse <- sampled_forested_rse$als_grp_key[sampled_forested_rse$rse_height > 0.04]

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

g0_data <- als_veg_data[als_veg_data$lai_grp == "g0", ]
g1_data <- als_veg_data[als_veg_data$lai_grp == "g1", ]

# ------------------------------------------------------------------------------
# g0 Forest Cover 
# ------------------------------------------------------------------------------

df1 <- g0_data
df1 <- df1[!df1$als_grp_key %in% remove_rse, ]

avg_fc <- df1 %>%
  group_by(age) %>%
  summarise(
    mean_fc = mean(fc, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc <- nls(
  mean_fc ~ SSlogis(age, Asym, xmid, scal),
  data = avg_fc
)

age_seq <- data.frame(
  age = seq(min(avg_fc$age), max(avg_fc$age), length.out = 300)
)

age_seq$fc_fit <- predict(fit_fc, newdata = age_seq)

ggplot(avg_fc, aes(age, mean_fc)) +
  geom_point(size = 2) +
  geom_line(color = "grey60") +
  geom_line(
    data = age_seq,
    aes(age, fc_fit),
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Mean g0 FC",
    title = "Logistic S-curve fit: FC vs Age"
  ) +
  theme_minimal()

summary(fit_fc)

# ------------------------------------------------------------------------------
# g1 Forest Cover 
# ------------------------------------------------------------------------------

df2 <- g1_data
df2 <- df2[!df2$als_grp_key %in% remove_rse, ]

avg_fc <- df2 %>%
  group_by(age) %>%
  summarise(
    mean_fc = median(fc, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc <- nls(
  mean_fc ~ SSlogis(age, Asym, xmid, scal),
  data = avg_fc
)

age_seq <- data.frame(
  age = seq(min(avg_fc$age), max(avg_fc$age), length.out = 300)
)

age_seq$fc_fit <- predict(fit_fc, newdata = age_seq)

ggplot(avg_fc, aes(age, mean_fc)) +
  geom_point(size = 2) +
  geom_line(color = "grey60") +
  geom_line(
    data = age_seq,
    aes(age, fc_fit),
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Mean FC g1",
    title = "Logistic S-curve fit: FC vs Age"
  ) +
  theme_minimal()

summary(fit_fc)

# ------------------------------------------------------------------------------
# g1 Height 
# ------------------------------------------------------------------------------

df3 <- g1_data
df3 <- df3[!df3$als_grp_key %in% remove_rse, ]
# df <- df[!df$age %in% c(30, 31, 32, 33, 34, 35, 36, 37, 39, 40), ]

avg_height <- df3 %>%
  group_by(age) %>%
  summarise(
    mean_h = mean(height, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc <- nls(
  mean_h ~ SSlogis(age, Asym, xmid, scal),
  data = avg_height
)

fit_fc <- nls(
  mean_h ~ L / (1 + exp(-k * (age - x0))),
  data = avg_height,
  start = list(
    L  = max(avg_height$mean_h, na.rm = TRUE),
    k  = 0.05,
    x0 = median(avg_height$age)
  )
)

fit_h <- nlsLM(
  mean_h ~ L / (1 + exp(-k * (age - x0))),
  data = avg_height,
  start = list(
    L  = max(avg_height$mean_h),
    k  = 0.03,
    x0 = median(avg_height$age)
  ),
  lower = c(L = 10, k = 0.005, x0 = min(avg_height$age)),
  upper = c(L = 80, k = 0.2,   x0 = max(avg_height$age)),
  control = nls.lm.control(maxiter = 500)
)

age_seq <- data.frame(
  age = seq(min(avg_fc$age), max(avg_fc$age), length.out = 300)
)

age_seq$fc_fit <- predict(fit_h, newdata = age_seq)

ggplot(avg_height, aes(age, mean_h)) +
  geom_point(size = 2) +
  geom_line(color = "grey60") +
  geom_line(
    data = age_seq,
    aes(age, fc_fit),
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Mean Height",
    title = "Logistic S-curve fit: FC vs Age"
  ) +
  theme_minimal()

summary(fit_fc)

# ------------------------------------------------------------------------------
# g1 Height 
# ------------------------------------------------------------------------------

df4 <- g0_data
df4 <- df4[!df4$als_grp_key %in% remove_rse, ]

avg_height <- df4 %>%
  group_by(age) %>%
  summarise(
    mean_h = mean(height, na.rm = TRUE),
    .groups = "drop"
  )

fit_fc <- nls(
  mean_h ~ SSlogis(age, Asym, xmid, scal),
  data = avg_height
)

age_seq <- data.frame(
  age = seq(min(avg_fc$age), max(avg_fc$age), length.out = 300)
)

age_seq$fc_fit <- predict(fit_h, newdata = age_seq)

ggplot(avg_height, aes(age, mean_h)) +
  geom_point(size = 2) +
  geom_line(color = "grey60") +
  geom_line(
    data = age_seq,
    aes(age, fc_fit),
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Age",
    y = "Mean Height",
    title = "Logistic S-curve fit: FC vs Age"
  ) +
  theme_minimal()

summary(fit_fc)

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
