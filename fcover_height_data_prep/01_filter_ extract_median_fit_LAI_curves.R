# ------------------------------------------------------------------------------
# filter, extract median, and fit quandtratic function to seasonal LAI curves 
# 
# february 13, 2026
# last updated february 19, 2026 to have more filtering on LAI values - removing physically unrealistic values
#
# ------------------------------------------------------------------------------

# library 
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)

# ------------------------------------------------------------------------------
# paths 
# ------------------------------------------------------------------------------

# fcover, veg height, and LAI
veg_params_path <- file.path(".", "data", "forest_params_by_age", 
                             "sampled_veg_params_byLAI_GRP_2015_2021_feb19.csv") # was feb 13 # can also try with all of the data instead

zone_path <- file.path(".", "data", "src", "ntems", "zone_key.csv")

# output path
result_path <- file.path(".", "data", "median_forest_params_LAI_by_age")

fig_path <- file.path(".", "data", "figs", "forest_params")
# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# BEC zone key 
zone_key <- read.csv(zone_path)
zone_key <- zone_key[c(1,2,5), ]

# read sampled forst param data (lai, height, fc)
veg <- read.csv(veg_params_path)
veg_data <- merge(veg, zone_key, by="id")

# ------------------------------------------------------------------------------
# change group label to string 
# ------------------------------------------------------------------------------
# forest params 
veg_data$lai_grp[veg_data$lai_grp == 0] <- "g0"
veg_data$lai_grp[veg_data$lai_grp == 1] <- "g1"

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
# data cleaning 
# ------------------------------------------------------------------------------

# assign class 

# assign age class
veg_data <- assign_ageclass(veg_data, elev_zoned_classes, "elev")

# filter out years > 100 bc low accuracy for age 
veg_data <- veg_data[veg_data$age <= 100, ]

veg_data_summary <- veg_data %>% 
  group_by(age_class) %>%
  summarize(q90_july = quantile(july_lai, probs = 0.9, na.rm = TRUE),
            q80_july = quantile(july_lai, probs = 0.8, na.rm = TRUE), 
            q60_july = quantile(july_lai, probs = 0.6, na.rm = TRUE))


veg_data_merge <- merge(veg_data, veg_data_summary, by="age_class")

veg_data_test <- veg_data_merge 

veg_data_test <- veg_data_test %>% 
  mutate(nov_i = if_else(nov_i > q60_july, NA, nov_i)) %>%
  mutate(dec_i = if_else(dec_i > q60_july, NA, dec_i)) %>%
  mutate(jan_lai = if_else(jan_lai > q60_july, NA, jan_lai)) %>%
  mutate(feb_lai = if_else(feb_lai > q60_july, NA, feb_lai)) %>%
  mutate(mar_lai = if_else(mar_lai > q60_july, NA, mar_lai)) %>%
  mutate(april_lai = if_else(april_lai > q60_july, NA, april_lai)) %>%
  mutate(may_lai = if_else(may_lai > q80_july, NA, may_lai))

veg_data_test1 <- veg_data_test %>%
  
  # age between 7 and 25 i.e. not disturbed from 2014-2021 so chm is indicative of height
  # if chm > 20 rm 
  dplyr::filter(!(age > 7 & age < 15 & X2014_height > 15)) %>% 
  dplyr::filter(!(age >= 15  & age < 25 & X2014_height > 16)) %>% 
  # height > 10m and less than 7 yrs old 
  dplyr::filter(!(age < 7 & X2014_height > 10)) 
  
  # dplyr::filter(!(age > 7 & age < 25 & jan_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & feb_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & mar_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & april_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & may_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & june_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & july_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & aug_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & sept_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & oct_lai > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & nov_i > 3)) %>% 
  # dplyr::filter(!(age > 7 & age < 25 & dec_i > 3)) %>%
  #
  
# rm entries only not whole rows where LAI is wrong 
veg_data_test1$jan_lai[veg_data_test1$age < 10 & veg_data_test1$jan_lai > 3] <- NA
veg_data_test1$feb_lai[veg_data_test1$age < 10 & veg_data_test1$feb_lai > 3] <- NA
veg_data_test1$mar_lai[veg_data_test1$age < 10 & veg_data_test1$mar_lai > 3] <- NA
veg_data_test1$april_lai[veg_data_test1$age < 10 & veg_data_test1$april_lai > 3] <- NA
veg_data_test1$may_lai[veg_data_test1$age < 10 & veg_data_test1$may_lai > 3] <- NA
veg_data_test1$june_lai[veg_data_test1$age < 10 & veg_data_test1$june_lai > 3] <- NA
veg_data_test1$july_lai[veg_data_test1$age < 10 & veg_data_test1$july_lai > 3] <- NA
veg_data_test1$aug_lai[veg_data_test1$age < 10 & veg_data_test1$aug_lai > 3] <- NA
veg_data_test1$sept_lai[veg_data_test1$age < 10 & veg_data_test1$sept_lai > 3] <- NA
veg_data_test1$oct_lai[veg_data_test1$age < 10 & veg_data_test1$oct_lai > 3] <- NA
veg_data_test1$nov_i[veg_data_test1$age < 10 & veg_data_test1$nov_i > 3] <- NA
veg_data_test1$dec_i[veg_data_test1$age < 10 & veg_data_test1$dec_i > 3] <- NA


ggplot(data = veg_data_test1[veg_data_test1$age_class == "R4_g0" & veg_data_test1$age == 20, ], aes(x = X2014_height, y = july_lai)) + 
  geom_point(alpha = 0.1)

ggplot(data = veg_data_test1[veg_data_test1$age_class == "R4_g1" & veg_data_test1$age == 20, ], aes(x = X2014_height, y = july_lai)) + 
  geom_point(alpha = 0.1)

# -------------------------------------------------------------------

# ------------------------------------------------------------------------------
# median calculation 
# ------------------------------------------------------------------------------
reclass_lai <- veg_data_test1  %>% 
  group_by(age_class) %>%
  summarize(#med_apr_lai = round(median(april_lai, na.rm = TRUE),2),
            med_may_lai = round(median(may_lai, na.rm = TRUE),2),
            med_june_lai = round(median(june_lai, na.rm = TRUE),2),
            med_july_lai = round(median(july_lai, na.rm = TRUE),2),
            med_aug_lai = round(median(aug_lai, na.rm = TRUE),2),
            med_sept_lai = round(median(sept_lai, na.rm = TRUE),2),
            med_oct_lai = round(median(oct_lai, na.rm = TRUE),2))
         #   med_nov_lai = round(median(nov_i, na.rm = TRUE),2))

winter_meds <- veg_data_test1 %>%
  pivot_longer(
    cols = c(nov_i, dec_i, jan_lai, feb_lai, mar_lai, april_lai),
    names_to = "month",
    values_to = "lai"
  ) %>%
  group_by(age_class) %>%
  summarise(
    med_winter_lai = round(median(lai, na.rm = TRUE), 2),
    .groups = "drop"
  )

reclass_lai2 <- merge(reclass_lai, winter_meds, by="age_class")

reclass_lai2 <- reclass_lai2 %>% 
  mutate(med_winter_lai = if_else(med_winter_lai > med_may_lai, med_may_lai, med_winter_lai))

# clean up some erouneously high april, and winter values or Missing values 
# if may is missing, use June 
# if may is > jul make may = october 
# if april is > may then make april = may 
# if winter is > oct then make make winter october 
# reclass_lai <- reclass_lai %>%
#   mutate(med_may_lai = if_else(is.na(med_may_lai), med_june_lai, med_may_lai)) %>%
#   
#   mutate(med_may_lai = if_else(is.na(med_may_lai) | med_may_lai > med_july_lai, med_oct_lai, med_may_lai)) %>%
#   
#   mutate(med_apr_lai = if_else(is.na(med_apr_lai) | med_apr_lai > med_may_lai, med_may_lai, med_apr_lai)) %>%
#   
#   mutate(med_winter_lai = if_else(is.na(med_winter_lai) | med_winter_lai > med_oct_lai, med_may_lai, med_winter_lai))

# ------------------------------------------------------------------------------
# plot median 
# ------------------------------------------------------------------------------
lai <- reclass_lai2
classes <- lai$age_class

lai$lai_grp <- NA
lai$lai_grp[grepl("g0", lai$age_class)] <- "g0"
lai$lai_grp[grepl("g1", lai$age_class)] <- "g1"

lai_df_long <- lai %>%
  pivot_longer(
    cols = ends_with("lai"),
    names_to = "month",
    values_to = "lai"
  )

# Clean the month names
lai_df_long <- lai_df_long %>%
  mutate(month_clean = stringr::str_replace(month, "_lai", "")) %>%
  mutate(month_clean = stringr::str_replace(month_clean, "med_", "")) %>%
  filter(month_clean %in% c("may", "june", "july", "aug", "sept", "oct", "winter"))

# Ensure correct order of months
ordered_months <- c("may", "june", "july", "aug", "sept", "oct", "winter")
lai_df_long$month_clean <- factor(lai_df_long$month_clean, levels = ordered_months, 
                                  labels = c("05", "06", "07", "08", "09", "10", "winter"))


mean_lai_from_recovery_classes <- ggplot(lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                                          color = age_class)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  geom_line(data = lai_df_long, aes(x = month_clean, y = lai, group = age_class, 
                                    color = age_class), linewidth = 1) +  
  facet_wrap(~lai_grp, ncol = 2) +
  labs(x = "Month", 
       y = "Median LAI", 
       color = "Recovery Class", 
       fill = "Recovery Class") +   
  theme_minimal(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12),
        legend.position = "right")

# ------------------------------------------------------------------------------
# fit quadratic curve to median 
# ------------------------------------------------------------------------------
veg_df_long1 <- lai_df_long %>% 
  mutate(
    month_num = case_when(
      month == "med_winter_lai"  ~ 7, 
      month == "med_may_lai" ~ 1, 
      month == "med_june_lai" ~ 2,
      month == "med_july_lai" ~ 3,
      month == "med_aug_lai" ~ 4,
      month == "med_sept_lai" ~ 5, 
      month == "med_oct_lai" ~ 6,
    )
  )

lai_models <- veg_df_long1 %>% 
  group_by(age_class, lai_grp) %>%
  nest() %>% 
  mutate(
    model = map(data, ~lm(lai ~poly(month_num, 2, raw = TRUE), data = .x)), 
    model_summary = map(model, tidy)
  )

lai_coeffs <- lai_models %>% 
  unnest(model_summary) %>%
  dplyr::select(lai_grp, age_class, term, estimate, std.error, p.value)

pred_data <- lai_models %>% 
  mutate(
    newdata = map(
      data, 
      ~ data.frame(month_num = seq(min(.x$month_num), max(.x$month_num), length.out = 100)),
    ), 
    preds = map2(model, newdata, predict)
  ) %>% 
  unnest(c(newdata, preds))

pred_df <- veg_df_long1 %>%
  group_by(lai_grp, age_class) %>%
  group_modify(~ {
    fit <- lm(lai ~ poly(month_num, 2, raw = TRUE), data = .x)
    
    new_months <- data.frame(
      month_num = seq(min(.x$month_num), max(.x$month_num), length.out = 100)
    )
    
    new_months$lai_fit <- predict(fit, newdata = new_months)
    new_months
  })

# ------------------------------------------------------------------------------
# plot 
# ------------------------------------------------------------------------------

env_df <- veg_df_long1 %>%
  group_by(lai_grp, age_class, month_num) %>% 
  summarise(lai_p10 = quantile(lai, 0.1, na.rm = TRUE), 
            lai_p90 = quantile(lai, 0.9, na.rm = TRUE), 
            .groups = "drop")

env_df$age_class <- stringr::str_replace(env_df$age_class, "_g0", "")
env_df$age_class <- stringr::str_replace(env_df$age_class, "_g1", "")

pred_df$age_class <- stringr::str_replace(pred_df$age_class, "_g0", "")
pred_df$age_class <- stringr::str_replace(pred_df$age_class, "_g1", "")

pred_df <- pred_df %>% 
  mutate(lai_grp = factor(lai_grp, 
                       levels = c("g1", "g0"), 
                       labels = c(expression(Elevation>=z[cp]),
                                  expression(Elevation<z[cp])))) 

pred_df <- pred_df %>%
  mutate(age_class = factor(age_class,
                            levels = c("D", "R1", "R2", "R3", "R4", "R5", "R6",
                                       "R7", "R8", "R9", "M"),
                            labels = c(expression(0-5), 
                                       expression(6-10),
                                       expression(11-15),
                                       expression(16-20),
                                       expression(21-25),
                                       expression(26-30),
                                       expression(31-35),
                                       expression(36-40),
                                       expression(41-45),
                                       expression(46-50),
                                       expression('>'*50))))
                                
ggplot() +
  # geom_ribbon(data = env_df, aes(x = month_num, ymin = lai_p10, ymax = lai_p90),
  #   fill = "red", alpha = 0.4) +
  geom_line(data = pred_df, aes(x = month_num, y = lai_fit, colour = age_class), linewidth = 1) +
  facet_grid(~lai_grp, labeller = as_labeller(label_parsed)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7),
                     labels = c( "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Winter")) +
  labs(x = "Month", 
       y = "LAI",
       title = "Seasonal LAI") +
  theme_minimal() +
  theme(legend.position = "bottom")

lai_plot <- ggplot(pred_df, aes(x = month_num, y = lai_fit, group = age_class)) +
  geom_line() +
  facet_grid(
    lai_grp ~ age_class,
    scales = "free_x",
    space  = "free_x",
    labeller = as_labeller(label_parsed)
  ) +
  scale_x_continuous(
    breaks = c(1, 4, 7),
    labels = c("M", "A", "W")
  ) +
  labs(
    x = "Month",
    y = "LAI"
  ) +
  theme_bw() + 
  theme(strip.background = element_rect("NA"))
  
ggsave(lai_plot, 
       filename = file.path(fig_path, "fitted_lai_mar16.png"),
       dpi = 1200, 
       units = "in", 
       width = 7,
       height = 3)

# ------------------------------------------------------------------------------
# final estimates 
# ------------------------------------------------------------------------------

pred_data_final <- lai_models %>% 
  mutate(
    newdata = map(
      data, 
      ~ data.frame(month_num = seq(min(.x$month_num), max(.x$month_num), length.out = 7)),
    ), 
    preds = map2(model, newdata, predict)
  ) %>% 
  unnest(c(newdata, preds))

pred_df_final <- pred_data_final %>% 
  select(c("age_class", "lai_grp", "preds", "month_num"))

names(pred_df_final) <- c("age_class", "lai_grp", "lai", "month_num")

write.csv(pred_df_final, 
          file.path(result_path, "median_LAI_ALS_ELEV_recovery_feb19.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

assign_ageclass <- function(df, zoned_classes, class_type) {
  
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
