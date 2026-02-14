# ------------------------------------------------------------------------------
# fit parametric curves to ALS-derived forest parameters 
#
# date: February 11, 2026 
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
  group_by(age, lai_grp) %>% 
  summarize(n = length(height),
            h_p10 = quantile(height, 0.1, na.rm = TRUE),
            h_p90 = quantile(height, 0.9, na.rm = TRUE),
            h_med = median(height, na.rm = TRUE),
            h_mean = mean(height, na.rm = TRUE),
            fc_p10 = quantile(fc, 0.1, na.rm = TRUE),
            fc_p90 = quantile(fc, 0.9, na.rm = TRUE),
            fc_med = median(fc, na.rm = TRUE),
            fc_mean = mean(fc, na.rm = TRUE),
            als_grp_key = first(als_grp_key)) %>% 
  subset(n > 100)

# get list of unique lai group - age IDs
age_veg_ids_keep <- unique(als_summary_df$als_grp_key)

# remove ones that are not kept with n < 100 filtering
als_veg_data <- als_veg_data[als_veg_data$als_grp_key %in% age_veg_ids_keep, ]

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

# ------------------------------------------------------------------------------
# get average, quantile, and N 
# ------------------------------------------------------------------------------

av_als_vals <- als_veg_data %>%
  group_by(age, lai_grp) %>%
  summarize(av_h = mean(height, na.rm = TRUE),
            n_h = length(height),
            h_p10 = quantile(height, 0.1, na.rm = TRUE), 
            h_p90 = quantile(height, 0.9, na.rm = TRUE),
            av_fc = mean(fc, na.rm = TRUE), 
            n_fc = length(fc),
            h_p10 = quantile(fc, 0.1, na.rm = TRUE), 
            h_p90 = quantile(fc, 0.9, na.rm = TRUE))

# ------------------------------------------------------------------------------
# Curve fitting 
# ------------------------------------------------------------------------------
# als_veg_data <- als_veg_data[als_veg_data$age != 0, ]
# keeping age 0 for now
g0_all <- subset(als_veg_data, lai_grp == "g0")
g0_all <- g0_all[!is.na(g0_all$age), ]
g0_all <- g0_all[!is.na(g0_all$height), ]
g0_all <- g0_all[!is.na(g0_all$fc), ]

g1_all <- subset(als_veg_data, lai_grp == "g1") 
g1_all <- g1_all[!is.na(g1_all$age), ]
g1_all <- g1_all[!is.na(g1_all$height), ]
g1_all <- g1_all[!is.na(g1_all$fc), ]

# ------------------------------------------------------------------------------
# g0 HEIGHT 
# ------------------------------------------------------------------------------

g0_grouped_all <- groupedData(height ~ age | lai_grp, data = g0_all)
g0_grouped <- groupedData(av_h ~ age | lai_grp, data = av_als_vals[av_als_vals$lai_grp == "g0", ])

g0_grouped_all <- g0_grouped_all[!g0_grouped_all$age == 0, ]

modseltable <- pn.mod.compare(g0_grouped_all$age, g0_grouped_all$height, g0_grouped_all$lai_grp, existing = FALSE, pn.options = "g0_h")

# model parameters 
modpar(g0_grouped_all$age, g0_grouped_all$height, pn.options = "g0_h")

# 23, 31, 32, 36, 25

richardsR31.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   RAsym = RAsym,
                                                   modno = 31, pn.options = "g0_h"), data = g0_grouped_all)

richardsR36.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   Rk = Rk, 
                                                   modno = 36, pn.options = "g0_h"), data = g0_grouped_all)

richardsR25.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   RM = RM, 
                                                   modno = 25, pn.options = "g0_h"), data = g0_grouped_all)

richardsR32.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   modno = 32, pn.options = "g0_h"), data = g0_grouped_all)

richardsR36.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   Rk = Rk, 
                                                   modno = 36, pn.options = "g0_h"), data = g0_grouped_all)

richardsR25.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   RM = RM, 
                                                   modno = 25, pn.options = "g0_h"), data = g0_grouped_all)

richardsR32.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                   Infl = Infl,
                                                   modno = 32, pn.options = "g0_h"), data = g0_grouped_all)

plot(g0_grouped$age , g0_grouped$av_h, xlab = "age", ylab = "average Heght")
plot(g0_grouped_all$age , g0_grouped_all$height, xlab = "age", ylab = "average Height")

par36 <- coef(richardsR36.nls)
curve(posnegRichards.eqn(x, Asym = par36[1], K = par36[2], Infl = par36[3],
                         Rk = par36[4],
                         modno = 36, pn.options = "g0_h"), add= TRUE)

par25 <- coef(richardsR25.nls)
curve(posnegRichards.eqn(x, Asym = par25[1], K = par25[2], Infl = par25[3],
                         RM = par25[4], 
                         modno = 25, pn.options = "g0_h"), add= TRUE)

par31 <- coef(richardsR31.nls)
curve(posnegRichards.eqn(x, Asym = par31[1], K = par31[2], Infl = par31[3],
                         RAsym = par31[4],
                         modno = 31, pn.options = "g0_h"), add= TRUE)

par32 <- coef(richardsR32.nls)
curve(posnegRichards.eqn(x, Asym = par32[1], K = par32[2], Infl = par32[3],
                         modno = 32, pn.options = "g0_h"), add= TRUE)

# ------------------------------------------------------------------------------
# g1 HEIGHT 
# ------------------------------------------------------------------------------
# 31, 36, 32 

g1_grouped_all <- groupedData(height ~ age | lai_grp, data = g1_all)
g1_grouped <- groupedData(av_h ~ age | lai_grp, data = av_als_vals[av_als_vals$lai_grp == "g1", ])

g1_grouped_all <- g1_grouped_all[!g1_grouped_all$age == 0, ]

modseltable <- pn.mod.compare(g1_grouped_all$age, g1_grouped_all$height, g1_grouped_all$lai_grp, existing = FALSE, pn.options = "g1_h")

# model parameters 
modpar(g1_grouped_all$age, g1_grouped_all$height, pn.options = "g1_h")

richardsR31.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 RAsym = RAsym,
                                                 modno = 31, pn.options = "g1_h"), data = g1_grouped_all)

richardsR36.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 Rk = Rk, 
                                                 modno = 36, pn.options = "g1_h"), data = g1_grouped_all)

richardsR32.nls <- nls(height ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 modno = 32, pn.options = "g1_h"), data = g1_grouped_all)

plot(g1_grouped$age , g1_grouped$av_h, xlab = "age", ylab = "average Heght g1")
plot(g1_grouped_all$age , g1_grouped_all$height, xlab = "age", ylab = "average Height")

par36 <- coef(richardsR36.nls)
curve(posnegRichards.eqn(x, Asym = par36[1], K = par36[2], Infl = par36[3],
                         Rk = par36[4],
                         modno = 36, pn.options = "g1_h"), add= TRUE)

par31 <- coef(richardsR31.nls)
curve(posnegRichards.eqn(x, Asym = par31[1], K = par31[2], Infl = par31[3],
                         RAsym = par31[4],
                         modno = 31, pn.options = "g1_h"), add= TRUE) 

par32 <- coef(richardsR32.nls)
curve(posnegRichards.eqn(x, Asym = par32[1], K = par32[2], Infl = par32[3],
                         modno = 32, pn.options = "g0_h"), add= TRUE) # this one

# ------------------------------------------------------------------------------
# g0 FCOVER  
# ------------------------------------------------------------------------------

g0_grouped_all <- groupedData(fc ~ age | lai_grp, data = g0_all)
g0_grouped <- groupedData(av_fc ~ age | lai_grp, data = av_als_vals[av_als_vals$lai_grp == "g0", ])

g0_grouped_all <- g0_grouped_all[!g0_grouped_all$age == 0, ]

modseltable <- pn.mod.compare(g0_grouped_all$age, g0_grouped_all$fc, g0_grouped_all$lai_grp, existing = FALSE, pn.options = "g0_fc")

plot(g0_grouped$age , g0_grouped$av_fc, xlab = "age", ylab = "average fc g0")
plot(g0_grouped_all$age , g0_grouped_all$fc, xlab = "age", ylab = "average Height")

richardsR11.nls <- nls(fc ~ SSposnegRichards(age, 
                                            M = M,
                                            Asym = Asym,
                                            Infl = Infl, 
                                            K = K, 
                                            RAsym = RAsym, 
                                            modno = 11, pn.options = "g0_fc"), data = g0_grouped_all)

richardsR12.nls <- nls(fc ~ SSposnegRichards(age, 
                                            M = M, 
                                            Asym = Asym,
                                            Infl = Infl, 
                                            K = K, 
                                            modno = 12, pn.options = "g0_fc"), data = g0_grouped_all)

par11 <- coef(richardsR11.nls)
curve(posnegRichards.eqn(x, Asym = par11[1], K = par11[2], Infl = par11[3],
                         M = par11[4], RAsym = par11[5],
                         modno = 11, pn.options = "g0_fc"), add= TRUE)

par12 <- coef(richardsR12.nls)
curve(posnegRichards.eqn(x, Asym = par12[1], K = par12[2], Infl = par12[3],
                         M = par12[4],
                         modno = 12, pn.options = "g0_fc"), add= TRUE) 
#  11, 12

# ------------------------------------------------------------------------------
# g1 FCOVER  
# ------------------------------------------------------------------------------

g1_grouped_all <- groupedData(fc ~ age | lai_grp, data = g1_all)
g1_grouped <- groupedData(av_fc ~ age | lai_grp, data = av_als_vals[av_als_vals$lai_grp == "g1", ])

g1_grouped_all <- g1_grouped_all[!g1_grouped_all$age == 0, ]

modseltable <- pn.mod.compare(g1_grouped_all$age, g1_grouped_all$fc, g1_grouped_all$lai_grp, existing = FALSE, pn.options = "g1_fc")

# 36 31 32!

richardsR31.nls <- nls(fc ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 RAsym = RAsym,
                                                 modno = 31, pn.options = "g1_fc"), data = g1_grouped_all)

richardsR36.nls <- nls(fc ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 Rk = Rk, 
                                                 modno = 36, pn.options = "g1_fc"), data = g1_grouped_all)

richardsR32.nls <- nls(fc ~ SSposnegRichards(age, Asym = Asym, K = K,
                                                 Infl = Infl,
                                                 modno = 32, pn.options = "g1_fc"), data = g1_grouped_all)

plot(g1_grouped$age , g1_grouped$av_fc, xlab = "age", ylab = "average fc g1")
plot(g1_grouped_all$age , g1_grouped_all$fc, xlab = "age", ylab = "average fc")

par36 <- coef(richardsR36.nls)
curve(posnegRichards.eqn(x, Asym = par36[1], K = par36[2], Infl = par36[3],
                         Rk = par36[4],
                         modno = 36, pn.options = "g1_fc"), add= TRUE)

par31 <- coef(richardsR31.nls)
curve(posnegRichards.eqn(x, Asym = par31[1], K = par31[2], Infl = par31[3],
                         RAsym = par31[4],
                         modno = 31, pn.options = "g1_fc"), add= TRUE) 

par32 <- coef(richardsR32.nls)
curve(posnegRichards.eqn(x, Asym = par32[1], K = par32[2], Infl = par32[3],
                         modno = 32, pn.options = "g0_fc"), add= TRUE) # this one
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
