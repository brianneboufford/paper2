# ------------------------------------------------------------------------------
# functions to make Raven Source Files 
# for Paper 2!
# ------------------------------------------------------------------------------

################################################################################
# : MAX_HT for : VegetationClasses
################################################################################
# get_max_height <- function(vegclasses){
#   
#   max_h <- c()
#   
#   for (i in 1:length(vegclasses)){
#     
#     if (grepl(pattern = "ALPINE", vegclasses[i]) |
#         grepl(pattern = "WET_LAND", vegclasses[i])){
#       
#       h <- 0.0 
#       
#     } else if (grepl(pattern = "SHRUB", vegclasses[i])){
#       
#       h <- 1.0
#       
#     } else if (grepl(pattern = "M_essf", vegclasses[i]) |
#                grepl(pattern = "M_idf", vegclasses[i]) |
#                grepl(pattern = "M_ms", vegclasses[i]) |
#                grepl(pattern = "R3_", vegclasses[i]) |
#                grepl(pattern = "R4_", vegclasses[i]) |
#                grepl(pattern = "R5_", vegclasses[i])){
#       
#       h <- 15.0 
#       
#     } else if (grepl(pattern = "D_", vegclasses[i])){
#       
#       h <- 0.5 
#       
#     } else if (grepl(pattern = "R1_", vegclasses[i]) | 
#                grepl(pattern = "R2_", vegclasses[i])){
#       
#       h <- 5.0
#       
#     }
#     max_h <- append(max_h, h)
#   }
#   
#   return(max_h)
#   
# }
################################################################################
# : SeasonalCanopyLAI ALL FORESTED 
################################################################################

# function to make seasonal canopy LAI that uses frst_cl not age_class to assign LAI
# for all forested run 
get_seasonal_canopy_lai_all_forest <- function(all_forested_lai){
  
  # SeaonalCanopyLAI
  name <- all_forested_lai$age_class
  jan <- as.numeric(all_forested_lai$med_winter_lai)
  
  seasonalcanopylai <- cbind(name, jan) %>%
    as.data.frame()
  seasonalcanopylai$jan <- as.numeric(seasonalcanopylai$jan)
  seasonalcanopylai$feb <- all_forested_lai$med_winter_lai 
  seasonalcanopylai$mar <- all_forested_lai$med_winter_lai
  seasonalcanopylai$apr <- all_forested_lai$med_winter_lai
  seasonalcanopylai$may <- all_forested_lai$med_may_lai
  seasonalcanopylai$jun <- all_forested_lai$med_june_lai
  seasonalcanopylai$jul <- all_forested_lai$med_july_lai
  seasonalcanopylai$aug <- all_forested_lai$med_aug_lai
  seasonalcanopylai$sep <- all_forested_lai$med_sept_lai
  seasonalcanopylai$oct <- all_forested_lai$med_oct_lai
  seasonalcanopylai$nov <- all_forested_lai$med_winter_lai
  seasonalcanopylai$dec <- all_forested_lai$med_winter_lai
  
  seasonalcanopylai[, -1] <- round(seasonalcanopylai[, -1] / seasonalcanopylai$jul , 2)
  
  return(seasonalcanopylai)
}

################################################################################
# : SeasonalCanopyLAI
################################################################################

# function to make seasonal canopy LAI for all distrubance runs
get_seasonal_canopy_lai <- function(all_forested_lai){
  
  # SeaonalCanopyLAI
  name <- all_forested_lai$age_class
  jan <- as.numeric(all_forested_lai$med_winter_lai)
  
  seasonalcanopylai <- cbind(name, jan) %>%
    as.data.frame()
  seasonalcanopylai$jan <- as.numeric(seasonalcanopylai$jan)
  seasonalcanopylai$feb <- all_forested_lai$med_winter_lai 
  seasonalcanopylai$mar <- all_forested_lai$med_winter_lai
  seasonalcanopylai$apr <- all_forested_lai$med_winter_lai
  seasonalcanopylai$may <- all_forested_lai$med_may_lai
  seasonalcanopylai$jun <- all_forested_lai$med_june_lai
  seasonalcanopylai$jul <- all_forested_lai$med_july_lai
  seasonalcanopylai$aug <- all_forested_lai$med_aug_lai
  seasonalcanopylai$sep <- all_forested_lai$med_sept_lai
  seasonalcanopylai$oct <- all_forested_lai$med_oct_lai
  seasonalcanopylai$nov <- all_forested_lai$med_winter_lai
  seasonalcanopylai$dec <- all_forested_lai$med_winter_lai
  
  seasonalcanopylai[, -1] <- round(seasonalcanopylai[, -1] / seasonalcanopylai$jul , 2)
  
  return(seasonalcanopylai)
}

################################################################################
# : Vegetation Parameter List
################################################################################

get_max_snow_c <- function(vegclasses){
  
  max_snow_c <- c()
  
  for (i in 1:length(vegclasses)){
    
    if (grepl(pattern = "ALPINE", vegclasses[i]) |
        grepl(pattern = "WET_LAND", vegclasses[i])){
      
      msc <- "_DEFAULT"
      
    } else if (grepl(pattern = "SHRUB", vegclasses[i])){
      
      msc <- 10
      
    } else if (grepl(pattern = "M_", vegclasses[i]) |
              grepl(pattern = "R5_", vegclasses[i]) |
              grepl(pattern = "R6_", vegclasses[i]) |
              grepl(pattern = "R7_", vegclasses[i]) |
              grepl(pattern = "R8_", vegclasses[i]) |
              grepl(pattern = "R9_", vegclasses[i])){
      
      msc <- "_DEFAULT"
      
    } else if (grepl(pattern = "D_", vegclasses[i]) |
               grepl(pattern = "R1_", vegclasses[i])){
      
      msc <- 5 
      
    } else if (grepl(pattern = "R2_", vegclasses[i]) | 
               grepl(pattern = "R3_", vegclasses[i]) | 
               grepl(pattern = "R4_", vegclasses[i])){
      
      msc <- 10
      
    }
    max_snow_c <- append(max_snow_c, msc)
  }
  
  return(max_snow_c)
  
}

################################################################################
# : LandUseClasses
################################################################################

# get_LU_classes <- function(vegclasses){
#   
#   fcover <- c()
#   imperm  <- c()
#   
#   for (i in 1:length(vegclasses)){
#     
#     if (grepl(pattern = "ALPINE", vegclasses[i])){
#       
#       im <- 0.0
#       fc <- 0.0
#       
#     } else if (grepl(pattern = "WET_LAND", vegclasses[i])){
#       
#       im <- 0.0
#       fc <- 0.5
#       
#     } else if (grepl(pattern = "SHRUB", vegclasses[i])){
#       
#       im <- 0.0
#       fc <- 0.6
#       
#     } else if (grepl(pattern = "R3_essf", vegclasses[i]) |
#                grepl(pattern = "R4_essf", vegclasses[i]) |
#                grepl(pattern = "R5_essf", vegclasses[i]) |
#                grepl(pattern = "M_essf", vegclasses[i])){
#       im <- 0.00
#       fc <- 0.75
#       
#     } else if (grepl(pattern = "R3_idf", vegclasses[i]) |
#                grepl(pattern = "R4_idf", vegclasses[i]) |
#                grepl(pattern = "R5_idf", vegclasses[i]) |
#                grepl(pattern = "M_idf", vegclasses[i]) |
#                grepl(pattern = "R3_ms", vegclasses[i]) |
#                grepl(pattern = "R4_ms", vegclasses[i]) |
#                grepl(pattern = "R5_ms", vegclasses[i]) |
#                grepl(pattern = "M_ms", vegclasses[i])){
#       
#       im <- 0.0 
#       fc <- 0.7
#       
#     } else if (grepl(pattern = "D_", vegclasses[i])){
#       
#       im <- 0.0
#       fc <- 0.6
#       
#     } else if (grepl(pattern = "R1_", vegclasses[i]) | 
#                grepl(pattern = "R2_", vegclasses[i])){
#       
#       im <- 0.00
#       fc <- 0.85
#       
#     } else if (grepl(pattern = "mature", vegclasses[i])){
#       im <- 0.0 
#       fc <- 0.7
#     }
#     imperm <- append(imperm, im)
#     fcover <- append(fcover, fc)
#   }
#   
#   r <- cbind(imperm, fcover)
#   return(r)
#   
# }

################################################################################
# : LandUseParameterList
################################################################################

get_LU_params <- function(vegclasses){
  
  hbv_melt_for_corr <- c()
  priestlytaylor  <- c()
  
  for (i in 1:length(vegclasses)){
    
    if (grepl(pattern = "ALPINE", vegclasses[i])){
      
      h <- "_DEFAULT"
      pt <- "_DEFAULT"
      
    } else if (grepl(pattern = "WET_LAND", vegclasses[i])){
      
      h <- "_DEFAULT"
      pt <- "_DEFAULT"
      
    } else if (grepl(pattern = "SHRUB", vegclasses[i])){
      
      h <- "_DEFAULT"
      pt <- "_DEFAULT"
      
    } else if (grepl(pattern = "R5_g1", vegclasses[i]) |
               grepl(pattern = "R6_g1", vegclasses[i]) |
               grepl(pattern = "R7_g1", vegclasses[i]) |
               grepl(pattern = "R8_g1", vegclasses[i]) |
               grepl(pattern = "R9_g1", vegclasses[i]) |
               grepl(pattern = "M_g1", vegclasses[i])){
      
      h <- 0.8
      pt <- 1.0
      
    } else if (grepl(pattern = "R5_g0", vegclasses[i]) |
               grepl(pattern = "R6_g0", vegclasses[i]) |
               grepl(pattern = "R7_g0", vegclasses[i]) |
               grepl(pattern = "R8_g0", vegclasses[i]) |
               grepl(pattern = "R9_g0", vegclasses[i]) |
               grepl(pattern = "M_g0", vegclasses[i])){
      
      h <- 0.8
      pt <- 1.1
      
      
    } else if (grepl(pattern = "D", vegclasses[i])){
      
      h <- "_DEFAULT"
      pt <- 1.25
      
    } else if (grepl(pattern = "R1_", vegclasses[i]) |
               grepl(pattern = "R2_", vegclasses[i]) | 
               grepl(pattern = "R3_", vegclasses[i]) |
               grepl(pattern = "R4_", vegclasses[i])){
      
      h <- 0.9
      pt <- 1.1
      
    }
    
    hbv_melt_for_corr <- append(hbv_melt_for_corr, h)
    priestlytaylor <- append(priestlytaylor, pt)
  }
  
  r <- cbind(hbv_melt_for_corr, priestlytaylor)
  return(r)
  
}
