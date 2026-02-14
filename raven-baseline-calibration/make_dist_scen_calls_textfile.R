# ------------------------------------------------------------------------------
# 
# make raven run calls text file
#
# Dec 3rd, 2025
# Brianne Boufford
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# libraries 

library(terra)
library(sf)
library(stringr)
library(tools)

# ------------------------------------------------------------------------------
# paths 

project_path <- file.path("C:", "Users", "blbouf", "Sync", "Paper2")
setwd(project_path)

# raven run  path
raven_run_path <- file.path(".", "raven-runs", "disturbance_recovery_scenarios")

calls_file <- file.path(raven_run_path, "Trapping_dist_scen_calls_dec4.txt")
calls_file_new <- file.path(raven_run_path, "Trapping_dist_scen_calls_dec4.txt")

# ---- READ FIRST LINE ----
lines <- readLines(calls_file)
template_line <- lines[1]

# ---- IDENTIFY .rvh FILES IN SAME DIRECTORY ----
dir_path <- dirname(calls_file)
rvh_files <- list.files(dir_path, pattern = "\\.rvh$", full.names = FALSE)
rvh_files <- rvh_files[rvh_files != "Trapping_all_forest.rvh"]

# Remove the template rvh file itself, if desired:
# rvh_files <- setdiff(rvh_files, "Trapping_all_forest.rvh")

# ---- GENERATE NEW CALLS ----
make_new_line <- function(f) {
  
  # base filename (no extension)
  base <- file_path_sans_ext(f)
  
  # replace the RVH argument
  new <- str_replace(template_line, "Trapping_all_forest\\.rvh", f)
  
  # replace the -h argument (if the rvh name appears there too)
  new <- str_replace(new, "Trapping_all_forest", base)
  
  # replace the Runs folder output path
  new <- str_replace(new, "Runs_dec4\\\\Trapping_all_forest", paste0("Runs_dec4\\", base))
  
  return(new)
}

list_lines <- lapply(rvh_files, make_new_line) %>%
  do.call(rbind,.)

# ---- WRITE BACK ----
# Option A: overwrite entire file with template + new lines
writeLines(c(template_line, list_lines), calls_file_new)
