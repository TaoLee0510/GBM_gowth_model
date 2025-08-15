
# Example call (adjust paths)
# result <- run_simulation("config.yaml","merged_mean.csv",".")

setwd('/Users/taolee/Documents/GitHub/GBM_gowth_model/')

source("GBM_abm.R")

run_ABM_simulation(
  cfg_file = "/Users/taolee/Documents/GitHub/GBM_gowth_model/config.yaml",
  karyolib_file = "/Users/taolee/Documents/GitHub/GBM_gowth_model/merged_mean.csv",
  base_output = '/Volumes/Protable Disk/Project/GBM/ABM', 
  workers = parallel::detectCores() - 2 
)