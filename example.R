
# Example call (adjust paths)
# result <- run_simulation("config.yaml","merged_mean.csv",".")

setwd('/Users/4482173/Library/CloudStorage/OneDrive-MoffittCancerCenter/GitHub/GBM_gowth_model/')

source("GBM_abm.R")

run_ABM_simulation(
  cfg_file = "/Users/4482173/Library/CloudStorage/OneDrive-MoffittCancerCenter/GitHub/GBM_gowth_model/config.yaml",
  karyolib_file = "/Users/4482173/Library/CloudStorage/OneDrive-MoffittCancerCenter/GitHub/GBM_gowth_model/merged_mean.csv",
  base_output = "/Users/4482173/Documents/Project/GBM_Model", 
  workers = parallel::detectCores() - 2 
)