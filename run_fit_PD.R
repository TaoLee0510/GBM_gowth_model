base_output <- "/Users/4482173/Documents/Project/GBM_Model"

source("/Users/4482173/Library/CloudStorage/OneDrive-MoffittCancerCenter/GitHub/GBM_gowth_model/fit_PD.R")
#
# 1) Fit all simulations under a specific supply together:
fit_PD_from_results(file.path(base_output, "Results"),
                    exclude_daily = TRUE,   # Exclude hours where step%%24==0
                    bins = 30,
                    fit_gam = TRUE,
                    fd = 1.0,               # Pass your Fd; use NA if unknown
                    rt = 0.02,          # Pass your Rt; use NA if you don't want to use it
                    death_mode = "resource_only")              # c("resource_only", "all", "random_only", "timeout_only")

# 2) Or manually read the events data for one supply and fit:
ev <- read_events_dir(file.path(base_output, "Results", "0.1", "Sim_001", "events"))
fit_PD_from_events(ev, out_prefix = file.path(base_output, "Results", "0.1", "Sim_001", "PD_fit_single"),bins = 100)



base_output <- "/Users/4482173/Documents/Project/GBM_Model"
summarize_PD_all(
  results_root   = file.path(base_output, "Results"),
  exclude_daily  = TRUE,
  bins           = 100,
  fit_gam        = TRUE,
  fd             = NA,
  rt             = 0.02,
  death_mode     = "resource_only",  # or "all" / "random_only" / "timeout_only"
  g_points       = 200,
  overwrite_fit  = FALSE             # Set to TRUE to rerun the _ALL fit for each supply
)

# ---- Diversity time series analysis ----
# This script orchestrates the diversity analysis and per-simulation plots.
# ---- Diversity & Time-Series Analysis (single call) ----
source("/Users/4482173/Documents/GitHub/GBM_gowth_model/Diversity_TimeSeries_Analysis.R")
run_diversity_timeseries(
  base_output = base_output,
  top_n_karyotypes = 12,
  sample_cells = 5000,
  do_diversity = TRUE,
  do_per_sim  = TRUE
)


# ---- Decomposition examples: Path A and Path B ----
# ev_one <- read_events_dir(file.path(base_output, "Results", "0.25", "Sim_001", "events"))
# Path A: baseline calibration (group by karyotype if available, else Label)
# PD_decompose_pathA(ev_one,
#   out_prefix = file.path(base_output, "Results", "0.25", "Sim_001", "PD_decomp_pathA"),
#   bins = 30, group_field = NULL, suff_threshold = 1.0,
#   exclude_daily = TRUE, death_mode = "resource_only")
# Path B: GAM-based decomposition
# PD_decompose_pathB_gam(ev_one,
#   out_prefix = file.path(base_output, "Results", "0.25", "Sim_001", "PD_decomp_pathB"),
#   group_field = NULL, k = 10, exclude_daily = TRUE, death_mode = "resource_only")