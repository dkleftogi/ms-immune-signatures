source("R/harmonise_cycombine.R")

input_dir    <- "outputs/cleanup/clean_fcs"               # ignored by git
metadata_csv <- "metadata/metadata.csv"
markers_file <- "metadata/panel_design.xlsx"
output_dir   <- "outputs/cycombine"

params <- default_cycombine_params()
# params$condition_col <- "Cohort"   # if you prefer Cohort
# params$norm_method <- "scale"

run_cycombine_harmonisation(input_dir, metadata_csv, markers_file, output_dir, params)
