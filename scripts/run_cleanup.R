# scripts/run_cleanup.R
# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway


source("R/cleanup_gmm.R")

# EDIT THESE TWO PATHS:
input_dir  <- "data/raw_fcs"        # local/SAFE folder with .fcs files (ignored by git)
output_dir <- "outputs/cleanup"     # local outputs (ignored by git)

#cleaning parameters are barcode-dependent
params <- default_cleanup_params()

# Example: tweak cutoffs if needed
# params$beads_cut <- 55

run_cleanup(input_dir, output_dir, params = params)
