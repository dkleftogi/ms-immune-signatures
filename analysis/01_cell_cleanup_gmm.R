#!/usr/bin/env Rscript
# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway


source("R/cleanup_gmm.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript analysis/01_cell_cleanup_gmm.R <input_fcs_dir> <output_dir>\n")
  quit(status = 1)
}

input_dir  <- args[[1]]
output_dir <- args[[2]]

run_cleanup(input_dir, output_dir, params = default_cleanup_params())
