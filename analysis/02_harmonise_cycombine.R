#!/usr/bin/env Rscript
# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway

source("R/harmonise_cycombine.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript analysis/02_harmonise_cycombine.R <input_dir> <metadata_csv> <markers_file> <output_dir>\n")
  quit(status = 1)
}

run_cycombine_harmonisation(
  input_dir = args[[1]],
  metadata_csv = args[[2]],
  markers_file = args[[3]],
  output_dir = args[[4]]
)
