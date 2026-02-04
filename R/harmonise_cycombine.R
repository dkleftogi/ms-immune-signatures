# R/harmonise_cycombine.R

# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway

suppressPackageStartupMessages({
  library(cyCombine)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(readxl)
  library(ggridges)
})

default_cycombine_params <- function() {
  list(
    # prepare_data
    fcs_pattern = "\\.fcs$",
    filename_col = "Filename",
    batch_col = "batch",
    condition_col = "condition",   # set NULL if not present
    derand = TRUE,
    down_sample = FALSE,
    asinh_cofactor = 5,
    
    # batch_correct
    norm_method = "scale",
    xdim = 8,
    ydim = 8,
    
    # QC plotting
    n_cells_plot = 20000,
    seed = 1
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Simple marker list reader: one marker per line
read_markers <- function(markers_file) {
  panel <- readxl::read_excel(markers_file)
  markers <- panel %>%
    filter(Type != "none") %>%
    pull(Marker) %>%
    str_remove_all("[ _-]")
  markers
}

# Ridge plots for before/after (optional but useful)
make_ridge_plots <- function(df, markers, group_col, out_pdf, n_cells = 20000, seed = 1) {
  set.seed(seed)
  df <- as.data.frame(df)
  
  # downsample for plotting
  if (nrow(df) > n_cells) {
    df <- df[sample(seq_len(nrow(df)), n_cells), ]
  }
  
  ensure_dir(dirname(out_pdf))
  pdf(out_pdf, width = 12, height = 12, onefile = TRUE)
  
  for (mk in markers) {
    if (!mk %in% colnames(df)) next
    tmp <- df %>%
      dplyr::select(all_of(group_col), all_of(mk)) %>%
      rename(group = all_of(group_col), expr = all_of(mk)) %>%
      mutate(group = as.factor(group))
    
    p <- ggplot(tmp, aes(x = expr, y = group, fill = group)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 1) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(title = mk, x = "expression", y = 'batch')
    
    print(p)
  }
  dev.off()
}

# Main function
run_cycombine_harmonisation <- function(input_dir,
                                        metadata_csv,
                                        markers_file,
                                        output_dir,
                                        params = default_cycombine_params(),
                                        write_outputs = TRUE) {
  
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  
  ensure_dir(output_dir)
  qc_dir <- file.path(output_dir, "qc")
  ensure_dir(qc_dir)
  
  # Read metadata + clean whitespace
  meta <- readr::read_csv(metadata_csv, show_col_types = FALSE) %>%
    mutate(across(where(is.character), ~ str_trim(.)))
  
  # Validate required columns
  req <- c(params$filename_col, params$batch_col)
  miss <- setdiff(req, colnames(meta))
  if (length(miss) > 0) stop("Metadata missing required columns: ", paste(miss, collapse = ", "))
  
  if (!is.null(params$condition_col) && !(params$condition_col %in% colnames(meta))) {
    stop("condition_col set to '", params$condition_col, "' but not found in metadata.")
  }
  
  markers <- read_markers(markers_file)
  if (length(markers) == 0) stop("No markers found in markers_file.")
  
  set.seed(params$seed)
  
  # Prepare data (cyCombine format)
  dfci <- cyCombine::prepare_data(
    data_dir = input_dir,
    metadata = meta,
    pattern = params$fcs_pattern,
    filename_col = params$filename_col,
    batch_ids = params$batch_col,
    condition = if (is.null(params$condition_col)) NULL else params$condition_col,
    markers = markers,
    derand = params$derand,
    down_sample = params$down_sample
  )
  
  # Optional: before-correction ridge plots by batch
  make_ridge_plots(
    df = dfci,
    markers = markers,
    group_col = "batch",
    out_pdf = file.path(qc_dir, "ridges_uncorrected_by_batch.pdf"),
    n_cells = params$n_cells_plot,
    seed = params$seed
  )
  
  # Batch correction
  corrected <- dfci %>%
    cyCombine::batch_correct(
      covar = if (is.null(params$condition_col)) NULL else "condition",
      xdim = params$xdim,
      ydim = params$ydim,
      norm_method = params$norm_method,
      markers = markers
    )
  
  # Optional: after-correction ridge plots by batch
  make_ridge_plots(
    df = corrected,
    markers = markers,
    group_col = "batch",
    out_pdf = file.path(qc_dir, "ridges_corrected_by_batch.pdf"),
    n_cells = params$n_cells_plot,
    seed = params$seed
  )
  
  # Save outputs (code-only repo; outputs directory ignored)
  if (isTRUE(write_outputs)) {
    saveRDS(corrected, file.path(output_dir, "cycombine_corrected.rds"))
    #write.csv(meta, file.path(output_dir, "metadata_used.csv"), row.names = FALSE)
    
    # record parameters + session
    saveRDS(params, file.path(output_dir, "params_used.rds"))
    sink(file.path(output_dir, "sessionInfo.txt"))
    sessionInfo()
    sink()
  }
  
  invisible(list(raw = dfci, corrected = corrected, metadata = meta, markers = markers))
}
