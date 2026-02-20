# R/io_zenodo.R
# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway


download_if_missing <- function(url, dest_path) {
  if (file.exists(dest_path)) return(invisible(dest_path))
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  
  message("Downloading data from Zenodo...")
  utils::download.file(url, destfile = dest_path, mode = "wb", quiet = FALSE)
  message("Saved to: ", dest_path)
  
  invisible(dest_path)
}
