# R/cleanup_gmm.R
# Core cleanup pipeline for CyTOF FCS (MaxPar-style QC + GMM CD45/CD66b lymph gate)

# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(mclust)
  library(fpc)
  library(concaveman)
  library(cowplot)
})

# usually the parameters change based on the specific barcode of the data
#it is highly advisable to learn barcode-specific parameters and modify the code below
default_cleanup_params <- function() {
  list(
    beads_cut = 60,
    residual_min = 5,
    residual_max = 120,
    center_min = 220,
    center_max = 700,
    offset_min = 0.05,
    offset_max = 20,
    width_min = 80,
    width_max = 500,
    eventlen_min = 19,
    eventlen_max = 55,
    dna1_min = 250,
    dna1_max = 1350,
    dna2_min = 250,
    dna2_max = 2500,
    asinh_cofactor = 5,
    gmm_components = 2,
    silhouette_n = 10000,
    seed = 1,
    # channel names (edit here if needed)
    ch_time = "Time",
    ch_beads = "Ce140Di",
    ch_residual = "Residual",
    ch_center = "Center",
    ch_offset = "Offset",
    ch_width = "Width",
    ch_eventlen = "Event_length",
    ch_dna1 = "Ir191Di",
    ch_dna2 = "Ir193Di",
    ch_cd66b = "Sm152Di",
    ch_cd45 = "Y89Di"
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

stop_if_missing_cols <- function(exprs_mat, required, file_tag) {
  missing <- setdiff(required, colnames(exprs_mat))
  if (length(missing) > 0) {
    stop(sprintf("Missing required channels in %s: %s", file_tag, paste(missing, collapse = ", ")))
  }
}

save_plot <- function(p, out_path, w = 8, h = 8, dpi = 300) {
  ggsave(out_path, p, width = w, height = h, dpi = dpi, bg = "white")
}

#' Run cleanup on all FCS files in a directory
#' @param input_dir folder with .fcs files
#' @param output_dir output folder (plots/, clean_fcs/, qc/)
#' @param params list of parameters; start from default_cleanup_params()
#' @return list with qc_summary and cluster_quality data.frames
run_cleanup <- function(input_dir,
                        output_dir,
                        params = default_cleanup_params(),
                        pattern = "\\.fcs$",
                        write_clean_fcs = TRUE) {
  
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  
  plots_dir <- file.path(output_dir, "plots")
  clean_dir <- file.path(output_dir, "clean_fcs")
  qc_dir    <- file.path(output_dir, "qc")
  
  ensure_dir(output_dir)
  ensure_dir(plots_dir)
  ensure_dir(clean_dir)
  ensure_dir(qc_dir)
  
  set.seed(params$seed)
  
  fcs_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(fcs_files) == 0) stop("No .fcs files found in: ", input_dir)
  
  qc_summary <- data.frame(
    file = character(),
    n_raw = integer(),
    n_after_beads = integer(),
    n_after_residual = integer(),
    n_after_center = integer(),
    n_after_offset = integer(),
    n_after_width = integer(),
    n_after_eventlen = integer(),
    n_after_dna1 = integer(),
    n_after_dna2 = integer(),
    n_lymph_final = integer(),
    stringsAsFactors = FALSE
  )
  
  cluster_quality <- data.frame(
    file = character(),
    avg_silwidth = numeric(),
    dunn2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  required <- c(
    params$ch_time, params$ch_beads, params$ch_residual, params$ch_center,
    params$ch_offset, params$ch_width, params$ch_eventlen,
    params$ch_dna1, params$ch_dna2, params$ch_cd66b, params$ch_cd45
  )
  
  for (f in fcs_files) {
    file_tag <- basename(f)
    message("Processing: ", file_tag)
    fcs_data <- read.FCS(f, transformation = FALSE)
    
    # Gate 0: beads (Ce140 vs Time)
    tmp1 <- fcs_data@exprs
    n0 <- nrow(tmp1)
    dt1 <- data.frame(Ce140=tmp1[,'Ce140Di'],
                      Time=tmp1[,"Time"])
    n_raw <- nrow(tmp1)
    o0 <- ggplot(dt1,aes(x=Time, y = log2(Ce140))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('Time')+
      ylab('Ce140 - log2')+
      geom_hline(yintercept = log2(params$beads_cut),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    indices_0 <- which(dt1[, "Ce140"] <= params$beads_cut )
    if (length(indices_0) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    fcs_data <- fcs_data[indices_0, ]
    tmp1 <- fcs_data@exprs
    n1 <- nrow(tmp1)
    dt1 <- data.frame(Residual=tmp1[,'Residual'],
                      Time=tmp1[,"Time"])
    
    # Gate 1: residual
    Residual_cut_min <- params$residual_min
    Residual_cut_max <- params$residual_max
    #visualise Residuals vs. Time
    o1 <- ggplot(dt1,aes(x=Time, y = log2(Residual))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('Time')+
      ylab('Residual - log2')+
      geom_hline(yintercept = log2(Residual_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(Residual_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    indices_1 <- which(dt1[, "Residual"] >= Residual_cut_min &
                         dt1[, "Residual"] <= Residual_cut_max)
    if (length(indices_1) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_1 <- fcs_data[indices_1, ]
    
    # Gate 2: center
    tmp1 <- filter_1@exprs
    n2 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      Center=tmp1[,'Center'])
    Center_cut_min <- params$center_min
    Center_cut_max <- params$center_max
    #visualise Center vs. Time
    o2 <- ggplot(dt1,aes(x=Time, y = log2(Center))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('Time')+
      ylab('Center - log2')+
      geom_hline(yintercept = log2(Center_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(Center_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    indices_2 <- which(dt1[, "Center"] >= Center_cut_min &
                         dt1[, "Center"] <= Center_cut_max)
    if (length(indices_2) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_2 <- filter_1[indices_2, ]
    
    # Gate 3: offset
    tmp1 <- filter_2@exprs
    n3 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      Offset=tmp1[,'Offset'])
    Offset_cut_min <- params$offset_min
    Offset_cut_max <- params$offset_max
    #visualise offset vs. Time
    o3 <- ggplot(dt1,aes(x=Time, y = log2(Offset))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('Time')+
      ylab('Offset - log2')+
      geom_hline(yintercept = log2(Offset_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(Offset_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    indices_3 <- which(dt1[, "Offset"] >= Offset_cut_min &
                         dt1[, "Offset"] <= Offset_cut_max)
    if (length(indices_3) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_3 <- filter_2[indices_3, ]
    
    # Gate 4: width
    tmp1 <- filter_3@exprs
    n4 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      Width=tmp1[,'Width'])
    
    Width_cut_min <- params$width_min
    Width_cut_max <- params$width_max
    #visualise width vs. Time
    o4 <- ggplot(dt1,aes(x=Time, y = log2(Width))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('Time')+
      ylab('Width - log2')+
      ylim(0,log2(5000))+
      geom_hline(yintercept = log2(Width_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(Width_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    indices_4 <- which(dt1[, "Width"] >= Width_cut_min &
                         dt1[, "Width"] <= Width_cut_max)
    if (length(indices_4) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_4 <- filter_3[indices_4, ]
    
    # Gate 5: event length
    tmp1 <- filter_4@exprs
    n5 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      Event_length=tmp1[,'Event_length'])
    Event_length_cut_min <- params$eventlen_min
    Event_length_cut_max <- params$eventlen_max
    #visualise Event_length vs. Time
    o5 <- ggplot(dt1,aes(x=Time, y = (Event_length))) +
      geom_point(size=0)+
      geom_hex(bins = 80) +
      xlab('Time')+
      ylab('Event_length')+
      geom_hline(yintercept = (Event_length_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = (Event_length_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    o5a <- ggplot(dt1, aes(x = Event_length)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 2, fill = "white",binwidth = 10) +
      geom_density(lwd = 0.5, colour = 6,
                   fill = 4, alpha = 0.25)+
      xlab('Event_length')+
      xlim(0,200)+
      geom_vline(xintercept =Event_length_cut_min,color='red',linetype='dashed',linewidth=0.2)+
      geom_vline(xintercept =Event_length_cut_max,color='red',linetype='dashed',linewidth=0.2)+
      #ggtitle(files_names[idx])+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),aspect.ratio = 1) 
    indices_5 <- which(dt1[, "Event_length"] >= Event_length_cut_min &
                         dt1[, "Event_length"] <= Event_length_cut_max)
    if (length(indices_5) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_5 <- filter_4[indices_5, ]
    
    # Gate 6: DNA1
    tmp1 <- filter_5@exprs
    n6 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      DNA1=tmp1[,'Ir191Di'])
    
    DNA1_cut_min <- params$dna1_min
    DNA1_cut_max <- params$dna1_max
    #visualise DNA1 vs. Time
    o6 <- ggplot(dt1,aes(x=Time, y = log2(DNA1))) +
      geom_point(size=0)+
      geom_hex(bins = 80) +
      xlab('Time')+
      ylab('DNA1 - log2')+
      ylim(0,log2(10000))+
      geom_hline(yintercept = log2(DNA1_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(DNA1_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    o6a <- ggplot(dt1, aes(x = DNA1)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 2, fill = "white",binwidth = 50) +
      geom_density(lwd = 0.5, colour = 6,
                   fill = 4, alpha = 0.25)+
      xlab('DNA1')+
      xlim(0,10000)+
      geom_vline(xintercept =DNA1_cut_min,color='red',linetype='dashed',linewidth=0.2)+
      geom_vline(xintercept =DNA1_cut_max,color='red',linetype='dashed',linewidth=0.2)+
      #ggtitle(files_names[idx])+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),aspect.ratio = 1) 
    indices_6 <- which(dt1[, "DNA1"] >= DNA1_cut_min &
                         dt1[, "DNA1"] <= DNA1_cut_max)
    if (length(indices_6) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    filter_6 <- filter_5[indices_6, ]
    
    # Gate 7: DNA2
    tmp1 <- filter_6@exprs
    n7 <- nrow(tmp1)
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      DNA1=tmp1[,'Ir191Di'],
                      DNA2=tmp1[,'Ir193Di'])
    
    DNA2_cut_min <- params$dna2_min
    DNA2_cut_max <- params$dna2_max
    #visualise DNA2 vs. Time
    o7 <- ggplot(dt1,aes(x=Time, y = (DNA2))) +
      geom_point(size=0)+
      geom_hex(bins = 80) +
      xlab('Time')+
      ylab('DNA1 - log2')+
      #ylim(0,log2(10000))+
      geom_hline(yintercept = (DNA2_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = (DNA2_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    o7a <- ggplot(dt1, aes(x = DNA2)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 2, fill = "white",binwidth = 50) +
      geom_density(lwd = 0.5, colour = 6,
                   fill = 4, alpha = 0.25)+
      xlab('DNA2')+
      xlim(0,10000)+
      geom_vline(xintercept =DNA2_cut_min,color='red',linetype='dashed',linewidth=0.2)+
      geom_vline(xintercept =DNA2_cut_max,color='red',linetype='dashed',linewidth=0.2)+
      #ggtitle(files_names[idx])+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),aspect.ratio = 1) 
    
    o7b <- ggplot(dt1,aes(x=log2(DNA1), y = log2(DNA2))) +
      geom_point(size=0)+
      geom_hex(bins = 200) +
      xlab('DNA1 - log2')+
      ylab('DNA2 - log2')+
      #ylim(0,log2(10000))+
      geom_hline(yintercept = log2(DNA2_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_hline(yintercept = log2(DNA2_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      geom_vline(xintercept = log2(DNA1_cut_min),color='red',linetype='dashed',linewidth=0.5)+
      geom_vline(xintercept = log2(DNA1_cut_max),color='red',linetype='dashed',linewidth=0.5)+
      coord_fixed(ratio = 1) +
      #ggtitle(files_names[idx])+
      scale_fill_viridis(option = "A")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    indices_7 <- which(dt1[, "DNA2"] >= DNA2_cut_min &
                         dt1[, "DNA2"] <= DNA2_cut_max)
    if (length(indices_7) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    #singlets
    filter_7 <- filter_6[indices_7, ]
    
    # Gate 8: GMM on asinh(CD66b/CD45) -> pick CD66blo CD45hi cluster
    tmp1 <- filter_7@exprs
    dt1 <- data.frame(Time=tmp1[,"Time"],
                      DNA1=tmp1[,'Ir191Di'],
                      DNA2=tmp1[,'Ir193Di'],
                      CD66b=tmp1[,'Sm152Di'],
                      CD45=tmp1[,'Y89Di'])
    
    ## Extract CD45 and CD66b
    cft <- params$asinh_cofactor
    #here the trick is to transform the data with the standard transformation
    df <- dt1 %>%
      mutate(CD45  = asinh(CD45 / cft),
             CD66b = asinh(CD66b / cft))
    
    # fit GMM to CD45 / CD66b
    gmm <- Mclust(df[,c("CD66b","CD45")], G=params$gmm_components)  # 2 components: lymphocytes vs granulocytes
    df$cluster <- gmm$classification
    
    means <- gmm$parameters$mean
    covs <- gmm$parameters$variance$sigma
    probs <- gmm$z
    ##############################
    cluster_means <- df %>%
      group_by(cluster) %>%
      summarise(
        mean_cd45  = mean(CD45, na.rm=TRUE),
        mean_cd66b = mean(CD66b, na.rm=TRUE)
      )
    # lymphocytes = low CD66b, high CD45
    lymph_cluster <- cluster_means %>%
      arrange(mean_cd66b, desc(mean_cd45)) %>%
      slice(1) %>%
      pull(cluster)
    
    df$selected <- df$cluster == lymph_cluster
    
    #after this command the lympocytes cluster is the TRUE and the granulocytes it the FALSE
    
    # subset the selected cells
    sub <- df %>% dplyr::filter(selected)
    sum_matrix <- as.matrix(sub[,c("CD66b","CD45")])
    # concave hull polygon around selected cluster
    hull_poly_lymphocytes <- concaveman(sum_matrix)
    
    a <- which(df$selected==FALSE)
    sum_matrix <- as.matrix(df[a,c("CD66b","CD45")])
    # concave hull polygon around selected cluster
    hull_poly_granulocytes <- concaveman(sum_matrix)
    
    #compute shillhouete
    ###############################################
    idx1 <- sample(nrow(df))
    idx1 <- idx1[1:10000]
    a <- df[idx1,c(1,2)]
    b <- gmm$classification[idx1]
    bd <- cluster.stats(dist(a), b)
    bd <- data.frame(Sample=file_tag,
                     avgSilwidths=bd$avg.silwidth,
                     Dune=bd$dunn2)
    cluster_quality <- rbind(cluster_quality,bd)
    
    ###############################################
    o8 <- ggplot(df,aes(x=(CD66b), y =(CD45))) +
      geom_point(size=0)+
      geom_hex(bins = 100) +
      xlab('CD66b - asinh(5) ')+
      ylab('CD45 - asinh(5)')+
      #xlim(0,log2(2000))+
      #ylim(0,log2(2000))+
      geom_polygon(data = hull_poly_granulocytes, aes(x = V1, y = V2),
                   fill = NA, color = "blue", size = 1)+
      geom_polygon(data = hull_poly_lymphocytes, aes(x = V1, y = V2),
                   fill = NA, color = "red", size = 1)+
      scale_color_manual(values = c("red", "blue")) +
      coord_fixed(ratio = 1) +
      scale_fill_viridis(option = "A")+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),
            aspect.ratio = 1)
    ################################################################
    a1 <- which(df$selected==TRUE)
    a2 <- which(df$selected==FALSE)
    
    lymphocytes <- filter_7[a1, ]
    granulocytes <- filter_7[a2, ]
    ###############################################################
    
    title <- ggdraw() + 
      draw_label(file_tag, fontface = 'bold', x = 0.5, hjust = 0.5, size = 16)
    
    combined_plot <- cowplot::plot_grid(o1,o2,o3,
                                        o4,o5,o5a,
                                        ncol = 3, 
                                        labels = c("A", "B",'C',
                                                   'D','E','F'),
                                        rel_heights = c(1,1,1,
                                                        1,1,1))
    final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
    
    #myfile<-paste('cleanup_plots/',files_names[idx],'_cleaning.png',sep ='')
    #ggsave(myfile, final_plot, width = 12, height = 12, dpi = 300,bg = "white")
    save_plot(final_plot, file.path(plots_dir, paste0(tools::file_path_sans_ext(file_tag), "_cleaning.png")))
    
    #####################################################################
    #finally test the expression of gated lymphocytes and granulocytes using CD3 and CD19
    data_lymph <- lymphocytes@exprs
    dt1 <- data.frame(CD3=data_lymph[,"Cd116Di"],
                      CD19=data_lymph[,'Cd114Di'],
                      CD66b=data_lymph[,'Sm152Di'],
                      CD45=data_lymph[,'Y89Di'],
                      CD16=data_lymph[,'Bi209Di'],
                      Event_length=data_lymph[,'Event_length'])
    
    o8b <- ggplot(dt1, aes(x = CD45)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 2, fill = "white",binwidth = 50) +
      geom_density(lwd = 0.5, colour = 6,
                   fill = 4, alpha = 0.25)+
      geom_vline(xintercept =quantile(dt1$CD45,probs =0.005),color='red',linetype='dashed',linewidth=0.2)+
      xlab('CD45 - Lymphocytes')+
      xlim(0,2000)+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),aspect.ratio = 1) 
    
    indices_lymp <- which(dt1[, "CD45"] > quantile(dt1$CD45,probs =0.005) )
    if (length(indices_lymp) == 0) {
      message("No singlets found in ", file, " Skipping.")
      next
    }
    data_lymph0 <- lymphocytes[indices_lymp, ]
    data_lymph <- data_lymph0@exprs
    dt1 <- data.frame(CD3=data_lymph[,"Cd116Di"],
                      CD19=data_lymph[,'Cd114Di'],
                      CD66b=data_lymph[,'Sm152Di'],
                      CD45=data_lymph[,'Y89Di'],
                      CD16=data_lymph[,'Bi209Di'],
                      Event_length=data_lymph[,'Event_length'])
    n_lymph <- nrow(data_lymph)
    
    o9 <- ggplot(dt1,aes(x=(CD3), y =(CD19))) +
      geom_point(size=0)+
      geom_hex(bins = 150) +
      xlab('CD3')+
      ylab('CD19')+
      xlim(0,(1000))+
      ylim(0,(200))+
      coord_fixed(ratio = 1) +
      scale_fill_viridis(option = "C")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    data_granul <- granulocytes@exprs
    dt1 <- data.frame(CD3=data_granul[,"Cd116Di"],
                      CD19=data_granul[,'Cd114Di'],
                      CD66b=data_granul[,'Sm152Di'],
                      CD45=data_granul[,'Y89Di'],
                      CD16=data_granul[,'Bi209Di'],
                      Event_length=data_granul[,'Event_length'])
    
    o10 <- ggplot(dt1,aes(x=(CD3), y =(CD19))) +
      geom_point(size=0)+
      geom_hex(bins = 150) +
      xlab('CD3')+
      ylab('CD19')+
      xlim(0,(1000))+
      ylim(0,(200))+
      coord_fixed(ratio = 1) +
      scale_fill_viridis(option = "C")+
      #scale_fill_gradientn(colours = c("#330000", "#FFFFCC", "#660000", "#FF0000"), trans = "sqrt") +
      theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 10), axis.title = element_text(size = 10),
                        aspect.ratio = 1)
    
    o8a <- ggplot(dt1, aes(x = CD66b)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 2, fill = "white",binwidth = 50) +
      geom_density(lwd = 0.5, colour = 6,
                   fill = 4, alpha = 0.25)+
      xlab('CD66b')+
      xlim(0,2000)+
      theme_bw()+ 
      theme(legend.position="none", axis.text = element_text(size = 10), 
            axis.title = element_text(size = 10),aspect.ratio = 1) 
    
    combined_plot <- cowplot::plot_grid(o6a,o7a,o7b,
                                        o8,o8b,o9,
                                        o8a,o10,
                                        ncol = 3, 
                                        labels = c("A", "B",'C',
                                                   'D','E','Lymphocytes',
                                                   'F','Granulocytes'),
                                        rel_heights = c(1,1,1,
                                                        1,1,1,
                                                        1,1))
    final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
    #myfile<-paste('cleanup_plots/',files_names[idx],'_lymph_gate.png',sep ='')
    #ggsave(myfile, final_plot, width = 12, height = 12, dpi = 300,bg = "white")
    save_plot(final_plot, file.path(plots_dir, paste0(tools::file_path_sans_ext(file_tag), "_lymph_gate.png")))
    
    
    qc_summary <- rbind(qc_summary, data.frame(
      file = file_tag,
      n_raw = n_raw,
      n_after_beads = n0,
      n_after_residual = n1,
      n_after_center = n2,
      n_after_offset = n3,
      n_after_width = n4,
      n_after_eventlen = n5,
      n_after_dna1 = n6,
      n_after_dna2 = n7,
      n_lymph_final = n_lymph
    ))
    
    # Write cleaned FCS
    if (isTRUE(write_clean_fcs)) {
      out_fcs <- file.path(clean_dir, paste0(tools::file_path_sans_ext(file_tag), "_lymphocytes.fcs"))
      write.FCS(data_lymph0, filename = out_fcs)
    }
    
  }
  
  # Save QC tables + session info
  write.csv(qc_summary, file.path(qc_dir, "qc_counts_per_file.csv"), row.names = FALSE)
  write.csv(cluster_quality, file.path(qc_dir, "cluster_quality_per_file.csv"), row.names = FALSE)
  
  sink(file.path(qc_dir, "sessionInfo.txt"))
  sessionInfo()
  sink()

  invisible(list(qc_summary = qc_summary, cluster_quality = cluster_quality))
}
