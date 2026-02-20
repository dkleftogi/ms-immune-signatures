# R/smartms_config.R
# Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)
# This code is licensed under the MIT License. Copyright 2026, University of Bergen (UiB) and Neuro-SysMed, Norway


get_smartms_config <- function() {
  col.names <- c("CD16","CD4","CD14","CD19","CD3","CD235ab","CD11c","CD33",
                 "CD133","CD123","CD162","CD185","CD45RA","CD278","CD194","CD161",
                 "CD184","CD27","CD44","CD127","CD10","CD73","HLADR","CD146",
                 "CD117","CD8a","CD34","CD105","CD49d","CD20","CD25","CD66b",
                 "CD49f","CD45RO","CD90","CD45","CD195","CD38","CD196","CD135","CD56")
  
  col.names.info <- c('type','type','type','type','type','state','type','type',
                      'state','type','state','state','type','state','state','state',
                      'state','state','state','state','state','state','type','state',
                      'state','type','type','state','state','type','state','type',
                      'state','type','state','type','state','type','state','state','type')
  
  list(
    marker_info = data.frame(Channel = NA,
                             marker_name = col.names,
                             marker_class = col.names.info),
    all_cellType_order = c('CD4_CM','CD4_EM','CD4_Naive','CD4_TEMRA','CD4_unclass',
                           'CD8_CM','CD8_EM','CD8_Naive','CD8_TEMRA','CD8_unclass',
                           'B_activ','B_Mem','B_Naive','B_CD34_hi','B_trans','B_pl',
                           'NK_br_CD16_low','NK_br_CD16_hi','NK_dim_CD16_low','NK_dim_CD16_hi','NKT',
                           'cDC','pDC','Mono_classical','Mono_NonClassical'),
    colors_all = c('slategray1','slategray3','steelblue1','steelblue3','navy',
                   'magenta','hotpink2','darkorchid1','darkorchid3','darkmagenta',
                   'pink','firebrick1','firebrick3','indianred2','indianred3','red4',
                   'gold','orange','orange3','peachpuff','cyan3',
                   'darkslategrey','palegreen4','green','green3'),
    patientLevels_SMART = paste0("S", 1:18)
  )
}
