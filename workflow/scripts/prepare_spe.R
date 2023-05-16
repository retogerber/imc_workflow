# set options
options(warn=2)
# logging
stdlog <- file(snakemake@log[["stdout"]], open="wt")
sink(stdlog, type = "output")
stderr <- file(snakemake@log[["stderr"]], open="wt")
sink(stderr, type = "message")

suppressPackageStartupMessages({
  library(magrittr)
  library(SpatialExperiment)
})

# read in parameters
#steinbock_dir <- "."
steinbock_dir <- snakemake@params[["steinbock_dir"]]

#channel_metadata_file <- "channel_metadata.csv"
channel_metadata_file <- snakemake@input[["channel_metadata_file"]]


marker_df <- read.csv(channel_metadata_file)
stopifnot("channel" %in% colnames(marker_df))
stopifnot("type" %in% colnames(marker_df))
marker_df <- dplyr::select(marker_df, channel, type)

#intensity_file <- "ECM_TMA_12_005.csv"
intensity_file <- snakemake@input[["intensities"]]
img_name <- stringr::str_replace(basename(intensity_file),"\\.csv$","")
# img_name <- "JG_11-22-2021_808106-2_PNGase_Bottom_resume_5"
#img_base_name_filt <- stringr::str_replace_all(img_name, "_ps05|_seg1234|_mespp_id[[:digit:]]{2}","")
#img_base_name <- stringr::str_replace_all(img_base_name_filt, "_filt(_thre[[:digit:]]{1,2}){0,1}","")
#img_folder_name <- paste0("img",ifelse(stringr::str_detect(img_name,"_filt"),"_filt",""))
#panel_file <- here::here(getwd(),steinbock_dir,"temp_summary_panels",paste0(img_name,"_summary.csv"))
panel_file <- snakemake@input[["panel_file"]]


# create spatialexperiment
spe <- imcRtools::read_steinbock(steinbock_dir, 
                                 pattern = paste0(img_name,".csv$"), 
                                 panel_file = panel_file)
colnames(spe) <- spe$ObjectNumber
# subset to given sample
spe <- spe[,spe$sample_id==img_name]

# add row metadata
rowData(spe) <- rowData(spe) %>% 
  as.data.frame() %>% 
  dplyr::mutate(channel = as.character(channel)) %>%
  dplyr::left_join(marker_df, by="channel") %>% 
  dplyr::mutate(type=ifelse(is.na(type),"unknown",type),
                type=factor(type)) %>% 
  DataFrame(row.names = rownames(spe))
rowData(spe)$use_for_clustering <- rowData(spe)$type == "type"


# transform area
colData(spe)$sqrt_area <- sqrt(colData(spe)$area)

# transformation
assay(spe, "exprs") <- asinh(counts(spe)+1)


saveRDS(spe, snakemake@output[[1]])
# saveRDS(spe, here::here("output","raw","spe",paste0("SPE_raw_",img_name,".rds")))

