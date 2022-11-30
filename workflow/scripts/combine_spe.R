## Setup
# load packages

suppressPackageStartupMessages({
  library(magrittr)
  library(SpatialExperiment)
})

# load used images

samplenamefile <- snakemake@input[["samcsv"]]
images_spec <- read.csv(here::here(samplenamefile),header=FALSE)[[1]]

images_pres <- snakemake@input[["spes"]]
iminds <- sapply(images_spec, function(im) which(stringr::str_detect(images_pres,im) ))
images_pres <- images_pres[iminds]

steinbock_dir <- here::here("results")




## Read SPE

#load spatialExperiment objects and combine

se_ls <- purrr::map(images_pres, function(img_name){
  readRDS(here::here(img_name))
})

# only keep marker that are present in all samples
se_marker_intersect <- purrr::map(se_ls,~rownames(.x)) %>%
  purrr::reduce(intersect)
se_ls_int <- se_ls %>% 
  purrr::map(~.x[rownames(.x) %in% se_marker_intersect]) 

# only keep row data that are present in all samples
se_rownames_intersect <- purrr::map(se_ls,~names(rowData(.x))) %>%
  purrr::reduce(intersect)
for (i in seq_along(se_ls_int)) {
  rowData(se_ls_int[[i]]) <- rowData(se_ls_int[[i]])[,se_rownames_intersect]
  #rowData(se_ls_int[[i]])$type <- rowData(se_ls_int[[length(se_ls_int)]])$type
}

# combine
sce <-
  se_ls_int %>%
  do.call(cbind, .)

# rename cells
colnames(sce) <- paste0(sce$sample_id,"_",colnames(sce))

# save 
saveRDS(sce, here::here(snakemake@output[["spe"]]))

