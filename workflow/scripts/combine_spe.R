## Setup
# logging
stdlog <- file(snakemake@log[["stdout"]], open="wt")
sink(stdlog, type = "output")
stderr <- file(snakemake@log[["stderr"]], open="wt")
sink(stderr, type = "message")

# load packages

suppressPackageStartupMessages({
  library(magrittr)
  library(SpatialExperiment)
})


## Read SPE

#load spatialExperiment objects and combine
images_pres <- snakemake@input[["spes"]]

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

