## Setup
# load packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(BiocParallel)
  library(magrittr)
})

# load data
sce <- readRDS(here::here(snakemake@input[["spe"]]))

clustering_used <- snakemake@params[["clustering_used"]]
clustering_annotation_file <- snakemake@input[["clustering_annotation_file"]]

annot_df <- read.csv(clustering_annotation_file)
stopifnot("cluster" %in% colnames(annot_df))
stopifnot("celltype" %in% colnames(annot_df))

# as factor
annot_df[["celltype"]] <- factor(annot_df[["celltype"]])
annot_df[["cluster"]] <- factor(annot_df[["cluster"]])

# additional cluster column
sce[["celltype_cluster"]] <- factor(sce[[clustering_used]])

# add annotation to data
colData(sce) <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::left_join(annot_df,by=c("celltype_cluster"="cluster")) %>% 
  DataFrame(row.names = rownames(colData(sce)))

# Save
saveRDS(sce, snakemake@output[["spe"]])
