## Setup
# load packages

suppressPackageStartupMessages({
  library(magrittr)
  library(SpatialExperiment)
})

# Read SPE
# sce <- readRDS(here::here("data/ECM/results/spe/raw/SPE_combined_markedfilt.rds"))
sce <- readRDS(here::here(snakemake@input[["spe"]]))

# Apply filter
sce <- sce[,sce[["to_keep"]]]


# n_workers <- 1
n_workers <- snakemake@threads
RhpcBLASctl::blas_set_num_threads(n_workers)

# metadata_file <- here::here("data/ECM/config/sample_metadata.csv")
metadata_file <- snakemake@input[["sample_metadata"]]

metadf <- read.csv(metadata_file)

# cols_for_batch <- c("sample","slide","tissue")
cols_for_batch <- snakemake@params[["cols_for_batch"]]

# check
stopifnot(all(cols_for_batch %in% colnames(metadf)))
stopifnot("sample" %in% colnames(metadf))

# add metadata
colData(sce) <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::left_join(metadf,
                   by=c("sample_id"="sample")) %>%
  DataFrame(row.names=rownames(colData(sce)))
cols_for_batch[cols_for_batch=="sample"] <- "sample_id"

# convert NA to "missing", otherwise harmony will fail
colData(sce)[,cols_for_batch] <- colData(sce)[,cols_for_batch] %>%
  as.data.frame() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), \(x) ifelse(is.na(x),"missing", x))) %>%
  DataFrame(row.names=rownames(colData(sce)))

# PCA
set.seed(123)
sce <- scater::runPCA(sce, subset_row=(rowData(sce)$type == "type") & rowData(sce)$use_for_clustering, exprs_values="exprs")

# UMAP
set.seed(123)
sce <- scater::runUMAP(sce, subset_row=(rowData(sce)$type == "type") & rowData(sce)$use_for_clustering, exprs_values="exprs", n_threads=n_workers)

## Batch correction

# number of pc's, here use all - 1, but no more than 10
n_pcs <- min(sum(rowData(sce)$type=="type") - 1, 10)

# prepare for harmony and run
mat <- t(assay(sce, "exprs")[(rowData(sce)$type=="type") & rowData(sce)$use_for_clustering,])
vars_use <- colData(sce)[,cols_for_batch, drop=FALSE]
suppressWarnings(harmony_emb <- harmony::HarmonyMatrix(mat, vars_use, 
                                                       vars_use = cols_for_batch,
                                                       npcs=min(dim(mat)[2]-1,n_pcs), do_pca=TRUE,
                                                       max.iter.harmony=20,
                                                       max.iter.cluster=500,
                                                       block.size=0.025,
                                                       plot_convergence = FALSE,
                                                       verbose = FALSE
))
reducedDim(sce, "harmony") <- harmony_emb

set.seed(123)
sce <- scater::runUMAP(sce, dimred = "harmony", name = "UMAP_harmony", n_threads=n_workers)

# Save
saveRDS(sce, snakemake@output[["spe"]])
