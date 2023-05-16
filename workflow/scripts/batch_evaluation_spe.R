## Setup
# logging
stdlog <- file(snakemake@log[["stdout"]], open="wt")
sink(stdlog, type = "output")
stderr <- file(snakemake@log[["stderr"]], open="wt")
sink(stderr, type = "message")

# load packages
suppressPackageStartupMessages({
  library(SpatialExperiment)
  #library(BiocParallel)
})


# load data
sce <- readRDS(here::here(snakemake@input[["spe"]]))

# number of threads
n_workers <- as.integer(snakemake@threads)
n_workers <- 1
print(n_workers)
#bpparam <- BiocParallel::MulticoreParam(workers=n_workers, RNGseed = 123)

# cols_for_batch <- c("sample","tissue")
cols_for_batch <- snakemake@params[["cols_for_batch"]]
cols_for_batch[cols_for_batch=="sample"] <- "sample_id"
stopifnot(all(cols_for_batch %in% colnames(colData(sce))))

k <- as.integer(snakemake@params[["k"]])

# calculate cell mixing score to see how good batch correction worked.
for(btc in cols_for_batch){
  set.seed(123)
  print(btc)
  #suppressWarnings(
    sce <- CellMixS::cms(sce, k=k, group = btc,
                         assay_name = "exprs", res_name = paste0(btc,"_uncorrected"),
                         unbalanced=TRUE)
  #)
  print(2) 
  #suppressWarnings(
    sce <- CellMixS::cms(sce, k=k, group = btc,
                         dim_red = "harmony", res_name = paste0(btc,"_harmony"),
                         n_dim=dim(reducedDim(sce,"harmony"))[2],
                         unbalanced=TRUE)
  #)
}

saveRDS(sce, here::here(snakemake@output[["spe"]]))

