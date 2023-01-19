## Setup
# load packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(future.apply)
})

# sce <- readRDS(here::here("data/ECM/results/spe/filt/SPE_combined_01.rds"))
sce <- readRDS(here::here(snakemake@input[["spe"]]))

# n_workers <- 1
n_workers <- snakemake@threads

# Clustering
# Run graph based clustering
# clustering with different numbers of shared neighbors (k) using walktrap and leiden.
resolutions <- c(0.01,0.05,0.1,0.15,0.2,0.25,0.5,1)
ksm <- list(leiden=rep(20,length(resolutions)),
            walktrap=c(20,25,30,35,40,45,50))
clu_algo <- c(rep("leiden",length(ksm[["leiden"]])),rep("walktrap",length(ksm[["walktrap"]])))
ks <- unlist(ksm)
clu_name <- paste0(clu_algo,"_k",ks)
clu_name[clu_algo=="leiden"] <- paste0(clu_name[clu_algo=="leiden"],"_res",resolutions)

cluster_args_ls <- lapply(seq_along(clu_name),function(i){
  if(stringr::str_detect(clu_name[i],"leiden")){
    list(resolution_parameter=as.numeric(stringr::str_replace(clu_name[i],"leiden_k20_res","")),n_iterations=4)
  } else{
    list()
  }
})


plan(multisession,workers=min(n_workers,length(ks)))

# options(future.globals.maxSize=+Inf)
sce_small <- sce
assay(sce_small,"counts") <- NULL
assay(sce_small,"exprs") <- NULL
colData(sce_small) <- NULL
reducedDim(sce_small,"PCA") <- NULL
reducedDim(sce_small,"UMAP") <- NULL
reducedDim(sce_small,"UMAP_harmony") <- NULL

set.seed(123)
clu_ls <- future_lapply(seq_along(clu_name), future.globals = list(sce_small=sce_small, clu_algo=clu_algo,ks=ks,cluster_args_ls=cluster_args_ls), future.seed=TRUE, function(i){
  suppressPackageStartupMessages(library(SpatialExperiment))
  scran::clusterCells(
    sce_small,
    use.dimred = "harmony",
    BLUSPARAM=bluster::TwoStepParam(
      first=bluster::KmeansParam(centers=ceiling(dim(sce_small)[2]^0.75)),
      second=bluster::SNNGraphParam(k=ks[i], type="rank", cluster.fun=clu_algo[i], cluster.args=cluster_args_ls[[i]])
    ), 
    full=TRUE)
})
rm(sce_small)
names(clu_ls) <- clu_name

for(c in clu_name){
  colData(sce)[[c]] <- clu_ls[[c]]$clusters
}

# Save
saveRDS(sce, snakemake@output[["spe"]])
