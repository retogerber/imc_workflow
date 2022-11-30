## Setup
# load packages

suppressPackageStartupMessages({
  library(magrittr)
  library(SpatialExperiment)
})


sce <- readRDS(here::here(snakemake@input[["spe"]]))

masks <- snakemake@input[["masks"]]

# Quality control

## Flag border cells

samples <- unique(sce$sample_id)

iminds <- sapply(samples, function(im) which(stringr::str_detect(masks,im) ))
masks <- masks[iminds]
bcdf <- purrr::map(masks,function(msk){
  mask <- EBImage::readImage(here::here(msk))*(2^16-1)
  data.frame(sample_id=sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(msk)),
             border_cell_ids=unique(c(mask[1,],mask[dim(mask)[1],],mask[,1],mask[,dim(mask)[2]]))
             )
})
bcdf <- do.call(rbind,bcdf)
bcdf <- bcdf %>% 
  dplyr::filter(border_cell_ids!=0) %>% 
  dplyr::mutate(cell_id=paste0(sample_id,"_",border_cell_ids))


sce[["is_border_cell"]] <- (colnames(sce) %in% bcdf$cell_id)
sce[["to_keep"]] <- !sce[["is_border_cell"]] 
#table(sce[["to_keep"]], useNA="ifany")

## image size

colData(sce)$no_pixels <- colData(sce)$width_px * colData(sce)$height_px

sce[["is_small_image"]] <- colData(sce)$no_pixels <= 1e4
sce[["to_keep"]] <- sce[["to_keep"]] & !sce[["is_small_image"]] 
#table(sce[["is_small_image"]], useNA="ifany")

## cell size

# flag cells with NA area

sce[["is_na_area"]] <- is.na(colData(sce)$area)
sce[["to_keep"]] <- sce[["to_keep"]] & !sce[["is_na_area"]] 
#table(sce[["is_na_area"]], useNA="ifany")

# flag cells with area smaller than 10 pixels

sce[["is_small_cell"]] <- colData(sce)$area < snakemake@params[["min_cell_area"]]
sce[["to_keep"]] <- na.omit(sce[["to_keep"]] & !sce[["is_small_cell"]])
#table(sce[["is_small_cell"]], useNA="ifany")


## celltype markers

# flag cells with mean cell type marker expression == 0 (all are equal to zero)

sce[["is_zero_celltype_expression"]] <- colMeans(counts(sce)[rownames(sce)[rowData(sce)$type=="type"],])==0
sce[["to_keep"]] <- sce[["to_keep"]] & !sce[["is_zero_celltype_expression"]] 

#table(sce[["is_zero_celltype_expression"]], useNA="ifany")



# Flag outlier cells (marker intensity)

mat <- assay(sce,"exprs")
sce$is_outlier_expression <- FALSE
n_to_remove <- c()
inds_to_remove <- list()
for(i in seq_len(dim(mat)[1])){
  x <- mat[i,]
  x_sub <- x[x>quantile(x,0.999)]
  a <- median(x_sub)
  b <- mad(x_sub)
  inds_to_remove[[i]] <- x>a+5*b
  n_to_remove <- c(n_to_remove,sum(x_sub>a+5*b))
  sce$is_outlier_expression <- sce$is_outlier_expression | (x > a+5*b)
}
indscomb <- do.call(rbind,inds_to_remove)
tmp <- apply(indscomb,2,sum)

sce[["to_keep"]] <- sce[["to_keep"]] & !sce[["is_outlier_expression"]] 
#table(sce[["is_outlier_expression"]])

saveRDS(sce, here::here(snakemake@output[["spe"]]))

