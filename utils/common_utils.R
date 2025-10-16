

require(anndataR)
read_sce <- function(f) read_h5ad(f, as = "SingleCellExperiment")
read_seurat <- function(f) read_h5ad(f, as = "Seurat")

require(rjson)
read_normmethod <- function(f) fromJSON(paste(readLines(f), 
                                              collapse=""))$normalize
