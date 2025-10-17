
load_pkgs <- function() {
  library(Seurat)
}

SeuratUMAP = function(args, n.pcs=20, n.cores = 10){
  message("Running SeuratUMAP")
  message(try(reticulate::use_python("/root/.virtualenvs/r-reticulate/bin/python")))
  z <- reticulate::py_config()
  message(z$python)
  message(z$pythonhome)
  message(z$version_string)
  seurat.obj <- read_seurat(args$integrate.ad)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:n.pcs, 
                       reduction="integrated",
                       umap.method = "umap-learn",
                       n_threads = n.cores)
  Embeddings(seurat.obj,"umap")
}
















