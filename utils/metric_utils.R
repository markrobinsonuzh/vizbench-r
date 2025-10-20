
load_pkgs <- function() {
  library(SeuratWrappers)
  library(Seurat)
  library(irlba)
  library(readr)
# library(parallel)
# library(distances)
# library(cluster)
# library(Rfast)
# library(dplyr)
}

  
celltype_shape = function(args) {
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  val = sapply(unique(batch), function(b){
    tapply(which(batch==b), celltype[batch==b], function(id){
      if(ncol(data)>2){
        svd_res = irlba(apply(data[id,], 2, function(d){d-mean(d)}), nv=2)
      }else{
        svd_res = svd(apply(data[id,], 2, function(d){d-mean(d)}))
      }
      as.numeric((svd_res$d[2] / svd_res$d[1]) >= 0.25)
    })
  })
  return(mean(val,na.rm=TRUE))
}


batch_mixture <- function(args, seed=42, n.cores = 10){
  
  # read embeddings
  data <- read_csv(args$visualize.csv.gz)
  # read SCE to get batch/celltype 
  sce <- read_sce(args$integrate.ad)
  batch <- sce$batch
  celltype <- sce$celltype
  rm(sce)
  
  set.seed(seed)
  n <- nrow(sce)
  if(n < 100000) { B <- 1 } else { B <- 100 }
  if(B == 1) n.cores <- 1

  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = FALSE)
    res = mean(sapply(unique(celltype[id]), function(c){
      ids = id[celltype[id]==c]
      if(length(unique(batch[ids]))>1){
        per = min(round(table(batch[ids])/2),30)
        lisi_res  = compute_lisi(data[ids,], data.frame(batch=batch[ids]), 
                                 "batch", perplexity = per)
        return(mean(lisi_res$batch))
      }else{
        return(NA)
      }
    }
    ),na.rm=TRUE)
    return(res)
  }, mc.cores = n.cores)
  val = unlist(val)

  return(mean(val,na.rm=T))
}


# original
# CelltypeShape = function(data, celltype, batch, seed=42){
#   val = sapply(unique(batch), function(b){
#     tapply(which(batch==b), celltype[batch==b], function(id){
#       if(ncol(data)>2){
#         svd_res = irlba(apply(data[id,], 2, function(d){d-mean(d)}), nv=2)
#       }else{
#         svd_res = svd(apply(data[id,], 2, function(d){d-mean(d)}))
#       }
#       as.numeric((svd_res$d[2] / svd_res$d[1]) >= 0.25)
#     })
#   })
#   return(mean(val,na.rm=T))
# }

# original
# BatchMixture = function(data, celltype, batch, seed=42, B, n, n.cores){
#   set.seed(seed)
#   if(B == 1){
#     n.cores=1
#   }
#   val = mclapply(1:B, FUN=function(i){
#     id = sample(1:nrow(data), n, replace = F)
#     res = mean(sapply(unique(celltype[id]), function(c){
#       ids = id[celltype[id]==c]
#       if(length(unique(batch[ids]))>1){
#         per = min(round(table(batch[ids])/2),30)
#         lisi_res  = compute_lisi(data[ids,], data.frame(batch=batch[ids]), "batch", perplexity = per)
#         return(mean(lisi_res$batch))
#       }else{
#         return(NA)
#       }
#     }
#     ),na.rm=T)
#     return(res)
#   },mc.cores = n.cores)
#   val = unlist(val)
#   
#   return(mean(val,na.rm=T))
# }
# 


# might need to implement this logic to prep datasets
# so far, implementing those without mean/var params from simulation
# instead of recomputing, this should be calculated at simulation time and read here.
# for(data in Dataset){
#     simulation = readRDS(file.path(data_dir,"Data/Simulation",
#                                    data,paste0(dataset_name, "_simulation.rds")))
#     para = readRDS(file.path(data_dir,"Data/Simulation",
#                              data,paste0(dataset_name, "_para.rds")))
#     celltype = simulation$meta$celltype
#     batch = simulation$meta$batch
#     ls = log10(colSums(simulation$counts))
#     zp <- colMeans(simulation$counts==0)
#     rm(simulation)
#     idx = lapply(unique(batch), function(b) tapply(which(batch==b), celltype[batch==b], function(id) id[1]))
#     names(idx) = unique(batch)
#     mean_par = lapply(idx, function(idx) para$mean_mat[idx,])
#     mean_par = lapply(mean_par, function(obj) {rownames(obj) = names(idx[[1]]); obj})
#     
#     variance_par = lapply(idx, function(idx)  para$mean_mat[idx,]^2 * para$sigma_mat[idx,])
#     variance_par = lapply(variance_par, function(obj) {rownames(obj) = names(idx[[1]]); obj})
#     
#     rm(para)
#     for(norm in Normalization){
#       for(integ in Integration){
#         for(visual in Visualization){
#           res = try(readRDS(paste0(path_visualization_results,"/",data,"_",norm,"+",integ,"+",visual,".rds")))
#           for(met in Metric){
#             if(!inherits(res, "try-error")){
#               val = Metric_eval(data =res, 
#                                 mean_par = mean_par, variance_par = variance_par,
#                                 library_size = ls, zero_proportion = zp,
#                                 celltype = celltype, batch = batch,
#                                 metric = met, n.cores = 20) 
#               metric_df[k, ] = c(data, norm, integ, visual, met, val)
#             }else{
#               metric_df[k, ] = c(data, norm, integ, visual, met, NA)
#             }
#             print(k)
#             saveRDS(metric_df,paste0(path_save_evaluation,"/",data,"_evaluation_metrics.rds"))
#             k=k+1
#           }
#         }
#       }
#     }
#   }


  

# Metric_eval = function(data, mean_par, variance_par, library_size, zero_proportion, celltype, batch, metric,
#                        seed = 42, B = 100, n = 10000, n.cores=10){
#   if(nrow(data)<100000){
#     B = 1
#     n = nrow(data)
#   }
#   val <- switch( metric,
#                  "DistancePreservation" = DistancePreservation(data, mean_par, celltype, batch, seed=seed),
#                  "VariancePreservation" = VariancePreservation(data, variance_par, celltype, batch, seed=seed),
#                  "VarianceSampleSize" = VarianceSampleSize(data, variance_par, celltype, batch, seed=seed),
#                  "LibrarySize" = LibrarySize(data, library_size, celltype, batch, seed=seed, B=B, n=n,n.cores=n.cores),
#                  "ZeroProportion" = ZeroProportion(data, zero_proportion, celltype, batch, seed=seed, B=B, n=n,n.cores=n.cores), 
#                  "CelltypeShape" = CelltypeShape(data, celltype, batch, seed=seed),
#                  "CelltypeSeparation" = CelltypeSeparation(data, celltype, batch, seed=seed, B=B, n=n,n.cores=n.cores),
#                  "BatchMixture" = BatchMixture(data,  celltype, batch, seed=seed, B=B, n=n,n.cores=n.cores),
#                  "Invalid Input!")
#   return(val)
# }
  

# distance_preservation <- function(data, mean_par, celltype, batch, seed=42, dist = "pca"){
#     val = sapply(unique(batch), function(b){
#       h = mean_par[[b]]
#       l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], mean))
#       h = h[complete.cases(h),]
#       l = l[complete.cases(l),]
#       if(dist == "pca"){
#         v = apply(h, 2, sd)
#         h = h[,v!=0]
#         pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
#         pca_scores <- pca_res$x
#         dh = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
#       }else{
#         dh = as.vector(distances::distance_matrix(distances::distances(h)))
#       }
#       if(dist == "both-pca"){
#         v = apply(l, 2, sd)
#         l = l[,v!=0]
#         pca_res <- prcomp(l, center = TRUE, scale. = TRUE)
#         pca_scores <- pca_res$x
#         dl = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
#       }else{
#         dl = as.vector(distances::distance_matrix(distances::distances(l))) 
#       }
#       cor(dh, dl, method = "spearman")
#     })
#     return(mean(val))
#   }

# 
# 
# DistancePreservation = function(data, mean_par, celltype, batch, seed=42, dist = "pca"){
#   val = sapply(unique(batch), function(b){
#     h = mean_par[[b]]
#     l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], mean))
#     h = h[complete.cases(h),]
#     l = l[complete.cases(l),]
#     if(dist == "pca"){
#       v = apply(h, 2, sd)
#       h = h[,v!=0]
#       pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
#       pca_scores <- pca_res$x
#       dh = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
#     }else{
#       dh = as.vector(distances::distance_matrix(distances::distances(h)))
#     }
#     if(dist == "both-pca"){
#       v = apply(l, 2, sd)
#       l = l[,v!=0]
#       pca_res <- prcomp(l, center = TRUE, scale. = TRUE)
#       pca_scores <- pca_res$x
#       dl = as.vector(distances::distance_matrix(distances::distances(pca_scores)))
#     }else{
#       dl = as.vector(distances::distance_matrix(distances::distances(l))) 
#     }
#     cor(dh, dl, method = "spearman")
#   })
#   return(mean(val))
# }

# VariancePreservation = function(data, variance_par, celltype, batch, seed=42, dist = "raw"){
#   val = sapply(unique(batch), function(b){
#     
#     h = mean_par[[b]]
#     l = apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var))
#     h = h[complete.cases(h),]
#     l = l[complete.cases(l),]
#     if(dist=="pca"){
#       v = apply(h, 2, sd)
#       h = h[,v!=0]
#       pca_res <- prcomp(h, center = TRUE, scale. = TRUE)
#       pca_scores <- pca_res$x
#       hv = apply(pca_scores, 1, var)
#     }else{
#       hv = rowSums(h)
#     }
#     if(dist=="both-pca"){
#       v = apply(l, 2, sd)
#       l = l[,v!=0]
#       pca_res <- prcomp(l, center = TRUE, scale. = TRUE)
#       pca_scores <- pca_res$x
#       lv = apply(pca_scores, 1, var)
#     }else{
#       lv = rowSums(l)
#     }
#     cor(hv, lv, method = "spearman")
#   })
#   return(mean(val))
# }
# 
# VarianceSampleSize = function(data, variance_par, celltype, batch, seed=42){
#   val = sapply(unique(batch), function(b){
#     h = tapply(which(batch==b), celltype[batch==b], length)
#     l = rowSums(apply(data[batch==b,],2, function(d) tapply(d, celltype[batch==b], var)))
#     h = h[complete.cases(h)]
#     l = l[complete.cases(l)]
#     1 - abs(cor(h, l, method = "spearman"))
#   })
#   return(mean(val))
# }
# 
# LibrarySize = function(data, library_size, celltype, batch, seed=42, B, n, n.cores){
#   set.seed(seed)
#   if(B == 1){
#     n.cores=1
#   }
#   val = mclapply(1:B, FUN=function(i){
#     id = sample(1:nrow(data), n, replace = F)
#     res = mean(sapply(unique(batch[id]), function(b){
#       ids = id[batch[id]==b]
#       if(length(unique(celltype[ids]))>1){
#         return(1 - tapply(ids, celltype[ids], function(idd) dcor(library_size[idd], data[idd,])$dcor))
#       }else{
#         return(NA)
#       }
#     }),na.rm=T)
#     return(res)
#   },mc.cores = n.cores)
#   val = unlist(val)
#   
#   return(mean(val))
# }
# 
# ZeroProportion = function(data, zero_proportion, celltype, batch, seed=42, B, n, n.cores){
#   set.seed(seed)
#   if(B == 1){
#     n.cores=1
#   }
#   val = mclapply(1:B, FUN=function(i){
#     id = sample(1:nrow(data), n, replace = F)
#     res = mean(sapply(unique(batch[id]), function(b){
#       ids = id[batch[id]==b]
#       if(length(unique(celltype[ids]))>1){
#         return(1 - tapply(ids, celltype[ids], function(idd) Rfast::dcor(zero_proportion[idd], data[idd,])$dcor))
#       }else{
#         return(NA)
#       }
#     }),na.rm=T)
#     return(res)
#   },mc.cores = n.cores)
#   val = unlist(val)
#   
#   return(mean(val))
# }
# 
# 
# CelltypeSeparation = function(data, celltype, batch, seed=42, B, n, n.cores){
#   set.seed(seed)
#   if(B == 1){
#     n.cores=1
#   }
#   start = Sys.time()
#   val = mclapply(1:B, FUN=function(i){
#     id = sample(1:nrow(data), n, replace = F)
#     res = mean(sapply(unique(batch[id]), function(b){
#       ids = id[batch[id]==b]
#       if(length(unique(celltype[ids]))>1){
#         dist = distances(data[ids,])
#         sil_res = silhouette(as.numeric(celltype[ids]), dist)
#         return(mean(tapply(sil_res[,"sil_width"], celltype[ids], mean)))
#       }else{
#         return(NA)
#       }
#     }
#     ), na.rm=T)
#     return(res)
#   },mc.cores = n.cores)
#   end = Sys.time()
#   val = unlist(val)
#   return(mean(val,na.rm=T))
# }

batch_mixture <- function(data, celltype, batch, seed=42, B, n, n.cores){
  set.seed(seed)
  if(B == 1){
    n.cores=1
  }
  val = mclapply(1:B, FUN=function(i){
    id = sample(1:nrow(data), n, replace = F)
    res = mean(sapply(unique(celltype[id]), function(c){
      ids = id[celltype[id]==c]
      if(length(unique(batch[ids]))>1){
        per = min(round(table(batch[ids])/2),30)
        lisi_res  = compute_lisi(data[ids,], data.frame(batch=batch[ids]), "batch", perplexity = per)
        return(mean(lisi_res$batch))
      }else{
        return(NA)
      }
    }
    ),na.rm=T)
    return(res)
  },mc.cores = n.cores)
  val = unlist(val)
  
  return(mean(val,na.rm=T))
}


