load_pkgs <- function() {
  library(scDesign3)
  library(SingleCellExperiment)
  library(Seurat)
  library(pbmcapply)
  library(parallel)
}

## NOTE: code below uses 40 cores!!


scdesign3 <- function(args) {
  
  sce <- read_sce(args$rawdata.ad)
  data = counts(sce)
  coldat = colData(sce)
  
  # QC
  batch = unique(coldat$batch)
  meta = data.frame(cbind(celltype=as.character(coldat$celltype),
                          batch=as.character(coldat$batch)))
  batch_filtered = colnames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,2]]
  celltype_filtered = rownames(table(meta))[which(table(meta)<5&table(meta)>0, arr.ind = T)[,1]]
  filtered_genes = which(apply(data,1,function(x) length(which(x!=0)))<10)
  if(length(filtered_genes)>0){
    data = data[-filtered_genes,]  
  }
  if(length(batch_filtered) > 0){
    filter = paste0(batch_filtered,"_", celltype_filtered)
    idx = which(paste0(meta$batch,"_", meta$celltype) %in% filter)
    data = data[, -idx]
    meta = meta[-idx, ]
    coldat = coldat[-idx,]
  }
  
  # HVG 
  seurat = CreateSeuratObject(data)
  seurat = FindVariableFeatures(seurat,nfeatures = 1000) 
  # MR: 1000 to make it faster?
  hvg = VariableFeatures(seurat)
  
  data = data[hvg,]
  if(args$verbose) message(paste0(dim(data), collapse=","))
  ## scDesign3
  sce <- SingleCellExperiment(assay = list(counts = data), 
                              colData = coldat)

  set.seed(123)
  data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "celltype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = "batch",
    corr_by = "celltype",
    ncell = 30000
    #ncell = 100000
  )
  
  if(args$verbose) message(mean(data$count_mat>0))
  if(args$verbose) message(paste0(dim(data$count_mat), collapse=","))
  
  marginal <- fit_marginal(
    data = data,
    predictor = "gene",
    #mu_formula = "celltype+s(batch,bs='re')",
    mu_formula = "celltype+batch",
    sigma_formula = "celltype",
    family_use = "nb",
    n_cores = 40,
    usebam = FALSE,
    parallelization = "mcmapply",
    # parallelization = "pbmcmapply",
    trace = TRUE
  )
  
  copula <- fit_copula(
    sce = sce,
    assay_use = "counts",
    marginal_list = marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 40,
    input_data = data$dat
  )
  
  para <- extract_para(
    sce = sce,
    marginal_list = marginal,
    n_cores = 40,
    family_use = "nb",
    new_covariate = data$newCovariate,
    parallelization = "mcmapply",
    data = data$dat#,
    #parallelization = "pbmcmapply"
  )
  
  # should we save the parameters? atm, just saving the dataset
  # dir.create(file.path(base_dir, "Data/Simulation"))
  # dir.create(file.path(base_dir, "Data/Simulation/Pilot"))
  # saveRDS(para,file.path(base_dir, "Data/Simulation/Pilot/MousePancreas_para.rds"))

  if(args$verbose) message(paste0(names(data), collapse=","))
  
  set.seed(123)
  newcounts = simu_new(
    sce = sce,
    filtered_gene = data$filtered_gene,
    mean_mat = para$mean_mat,
    sigma_mat = para$sigma_mat,
    zero_mat = para$zero_mat,
    quantile_mat = NULL,
    copula_list = copula$copula_list,
    n_cores = 40,
    family_use = "nb",
    input_data = data$dat,
    new_covariate = data$newCovariate,
    important_feature = copula$important_feature,
    parallelization = "pbmcmapply"
  )

  seurat.obj = CreateSeuratObject(newcounts, 
				  meta.data=data$newCovariate)
  return(seurat.obj)
}

