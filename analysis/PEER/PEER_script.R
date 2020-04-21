#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-e", "--expr", help="Expressions file (tsv)")
parser$add_argument("-r", "--resdir", help="Directory to save the clustering results to")
parser$add_argument("-d", "--datapath", help="Directory of metadata")
parser$add_argument("-s", "--savepath", help="Directory to save expressions data")
parser$add_argument("-sep", "--dosep", help="Cluster cell types separately", action="store_true")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if (!is.null(args$expr)) {
  if (is.null(args$resdir)){
    args$resdir = "."
  } else {
    if (!dir.exists(args$resdir)) { dir.create(args$resdir) }
    if (!dir.exists(args$savepath)){ dir.create(args$savepath) }
  }
  print(paste("Expressions are in:", args$expr))
  print(paste("The results will be saved to:", args$resdir))
}else{
  args$print_help
  stop("At least expressions file should be provided", call.=FALSE)
}


library(dplyr)
library(PLIER)
library(peer)

respath = args$resdir
datapath = args$datapath
savepath = args$savepath

mat_norm = read.table(args$expr, header=TRUE, sep="\t", check.names=F, stringsAsFactors = F)

# Read in sample metadata
sample_metadata = read.table(file.path(datapath, "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)
row.names(sample_metadata) = sample_metadata$sample_id
# Create sample groups
sample_list = sample_metadata %>% select(sample_id, qtl_group) %>% group_split(qtl_group)
sample_list2 = lapply(sample_list, function(x){return(x$sample_id)})
names(sample_list2) = unlist(lapply(sample_list, function(x){return(x$qtl_group[1])}))

if(args$dosep){
  # split expressions by cell type
  print("Separate for cell types")
  for(name in names(sample_list2)){
    print(name)
    idx = sample_list2[[name]]
    idx = intersect(idx, colnames(mat_norm))
    subdf = mat_norm[,idx]
    # Remove zero-std genes
    std = apply(subdf, 1, sd)
    subdf = subdf[std!=0,]
    
    subdf = as.matrix(subdf)
    mat_scale = rowNorm(subdf)
    nr.pc = num.pc(mat_scale)
    
    mat_scale = as.matrix(t(mat_scale))
    
    # Create the model object
    print("Creating the model object")
    model = PEER()
    
    # set the observed data
    
    PEER_setPhenoMean(model, mat_scale)
    print(dim(PEER_getPhenoMean(model)))
    
    # we want to infer K hidden confounders
    PEER_setNk(model, nr.pc)
    
    # include covariate to account for the mean expression
    print("Account for mean expression")
    PEER_setAdd_mean(model, TRUE)
    
    # perform the inference
    print("Performing the inference")
    PEER_update(model)
    
    # observing output
    print("Get factors")
    factors = PEER_getX(model)
    
    print("Get weights")
    weights = PEER_getW(model)
    
    print("Get precision")
    precision = PEER_getAlpha(model)
    
    print("Get residuals")
    residuals = PEER_getResiduals(model)
    
    print("Merge the results")
    res = list()
    res[["factors"]] = factors
    res[["weights"]] = weights
    res[["precision"]] = precision
    res[["residuals"]] = residuals
    
    print("PEER is done")
    # save the result
    group = name
    resfile = paste0("PEER_", group, ".rds")
    saveRDS(res, file=file.path(respath, resfile))
    
    # Find the components: factors are in columns and samples in rows
    res_mat = data.frame(res$factors)
    
    # First column is the mean factor
    res_mat = res_mat[,-1]
    
    colnames(res_mat) = paste0("F", seq(1, ncol(res_mat)))
    
    sample_order = colnames(subdf)
    row.names(res_mat) = sample_order
    
    res_mat = base::merge(res_mat, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
    rownames(res_mat) = res_mat$Row.names
    res_mat$Row.names = NULL
    write.table(res_mat, file.path(savepath, paste("PEER", name, "separate_eigen.tsv", sep="_")), quote = F, sep = "\t")
  }
} else {
  print("Cell types together")
  # Remove zero-std genes
  std = apply(mat_norm, 1, sd)
  mat_norm = mat_norm[std!=0,]
  mat_norm = as.matrix(mat_norm)
  mat_scale = rowNorm(mat_norm)
  
  nr.pc = num.pc(mat_scale)
  
  mat_scale = as.matrix(t(mat_scale))
  
  # Create the model object
  print("Creating the model object")
  model = PEER()
  
  # set the observed data
  
  PEER_setPhenoMean(model, mat_scale)
  print(dim(PEER_getPhenoMean(model)))
  
  # we want to infer K hidden confounders
  PEER_setNk(model, nr.pc)
  
  # include covariate to account for the mean expression
  print("Account for mean expression")
  PEER_setAdd_mean(model, TRUE)
  
  # perform the inference
  print("Performing the inference")
  PEER_update(model)
  
  # observing output
  print("Get factors")
  factors = PEER_getX(model)
  
  print("Get weights")
  weights = PEER_getW(model)
  
  print("Get precision")
  precision = PEER_getAlpha(model)
  
  print("Get residuals")
  residuals = PEER_getResiduals(model)
  
  print("Merge the results")
  res = list()
  res[["factors"]] = factors
  res[["weights"]] = weights
  res[["precision"]] = precision
  res[["residuals"]] = residuals
  
  print("PEER is done")
  
  # save the result
  resfile = paste0("PEER_integrated",".rds")
  saveRDS(res, file=file.path(respath, resfile))
  
  # Script to create separate expression matrices that contain the metadata as well
  # This is needed for the eQTL analysis
  
  # Find the components: factors are in columns and samples in rows
  res_mat = data.frame(res$factors)
  
  # First column is the mean factor
  res_mat = res_mat[,-1]
  
  colnames(res_mat) = paste0("F", seq(1, ncol(res_mat)))
  
  sample_order = colnames(mat_norm)
  row.names(res_mat) = sample_order
  
  # Split eigen-expressions to separate groups and write them to files that are expression matrices
  # Row is the sample_id and columns are the component indexes
  
  for(name in names(sample_list2)){
    idx = sample_list2[[name]]
    if(any(idx %in% row.names(res_mat))){
      print(name)
      idx = intersect(idx, row.names(res_mat))
      subdf = res_mat[idx,]
      subdf = base::merge(subdf, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
      rownames(subdf) = subdf$Row.names
      subdf$Row.names = NULL
      write.table(subdf, file.path(savepath, paste("PEER", name, "full_eigen.tsv", sep="_")), quote = F, sep = "\t")
    }else{
      next
    }
  }
}

print("Done!")