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

library("SummarizedExperiment")
library("devtools")
#install_github("wgmao/PLIER")
library(PLIER)
library(dplyr)


respath = args$resdir
datapath = args$datapath
savepath = args$savepath

mat_norm = read.table(args$expr, header=TRUE, sep="\t", check.names=F, stringsAsFactors = F)

# load the pathway data for PLIER
allPaths = readRDS(file.path(respath, "PLIERpaths.rds"))

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
    
    # run PLIER
    cm.genes = commonRows(allPaths, subdf)
    
    print(paste0("Common genes: ", length(cm.genes)))
    
    print("Run svd")
    mat.svd = svd(subdf[cm.genes,])
    print("Run chat")
    Chat = computeChat(allPaths[cm.genes,])
    k = num.pc(mat.svd)
    # run PLIER
    print("Run PLIER")
    # genes in rows and samples in columns
    res = PLIER(data = subdf[cm.genes,], priorMat = allPaths[cm.genes,], k = k, doCrossval = FALSE, trace=TRUE, max.iter=100, svdres = mat.svd, allGenes = F, Chat=Chat, scale=T)
    print("PLIER is ready")
    
    # save the result
    group = name
    resfile = paste0("PLIER_", group, ".rds")
    saveRDS(res, file=file.path(respath, resfile))
    
    # Find the components: factors are in columns and samples in rows
    
    # Find the components: LVS are in rows and samples in columns
    res_mat = res$B
    res_mat = t(res_mat)
    # this is a hack to get correct sample ids
    row.names(res_mat) = gsub("\\.", " ", row.names(res_mat))
    
    
    res_mat = base::merge(res_mat, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
    rownames(res_mat) = res_mat$Row.names
    res_mat$Row.names = NULL
    write.table(res_mat, file.path(savepath, paste("PLIER", name, "separate_eigen.tsv", sep="_")), quote = F, sep = "\t")
  }
} else {
  print("Cell types together")
  # Remove zero-std genes
  std = apply(mat_norm, 1, sd)
  mat_norm = mat_norm[std!=0,]
  mat_norm = as.matrix(mat_norm)
  
  # run PLIER
  cm.genes = commonRows(allPaths, mat_norm)
  #nr.pc = num.pc(mat_scale)
  
  print(paste0("Common genes: ", length(cm.genes)))
  
  print("Run svd")
  mat.svd = svd(mat_norm[cm.genes,])

  print("Run chat")
  Chat = computeChat(allPaths[cm.genes,])
  k = num.pc(mat.svd)
  # run PLIER
  print("Run PLIER")
  res = PLIER(data = mat_norm[cm.genes,], priorMat = allPaths[cm.genes,], k = k, doCrossval = FALSE, trace=TRUE, max.iter=100, svdres = mat.svd, allGenes = F, Chat=Chat, scale=T)

  print("PLIER is ready")
  
  # save the result
  resfile = paste0("PLIER_integrated", ".rds")
  saveRDS(res, file=file.path(respath, resfile))
  
  # Script to create separate expression matrices that contain the metadata as well
  # This is needed for the eQTL analysis
  
  # Find the components: LVS are in rows and samples in columns
  res_mat = res$B
  res_mat = t(res_mat)
  # this is a hack to get correct sample ids
  row.names(res_mat) = gsub("\\.", " ", row.names(res_mat))
  
  # Split eigen-expressions to separate groups and write them to files that are expression matrices
  # Row is the sample_id and columns are the component indexes
  
  #reslist = list()
  for (name in names(sample_list2)){
    idx = sample_list2[[name]]
    if(any(idx %in% row.names(res_mat))){
      print(name)
      idx = intersect(idx, row.names(res_mat))
      subdf = res_mat[idx,]
      subdf = base::merge(subdf, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
      rownames(subdf) = subdf$Row.names
      subdf$Row.names = NULL
      write.table(subdf, file.path(savepath, paste("PLIER", name, "full_eigen.tsv", sep="_")), quote = F, sep = "\t")
    }else{
      next
    }
  }
}

print("Done!")