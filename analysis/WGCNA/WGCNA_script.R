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
library(dplyr)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 5)

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

if( args$dosep ){
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
    
    # columns are genes and rows are samples
    wgcnares = blockwiseModules(t(subdf))
    print("WGCNA is done")
    # save the result
    group = name
    resfile = paste0("WGCNA_", group, ".rds")
    saveRDS(wgcnares, file=file.path(respath, resfile))
    
    # Find the components: samples are in rows and modules in columns
    res_mat = data.frame(wgcnares$MEs)
    sample_order = colnames(subdf)
    row.names(res_mat) = sample_order
    # exclude the grey module that contains genes that are not in any clusters
    res_mat = res_mat[,colnames(res_mat)!="MEgrey"]

    res_mat = base::merge(res_mat, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
    rownames(res_mat) = res_mat$Row.names
    res_mat$Row.names = NULL
    write.table(res_mat, file.path(savepath, paste("WGCNA", name, "separate_eigen.tsv", sep="_")), quote = F, sep = "\t")
  }
} else {
  print("Cell types together")
  # Remove zero-std genes
  std = apply(mat_norm, 1, sd)
  mat_norm = mat_norm[std!=0,]
  wgcnares = blockwiseModules(t(mat_norm))
  print("WGCNA is done")
  # save the result
  #group = sub('\\..*$', '', basename(args$expr))
  resfile = paste0("WGCNA_integrated", ".rds")
  saveRDS(wgcnares, file=file.path(respath, resfile))
  
  # Script to create separate expression matrices that contain the metadata as well
  # This is needed for the eQTL analysis
  
  # Find the components: samples are in rows and modules in columns
  res_mat = data.frame(wgcnares$MEs)
  sample_order = colnames(mat_norm)
  row.names(res_mat) = sample_order
  
  # exclude the grey module that contains genes that are not in any clusters
  res_mat = res_mat[,colnames(res_mat)!="MEgrey"]
  
  # Split eigen-expressions to separate groups and write them to files that are expression matrices
  # Row is the sample_id and columns are the component indexes
  
  for(name in names(sample_list2)){
    idx = sample_list2[[name]]
    if( any(idx %in% row.names(res_mat)) ){
      print(name)
      idx = intersect(idx, row.names(res_mat))
      subdf = res_mat[idx,]
      subdf = base::merge(subdf, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
      rownames(subdf) = subdf$Row.names
      subdf$Row.names = NULL
      write.table(subdf, file.path(savepath, paste("WGCNA", name, "full_eigen.tsv", sep="_")), quote = F, sep = "\t")
    }else{
      next
    }
  }
}

print("Done!")
