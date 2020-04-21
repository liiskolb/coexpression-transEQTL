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
library(picaplot)

respath = args$resdir
datapath = args$datapath
savepath = args$savepath

mat_norm = read.table(args$expr, header=TRUE, sep="\t", stringsAsFactors = F)

# Read in sample metadata
sample_metadata = read.table(file.path(datapath, "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)
row.names(sample_metadata) = sample_metadata$sample_id
# Create sample groups
sample_list = sample_metadata %>% select(sample_id, qtl_group) %>% group_split(qtl_group)
sample_list2 = lapply(sample_list, function(x){return(x$sample_id)})
names(sample_list2) = unlist(lapply(sample_list, function(x){return(x$qtl_group[1])}))

if (args$dosep){
  # split expressions by qtl group
  print("Separate for qtl groups")
  for (name in names(sample_list2)){
    print(name)
    idx = sample_list2[[name]]
    idx = intersect(idx, colnames(mat_norm))
    subdf = mat_norm[,idx]
    sample_info = sample_metadata[idx, c("genotype_id", "sex", "cell_type_label", "batch", "smoker", "age", "condition_label")]
    icares = runICA(subdf, n_runs = 15, n_cores = 5, scale_pheno = T, max_iter = 10, var_cutoff = 70) # run ICA
    #icares = detectClusters(icares)
    print("ICA is done")
    # save the result
    group = name
    resfile = paste0("ICA_", group, ".rds")
    saveRDS(icares, file=file.path(respath, resfile))
    # This is a matrix that contains the components
    # sample_ids are in the rows and component id-s in the columns
    
    #ic_covar_mx = picaplot::getCovariateMx(icares)
    
    ic_covar_mx = t(icares$A) # sample_ids are in the rows and component id-s in the columns
    ic_covar_mx = base::merge(ic_covar_mx, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
    rownames(ic_covar_mx) = ic_covar_mx$Row.names
    ic_covar_mx$Row.names = NULL
    write.table(ic_covar_mx, file.path(savepath, paste("ICA", name, "separate_eigen.tsv", sep="_")), quote = F, sep = "\t")
  }
} else {
  print("Cell types together")
  sample_info = sample_metadata[, c("genotype_id", "sex", "cell_type_label", "batch", "smoker", "age", "condition_label")]
  icares = runICA(mat_norm, n_runs = 15, n_cores = 5, max_iter = 20, scale_pheno = T, var_cutoff = 70) # run ICA
  #icares = detectClusters(icares)
  print("ICA is done")
  # save the result
  #group = sub('\\..*$', '', basename(args$expr))
  resfile = paste0("ICA_integrated", ".rds")
  saveRDS(icares, file=file.path(respath, resfile))
  
  # Script to create separate expression matrices that contain the metadata as well
  # This is needed for the eQTL analysis
  
  # This is a matrix that contains clustered components
  # sample_ids are in the rows and component id-s in the columns
  
  #ic_covar_mx = picaplot::getCovariateMx(icares)
  ic_covar_mx = t(icares$A) # sample_ids are in the rows and component id-s in the columns
  
  # Split eigen-expressions to separate groups and write them to files that are expression matrices
  # Row is the sample_id and columns are the component indexes
  
  for (name in names(sample_list2)){
    idx = sample_list2[[name]]
    if ( any(idx %in% row.names(ic_covar_mx)) ){
      print(name)
      idx = intersect(idx, row.names(ic_covar_mx))
      subdf = ic_covar_mx[idx,]
      subdf = base::merge(subdf, sample_metadata[,c("sample_id", "genotype_id")], by.x = 0, by.y = "sample_id")
      rownames(subdf) = subdf$Row.names
      subdf$Row.names = NULL
      write.table(subdf, file.path(savepath, paste("ICA", name, "full_eigen.tsv", sep="_")), quote = F, sep = "\t")
    }else{
      next
    }
  }
}

print("Done!")
