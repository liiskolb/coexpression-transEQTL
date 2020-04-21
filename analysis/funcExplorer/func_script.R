#!/usr/bin/env Rscript

# Script to create separate expression matrices that contain the metadata as well
# This is needed for the eQTL analysis

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-c", "--clres", help="Cluster result file (tsv)")
parser$add_argument("-s", "--savepath", help="Directory to save expressions data")
parser$add_argument("-d", "--datapath", help="Directory of metadata")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if (!is.null(args$c)) {
  if(is.null(args$savepath)){
    args$resdir = "."
  }else{
    if (!dir.exists(args$savepath)) { dir.create(args$savepath) }
  }
  print(paste("Clusters are in:", args$clres))
  print(paste("The results will be saved to:", args$savepath))
}else{
  args$print_help
  stop("At least expressions file should be provided", call.=FALSE)
}

#respath = "/gpfs/hpc/home/liiskolb/transqtl/results/integrated_coexpr/funcExplorer/clusters"
#savepath = "/gpfs/hpc/home/liiskolb/transqtl/results/integrated_coexpr/funcExplorer/expres"
#datapath = "/gpfs/hpc/home/liiskolb/transqtl/data"

datapath = args$datapath
resfile = args$clres
savepath = args$savepath

# funcexplorer results
res = read.table(resfile, sep="\t", header=T, stringsAsFactor=F, row.names = 1, check.names = F)
rownames(res) = paste0("Cluster_", rownames(res))

# Make sure to have Naranbhai sample names as CD16_1
colnames(res)[grepl("Sample ", colnames(res))] = paste0("CD16_", gsub("[^\\d]+", "", colnames(res)[grepl("Sample ", colnames(res))], perl=TRUE))

# Read in sample metadata

sample_metadata = readRDS(file.path(datapath, "metadata", "sample_metadata_filt.rds"))

# Read in sample groups
sample_list = readRDS(file.path(datapath, "metadata", "sample_list_by_group.rds"))

# Find the components: samples are in rows and clusters in columns
res_mat = t(res)

# Split eigen-expressions to separate groups and write them to a rds file that contains list of expression matrices
# Row is the sample_id and columns are the component indexes

#reslist = list()
for(name in names(sample_list)){
  idx = unlist(strsplit(sample_list[[name]],", "))
  idx = intersect(idx, colnames(res))
  if(length(idx) > 0){
    print(name)
    subdf = res_mat[as.character(idx),]
    subdf = base::merge(subdf, sample_metadata[,c("sample_id", "genotype_id", "study", "marker", "sex")], by.x = 0, by.y = "sample_id")
    rownames(subdf) = subdf$Row.names
    subdf$Row.names = NULL
    write.table(subdf, file.path(savepath, paste("funcExplorer", name, "exprs.tsv", sep="_")), quote = F, sep = "\t")
  }
  else{
    next
  }
}

print("Done!")
