#!/usr/bin/env Rscript

# Run trans-eqtl on meta lead snps and all genes

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs"

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-e", "--expr", help="Expressions file (tsv)")
parser$add_argument("-g", "--geno", help="Genotype file (GDS)") 
parser$add_argument("-c", "--cov", action="store_true", default=TRUE, help="Add if should include covariates")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


if (!is.null(args$expr) & !is.null(args$geno)) {
  print(paste("Expressions are in:", args$expr))
  print(paste("Genotypes are in:", args$geno))
  print(paste("Covariates usage:", args$cov))
  print(paste("The results will be saved to:", resdir))
}else{
  args$print_help
  stop("At least expressions file and genotypes should be provided", call.=FALSE)
}

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("GWASTools"))
suppressPackageStartupMessages(library("MatrixEQTL"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("gdsfmt"))
suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("reshape2"))


# Read in sample metadata
sample_metadata = read.table(file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data", "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)
print("sample metadata in")

# Genotype data is in the 'Genomic Data Structure' (GDS) format

gdsToMatrix <- function(gds_file){
  #Extract genotypes
  gds <- GWASTools::GdsGenotypeReader(gds_file)
  genotypes = GWASTools::getGenotype(gds)
  sample_ids = GWASTools::getVariable(gds, "sample.id")
  snp_rs_ids = GWASTools::getVariable(gds, "snp.rs.id")
  snp_ids = GWASTools::getVariable(gds, "snp.id")
  
  #Invent id for snps that do not have an rs id
  new_snp_ids = paste("snp",snp_ids[snp_rs_ids == ""], sep = "")
  snp_rs_ids[snp_rs_ids == ""] = new_snp_ids
  colnames(genotypes) = sample_ids
  rownames(genotypes) = snp_rs_ids
  
  #Extract SNP coordinates
  snpspos = dplyr::data_frame(snpid = snp_rs_ids, 
                              chr = GWASTools::getVariable(gds, "snp.chromosome"), 
                              pos = GWASTools::getVariable(gds, "snp.position"))
  GWASTools::close(gds)
  return(list(snpspos = snpspos, genotypes = genotypes))
}

# Run trans-MatrixEQTL

runMatrixEQTL <- function(exp_data, geno_data, covariates = NULL, pvOutputThreshold = 1e-5, model = modelLINEAR, resfile){
  # Run matrixeQTL on a prepared data set
  
  #Perform some sanity checks
  if(!all(colnames(exp_data) == colnames(geno_data))){
    print(colnames(exp_data))
    print(colnames(geno_data))
    stop("Column names of expression and genotype data are not equal.")
  }
  
  #Construct a SlicedData object of the expression data
  expression_sliced = SlicedData$new()
  expression_sliced$CreateFromMatrix(exp_data)
  expression_sliced$ResliceCombined(sliceSize = 2000)
  
  #Create a SlicedData object for the genotypes
  snps = SlicedData$new()
  snps$CreateFromMatrix(geno_data)
  snps$ResliceCombined()
  
  #Add covariates
  cvrt = SlicedData$new()
  
  if (!is.null(covariates)){
    if(!all(colnames(exp_data) == colnames(covariates))){
      print(colnames(exp_data))
      print(colnames(covariates))
      stop("Column names of expression and covariates data are not equal.")
    }
    cvrt$CreateFromMatrix(covariates)
    cvrt$ResliceCombined()
  }
  
  #RUN
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = expression_sliced,
    cvrt = cvrt,
    output_file_name = resfile,
    pvOutputThreshold = pvOutputThreshold,
    useModel = model,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  return(me)
}


# Helper functions

#' make a tabix file from a data.frame or data.table part of a 
#' methylKit object
#' 
#' @param df data part of a methylKit object, either data.frame or data.table
#' @param outfile name of the output file
#' 
#' @usage df2tabix(df,outfile)
#' @noRd
df2tabix <- function(df, outfile){
  
  if(file.exists(outfile)){
    message("overwriting ", outfile)
    unlink(outfile)
  }
  
  # write the file to disk
  message("writing file to disk")
  numCores <- 8
  print(paste("number of cores:", numCores))
  data.table::fwrite(df, file = outfile, quote = F, sep = "\t", col.names = F, nThread = numCores, verbose = TRUE, showProgress = TRUE)
  #readr::write_delim(df, path = outfile, delim="\t", col_names = F)
  #write.table(df, gzfile(outfile), quote=FALSE, sep="\t",
  #            col.names=FALSE, row.names=FALSE)
  
  # compressing
  message("compressing the file with bgzip...")
  command = paste("bgzip", outfile)
  system(command, wait = T)
  rm(df)
  gc(TRUE)
  # make tabix out if the file
  message("making tabix index...")
  zipped = paste0(outfile, ".gz")
  idx_file = Rsamtools::indexTabix(zipped,
                                   seq=1, start=2, end=3,
                                   skip=0, comment="#", zeroBased=FALSE)
  
  # delete initial tab file
  #message("removing file ..")
  #file.remove(outfile)
  return(zipped)
}

# Quantile normalization
tRank <- function(x) {
  t1 <- x
  t2 <- x[!is.na(x)]
  r <- qnorm(rank(t2)/(length(rank(t2))+1))
  t1[!is.na(t1)] <- r
  t1
}


# Wrapper to run MatrixEQTL on all the expression groups

eQTL_wrapper <- function(gen_file, exprs_file, snps = NULL, cov_data = TRUE, resdir){
  # Get genotypes data
  genotypes = gdsToMatrix(gen_file)
  gen_data = genotypes$genotypes # Matrix with full genotype data (id-s in columns)
  print("Genotypes are in")
  print(paste("Dimensions are", ncol(gen_data)))
  
  # Get expressions data
  exprs_data = read.table(exprs_file, sep="\t", header = TRUE) 
  # remember sample_id
  keep_ids = row.names(exprs_data)
  # use genotype id as row id
  rownames(exprs_data) = exprs_data$genotype_id
  
  # get the group name of expr data
  group_name = sub('\\.tsv$', '', basename(exprs_file))
  print(paste("The group is:", group_name))
  
  # expr mat
  exprs = exprs_data[,!colnames(exprs_data) %in% c("genotype_id", "sample_id")]
  exprs = t(exprs) # genotype ids in columns
  print("Expressions read in")
  print(paste("Dimensions: ", dim(exprs)))
  
  # Make sure the genotype ids match everywhere
  
  sample_ids = intersect(colnames(exprs), colnames(gen_data)) # genotypes ids
  
  if(length(sample_ids)==0){
    stop("Sample ids don't match at all")
  }else{
    print(paste("Number of mutual genotype ids:", length(sample_ids)))
    # keep data for overlapping genotype ids
    exprs = exprs[,sample_ids]
    
    # Quantile normalize eigenvectors
    #print("Quantile normalization")
    #exprs = t(data.frame(apply(exprs, 1, tRank))) # genotype IDs in columns
    #print("Quantile done")
    # Standardize gene profiles
    print("Standardize gene profiles")
    exprs = t(data.frame(apply(exprs, 1, tRank))) # genotype IDs in columns
    #exprs = t(scale(t(exprs)))
    print("Standardization done")
    
    gens = gen_data[,sample_ids]
    print(paste("Genotypes", dim(gens)))
    print(class(gens))
  }
  
  # Get covariates data
  
  if (cov_data == TRUE){
    print("Detecting covariates data")
    # get covariates from sample_metadata
    covs = sample_metadata %>% filter(genotype_id %in% sample_ids, sample_id %in% keep_ids) %>% select(genotype_id, batch, sex, marker)
    row.names(covs) = covs$genotype_id
    covs$genotype_id = NULL
    covs = covs[sample_ids,]
    
    # Exclude covariates that don't vary
    not_consts = apply(covs, 2, function(col) { length(unique(col)) > 1 })
    covs = covs[ ,not_consts, drop=FALSE]
    
    # Factors to numeric
    indx <- names(sapply(covs, is.factor))
    covs[indx] <- lapply(covs[indx], function(x) as.integer(x))
    
    # Detect colinear covariates and keep only one of them
    cors = cor(covs)
    cors[lower.tri(cors, diag=TRUE)] = NA
    is_col = base::subset(melt(cors, na.rm = TRUE), round(value) == 1)
    
    if (dim(is_col)[1] > 0){
      # For every row in here choose name of Var2
      col_vars = is_col$Var1
      not_col = setdiff(c("batch", "sex", "marker"), col_vars)
      covs = covs[, not_col, drop=FALSE]
    }
    
    # Run PCA on genotypes
    print("PCA on genotypes")
    # Run PCA
    gen_obj = snpgdsOpen(gen_file)
    # Try different LD thresholds for sensitivity analysis
    snpset <- SNPRelate::snpgdsLDpruning(gen_obj, sample.id = sample_ids, ld.threshold=0.3)
    # Get all selected snp id
    snpset.id <- unlist(unname(snpset))
    pca <- SNPRelate::snpgdsPCA(gen_obj, sample.id = sample_ids, snp.id = snpset.id, num.thread=3)
    # make a data.frame of first 3 PCs
    pcs <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],    # the third eigenvector
                      stringsAsFactors = FALSE)
    covs = base::merge(covs, pcs, by.x = 0, by.y = "sample.id")
    print(paste("Covariates dim: ", paste0(dim(covs), collapse=",")))
    rownames(covs) = covs$Row.names
    covs$Row.names = NULL
    covs = t(data.matrix(covs))
    covs = covs[,sample_ids]
    # remove rows with any NA-s
    covs = covs[complete.cases(covs),]
    print(paste("Covariates dim: ", paste0(dim(covs), collapse=",")))
  }else{
    covs = NULL
  }
  
  resfile = file.path(resdir, paste0(group_name, "_eQTLres"))
  print(paste("The results will be in the file:", resfile))
  
  # filter meta_lead snps from genotypes data
  gens = gens[snps,]
  print(paste("Genotypes dimensions are:", dim(gens)))
  
  # Run the trans-eQTL analysis
  me = runMatrixEQTL(exprs, gens, covs, pvOutputThreshold = 1, model = modelLINEAR, resfile = NULL)
  # Save file as gzip
  # add chr and location columns and order by that
  me_res = me$all$eqtls
  print("Clear memory")
  rm(me)
  rm(covs)
  rm(cors)
  gc(TRUE)
  print(dim(me_res))
  me_res$chr = as.numeric(sub("chr", "", sub("\\_.*", "", me_res$snps)))
  print("chr")
  me_res$start = as.integer(gsub("(.*_){1}(\\d+)_.+", "\\2", as.character(me_res$snps)))
  print("start")
  me_res$end = me_res$start
  print("end")
  me_res = me_res %>% dplyr::arrange(chr, start, end) %>% dplyr::select(chr, start, end, snps, gene, statistic, pvalue, FDR, beta) 
  zipped = df2tabix(me_res, resfile)
  print(zipped)
  return(TRUE)
}


# Get the metacs lead ids
metares = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
meta_leads = unique(metares$meta_id)


# Run the analysis
res = eQTL_wrapper(gen_file = args$geno, exprs_file = args$expr, snps = meta_leads, cov_data = args$cov, resdir = resdir)
print("Everything is ready!")