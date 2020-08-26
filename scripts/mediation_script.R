#!/usr/bin/env Rscript

library("dplyr")
library("SummarizedExperiment")
library("GWASTools")
library("Rsamtools")
library("data.table")
#library(gdsfmt)
library("SNPRelate")
library("reshape2")
#library("GMAC")
library(mediation)
#library("MatrixEQTL")
library("Matrix")

### Helper functions ####
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
      
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
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
  expression_sliced$fileOmitCharacters = "NA"
  expression_sliced$ResliceCombined(sliceSize = 2000)
  
  #Create a SlicedData object for the genotypes
  snps = SlicedData$new()
  snps$CreateFromMatrix(geno_data)
  snps$fileOmitCharacters = "NA"
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
    cvrt$fileOmitCharacters = "NA"
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
  print(2)
  
  return(me)
}

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


# Rank normalization
tRank <- function(x) {
  t1 <- x
  t2 <- x[!is.na(x)]
  r <- qnorm(rank(t2)/(length(rank(t2))+1))
  t1[!is.na(t1)] <- r
  t1
}
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

#####

# Mediation analysis for ARHGEF3, SLC39A8 and LYZ
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004818#s4

# 1. First run MatrixeQTL so that the cis gene expression is a covariate in the model
# and save the results to files
# 2. Estimate the differences in effects using Sobel test or something

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation"
geno = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.gds"
all_clusters_filt = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters_filt.rds")

# Read in sample metadata
sample_metadata = read.table(file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data", "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)
print("sample metadata in")

# ARHGEF3
trans_lead = "chr3_56815721_T_C"
cis_gene = "ENSG00000163947" # ARHGEF3
cl_ids = names(all_clusters_filt[[trans_lead]])
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
#expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr3_56815721_T_C") %>% select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs = data.frame()
# for(i in 1:nrow(expr_meta)){
#   approach = expr_meta[i,][["approach"]]
#   cl_method = expr_meta[i,][["cl_method"]]
#   cl_id = expr_meta[i,][["cl_id"]]
#   qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#   eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#   eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", 
#                         ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                         cl_method, "eigens", eigenfile)
#   cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#   cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#   cl_eigen$sample_id = row.names(cl_eigen)
#   if(nrow(exprs) == 0){
#     print(cl_id)
#     exprs = cl_eigen
#   }else{
#     exprs = merge(exprs, cl_eigen, by = c("sample_id", "genotype_id"))
#   }
# }

#write.table(exprs, file = file.path(resdir, "ARHGEF3_eigens.tsv"), sep = "\t", quote = F, row.names = F)
exprs = read.table(file = file.path(resdir, "ARHGEF3_eigens.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs) = exprs$sample_id

# get ARHGEF gene expression
gene_expr = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/CL_0000233_naive.tsv", sep = "\t", stringsAsFactors = F, header = T)
cis_gene_expr = t(gene_expr[cis_gene, , drop = F])

mediation_wrapper <- function(gen_file, exprs_data, snps = NULL, cis_gene_expr = NULL, cov_data = TRUE, resdir){
  # Get genotypes data
  genotypes = gdsToMatrix(gen_file)
  gen_data = genotypes$genotypes # Matrix with full genotype data (id-s in columns)
  print("Genotypes are in")
  print(paste("Dimensions are", ncol(gen_data)))
  
  # Get expressions data
  #exprs_data = read.table(exprs_file, sep="\t", header = TRUE) 
  # remember sample_id
  keep_ids = row.names(exprs_data)
  # use genotype id as row id
  rownames(exprs_data) = exprs_data$genotype_id
  
  # expr mat
  exprs = exprs_data[,!colnames(exprs_data) %in% c("genotype_id", "sample_id"), drop = F]
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
    exprs = exprs[,sample_ids, drop = F]
    
    # Quantile normalize eigenvectors
    print("Standardize gene profiles")
    exprs = data.frame(apply(exprs, 1, tRank)) # genotype IDs in rows
    
    gens = gen_data[,sample_ids]
    print(paste("Genotypes", dim(gens)))
    print(class(gens))
  }
  
  # Get covariates data
  
  if (cov_data == TRUE){
    print("Detecting covariates data")
    # get covariates from sample_metadata
    print(head(sample_metadata))
    covs = sample_metadata %>% dplyr::filter(genotype_id %in% sample_ids, sample_id %in% keep_ids) %>% dplyr::select(genotype_id, batch, sex, marker)
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
    # close the file
    snpgdsClose(gen_obj)
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
    
    # add cis gene as covariate
    #print("Including cis gene expression as covariate")
    # include genotype_ids
    #cis_gene_expr = merge(cis_gene_expr, unique(sample_metadata[sample_metadata$sample_id %in% keep_ids,c("sample_id", "genotype_id")]), by.x = 0, by.y = "sample_id", all.x = TRUE)
    #row.names(cis_gene_expr) = cis_gene_expr$genotype_id
    #cis_gene_expr$genotype_id = NULL
    #cis_gene_expr$Row.names = NULL
    #cis_gene_expr = cis_gene_expr[sample_ids,, drop = F]
    # trank normalize
    #cis_gene_expr = t(cis_gene_expr)
    #cis_gene_expr = data.frame(apply(cis_gene_expr, 1, tRank)) # sample IDs to rows
    
    #covs = base::merge(covs, cis_gene_expr, by.x = 0, by.y = 0)
    #rownames(covs) = covs$Row.names
    #covs$Row.names = NULL
    #print(paste("Covariates dim2: ", paste0(dim(covs), collapse=",")))
    
    covs = covs[sample_ids, , drop = F]
    covs = t(data.matrix(covs))
    
    # remove rows with any NA-s
    covs = covs[complete.cases(covs),,drop = F]
    covs = t(covs)
    print(head(covs))
    
  }else{
    covs = NULL
  }
  
  # cis gene expression as mediatior
  cis_gene_expr = merge(cis_gene_expr, unique(sample_metadata[sample_metadata$sample_id %in% keep_ids,c("sample_id", "genotype_id")]), by.x = 0, by.y = "sample_id", all.x = TRUE)
  row.names(cis_gene_expr) = cis_gene_expr$genotype_id
  cis_gene_expr$genotype_id = NULL
  cis_gene_expr$Row.names = NULL
  cis_gene_expr = cis_gene_expr[sample_ids,, drop = F]
  # trank normalize
  cis_gene_expr = t(cis_gene_expr)
  cis_gene_expr = data.frame(apply(cis_gene_expr, 1, tRank)) # genotype IDs in rows
  
  group_name = paste(colnames(cis_gene_expr), "mediation", sep = "_")
  resfile = file.path(resdir, paste0(group_name, "_eQTLres"))
  print(paste("The results will be in the file:", resfile))
  
  # filter snps from genotypes data
  gens = gens[snps,, drop = F]
  print(paste("Genotypes dimensions are:", dim(gens)))
  #print(head(gens))
  
  # Run the mediation analysis
  # add cis gene expression to expr matrix
  exprs = merge(exprs, cis_gene_expr, by.x = 0, by.y = 0)
  row.names(exprs) = exprs$Row.names
  exprs$Row.names = NULL
  
  # trios - matrix of selected trios indexes (row numbers) for mediation tests
  # Each row consists of the index (i.e., row number) of the eQTL in snp.dat.cis, the index of cis-gene transcript in exp.dat, and the index of trans-gene in exp.dat.
  #trios = data.frame("idx.cisSNP" = 1, "idx.cisGene" = which(row.names(exprs) %in% colnames(cis_gene_expr)),
  #                   "idx.transGene" = which(!(row.names(exprs) %in% colnames(cis_gene_expr))))
  #print(trios)
  #print(head(exprs))
  #me = gmap.ac.gpd(snp.dat = gens, fea.dat = exprs, known.conf = covs, trios.idx = trios, pc.num = 0)
  
  #me = GMAC::gmac(known.conf = covs, cov.pool = NULL, exp.dat = exprs, 
  #                snp.dat.cis = gens, trios.idx = trios, nominal.p = TRUE, nperm = 50)
  mediation_analysis_fun <- function(data, dependent_v_string, independent_v_string, mediator_string, covs = NULL) {  
    require(mediation)
    
    #data <- data[complete.cases(data[, c(dependent_v_string, independent_v_string, mediator_string)]),]  
    
    dependent_v <- data[, dependent_v_string]
    independent_v <- data[, independent_v_string]
    mediator <- data[, mediator_string]  
    
    if(is.null(covs)){
      # step1 (total effect of independent variable on the dependent variable) good to have significant, but doesn't
      # necessarily have to be
      # eval(as.formula(paste(dependent_v, '~', independent_v)))
      #fit_step1 <- lm(as.formula(paste(dependent_v_string, '~', independent_v_string)), data = data)
      fit_step1 <- lm(dependent_v ~ independent_v)  
      
      # step2 (effect of the independent variable on mediator) must be significant!
      # eval(as.formula(paste(mediator, '~', independent_v))
      #fit_step2 <- lm(as.formula(paste(mediator_string, '~', independent_v_string)), data = data)
      fit_step2 <- lm(mediator ~ independent_v)
      
      # step3 (effect of the mediator on dependent viariable while controlling for the independent variable)
      # eval(as.formula(paste(dependent_v, '~', independent_v, '+', mediator)))
      #fit_step3 <- lm(as.formula(paste(dependent_v_string, '~', independent_v_string, '+', mediator_string)), data = data)
      fit_step3 <- lm(dependent_v ~ independent_v + mediator)  
    }else{
      # add covariates
      #data = merge(data, covs, by = 0)
      #print(head(data))
      covnames = colnames(covs)
      print(covnames)
      for (i in covnames) {
        assign(i, covs[, i]) # covariates as variables
      }
      
      # step1 (total effect of independent variable on the dependent variable) good to have significant, but doesn't
      # necessarily have to be
      # eval(as.formula(paste(dependent_v, '~', independent_v)))
      #fit_step1 <- lm(as.formula(paste(dependent_v_string, '~', independent_v_string, "+", paste(covnames, collapse = "+"))), data = data)
      #fit_step1 <- lm(dependent_v ~ independent_v + covs)  
      fit_step1 <- lm(as.formula(paste("dependent_v", "~", "independent_v", "+", paste(covnames, collapse = "+")) )) 
      
      # step2 (effect of the independent variable on mediator) must be significant!
      # eval(as.formula(paste(mediator, '~', independent_v))
      #fit_step2 <- lm(as.formula(paste(mediator_string, '~', independent_v_string, "+", paste(covnames, collapse = "+"))), data = data)
      #fit_step2 <- lm(mediator ~ independent_v + covs)
      fit_step2 <- lm(as.formula(paste("mediator", "~", "independent_v", "+", paste(covnames, collapse = "+")) ))
      print("step2")
      
      # step3 (effect of the mediator on dependent viariable while controlling for the independent variable)
      # eval(as.formula(paste(dependent_v, '~', independent_v, '+', mediator)))
      #fit_step3 <- lm(as.formula(paste(dependent_v_string, '~', independent_v_string, '+', mediator_string, "+", paste(covnames, collapse = " + "))), data = data)
      print("step3")
      #fit_step3 <- lm(dependent_v ~ independent_v + mediator + covs)  
      fit_step3 <- lm(as.formula(paste("dependent_v", "~", "independent_v", "+", "mediator", "+", paste(covnames, collapse = "+"))))  
    }
   
    # step4
    results <- mediate(model.m = fit_step2,
                       model.y = fit_step3,
                       treat = "independent_v",
                       mediator = "mediator",
                       boot = T)  
    
    #print(paste(dependent_v_string, '=', independent_v_string, '+', mediator_string))  
    
    print(summary(results))  
    
    return(results)  
  }
  
  mediation_mat = merge(exprs, t(gens), by = 0) # has variables in columns
  row.names(mediation_mat) = mediation_mat$Row.names
  mediation_mat$Row.names = NULL
  print(head(mediation_mat))
  covnames = colnames(covs)
  
  mediation_mat = merge(mediation_mat, covs, by = 0)
  
  mediator_cis = colnames(cis_gene_expr)
  med_res = data.frame()
  # linear regressions
  for(cl in colnames(exprs)){
    if(cl %in% c(mediator_cis, snps)){
      next
    }else{
      print(cl)
      mediator = mediation_mat[,mediator_cis]
      snp_v = mediation_mat[,snps]
      trans_v = mediation_mat[,cl]
      EV1 = covs[,"EV1"]
      EV2 = covs[,"EV2"]
      EV3 = covs[,"EV3"]
      print("step2")
      if("sex" %in% colnames(covs)){
        print("sex included")
        sex = covs[,"sex"]
        fit_step2 <- lm(mediator ~ snp_v + EV1 + EV2 + EV3 + sex)
        fit_step3 = lm(trans_v ~ snp_v + mediator + EV1 + EV2 + EV3 + sex)
      }else {
        fit_step2 <- lm(mediator ~ snp_v + EV1 + EV2 + EV3)
        fit_step3 = lm(trans_v ~ snp_v + mediator + EV1 + EV2 + EV3)
      }
      print("step3")
      me <- mediate(model.m = fit_step2,
                         model.y = fit_step3,
                         treat = "snp_v",
                         mediator = "mediator",
                         boot = T) 
      print("mediation done")
      #me = mediation_analysis_fun(mediation_mat, dependent_v_string = cl, independent_v_string = snps, mediator_string = mediator_cis, covs = covs)
      print(summary(me))
      mediation_res_summary <- me   
      subres = data.frame("snp" = snps, "trans" = cl, "cis" = mediator_cis, 
                          "ACME_estimate" = mediation_res_summary$d0,
                          "ACME_95CI" = paste(formatC(mediation_res_summary$d0.ci, digits = 2, format = 'e'), collapse = '...'),
                          "ACME_p_val" = mediation_res_summary$d0.p,      
                          
                          "ADE_estimate" = mediation_res_summary$z0,
                          "ADE_95CI" = paste(formatC(mediation_res_summary$z0.ci, digits = 2, format = 'e'), collapse = '...'),
                          "ADE_p_val" = mediation_res_summary$z0.p,      
                          
                          "total_effect_estimate" = mediation_res_summary$tau.coef,
                          "total_effect_95CI" = paste(formatC(mediation_res_summary$tau.ci, digits = 2, format = 'e'), collapse = '...'),
                          "total_effect_p_val" = mediation_res_summary$tau.p,      
                          
                          "prop_mediated_estimate" = mediation_res_summary$n0,
                          "prop_mediated_95CI" = paste(formatC(mediation_res_summary$n0.ci, digits = 2, format = 'e'), collapse = '...'),
                          "prop_mediated_p_val" = mediation_res_summary$n0.p,      
                          
                          "sample_size" = mediation_res_summary$nobs,
                          "simulations" = mediation_res_summary$sims)
      med_res = rbind(med_res, subres)
    }
  }
  #write.table(med_res, file = resfile, sep = "\t", quote = F, row.names = F)
  
  ###### old code #####
  #me = runMatrixEQTL(exprs, gens, covs, pvOutputThreshold = 1, model = modelLINEAR, resfile = NULL)
  # Save file as gzip
  # add chr and location columns and order by that
  #me_res = me$all$eqtls
  #print("Clear memory")
  #rm(me)
  #rm(covs)
  #rm(cors)
  #gc(TRUE)
  #print(dim(me_res))
 # me_res$chr = as.numeric(sub("chr", "", sub("\\_.*", "", me_res$snps)))
  #print("chr")
  #me_res$start = as.integer(gsub("(.*_){1}(\\d+)_.+", "\\2", as.character(me_res$snps)))
  #print("start")
  #me_res$end = me_res$start
  #print("end")
  #me_res = me_res %>% dplyr::arrange(chr, start, end) %>% dplyr::select(chr, start, end, snps, gene, statistic, pvalue, FDR, beta) 
  #zipped = df2tabix(me_res, resfile)
  #print(zipped)
  ########
  return(med_res)
}

# Run the analysis for ARHGEF3
print("Start mediation")
#res = mediation_wrapper(gen_file = geno, exprs_data = exprs, snps = trans_lead, cis_gene_expr = cis_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")
#write.table(res, file = file.path(resdir, paste0(cis_gene, "_mediation.tsv")), sep = "\t", quote = F, row.names = F)

# LYZ (integrated, monocytes naive)
trans_lead = "chr12_69344099_A_G"
cis_gene = "ENSG00000090382" # LYZ
cl_ids = names(all_clusters_filt[[trans_lead]])
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr12_69344099_A_G", qtl_group == "CL_0002057_naive", approach == "integrated") %>% dplyr::select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs_lyz = data.frame()
# for(i in 1:nrow(expr_meta)){
#    approach = expr_meta[i,][["approach"]]
#    cl_method = expr_meta[i,][["cl_method"]]
#    cl_id = expr_meta[i,][["cl_id"]]
#    qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#    eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#    eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", 
#                          ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                          cl_method, "eigens", eigenfile)
#    cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#    cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#    cl_eigen$sample_id = row.names(cl_eigen)
#    if(nrow(exprs_lyz) == 0){
#      print(cl_id)
#      exprs_lyz = cl_eigen
#    }else{
#      exprs_lyz = merge(exprs_lyz, cl_eigen, by = c("sample_id", "genotype_id"))
#    }
# }
# 
# write.table(exprs_lyz, file = file.path(resdir, "LYZ_eigens_integrated.tsv"), sep = "\t", quote = F, row.names = F)

# separate

# expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr12_69344099_A_G", qtl_group == "CL_0002057_naive", approach == "separate") %>% dplyr::select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs_lyz = data.frame()
# for(i in 1:nrow(expr_meta)){
#   approach = expr_meta[i,][["approach"]]
#   cl_method = expr_meta[i,][["cl_method"]]
#   cl_id = expr_meta[i,][["cl_id"]]
#   qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#   eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#   eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", 
#                         ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                         cl_method, "eigens", eigenfile)
#   cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#   cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#   cl_eigen$sample_id = row.names(cl_eigen)
#   if(nrow(exprs_lyz) == 0){
#     print(cl_id)
#     exprs_lyz = cl_eigen
#   }else{
#     exprs_lyz = merge(exprs_lyz, cl_eigen, by = c("sample_id", "genotype_id"))
#   }
# }
# 
# write.table(exprs_lyz, file = file.path(resdir, "LYZ_eigens_separate.tsv"), sep = "\t", quote = F, row.names = F)


# get LYZ gene expression in naive
gene_expr_lyz = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/CL_0002057_naive.tsv", sep = "\t", stringsAsFactors = F, header = T)
lyz_gene_expr = t(gene_expr_lyz[cis_gene, , drop = F])

exprs_lyz = read.table(file = file.path(resdir, "LYZ_eigens_integrated.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs_lyz) = exprs_lyz$sample_id

# Run the analysis for LYZ
print("Start mediation")
#res1 = mediation_wrapper(gen_file = geno, exprs_data = exprs_lyz, snps = trans_lead, cis_gene_expr = lyz_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")

exprs_lyz = read.table(file = file.path(resdir, "LYZ_eigens_separate.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs_lyz) = exprs_lyz$sample_id

# Run the analysis for LYZ
print("Start mediation")
#res2 = mediation_wrapper(gen_file = geno, exprs_data = exprs_lyz, snps = trans_lead, cis_gene_expr = lyz_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")

#write.table(res1, file = file.path(resdir, paste0(cis_gene, "_integrated_mediation.tsv")), sep = "\t", quote = F, row.names = F)
#write.table(res2, file = file.path(resdir, paste0(cis_gene, "_separate_mediation.tsv")), sep = "\t", quote = F, row.names = F)


# LYZ (integrated, monocytes IFNg 24h)
trans_lead = "chr12_69344099_A_G"
cis_gene = "ENSG00000090382" # LYZ
cl_ids = names(all_clusters_filt[[trans_lead]])
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr12_69344099_A_G", qtl_group == "CL_0002057_IFNg_24h", approach == "integrated") %>% dplyr::select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs_lyz = data.frame()
# for(i in 1:nrow(expr_meta)){
#    approach = expr_meta[i,][["approach"]]
#    cl_method = expr_meta[i,][["cl_method"]]
#    cl_id = expr_meta[i,][["cl_id"]]
#    qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#    eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#    eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results",
#                          ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                          cl_method, "eigens", eigenfile)
#    cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#    cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#    cl_eigen$sample_id = row.names(cl_eigen)
#    if(nrow(exprs_lyz) == 0){
#      print(cl_id)
#      exprs_lyz = cl_eigen
#    }else{
#      exprs_lyz = merge(exprs_lyz, cl_eigen, by = c("sample_id", "genotype_id"))
#    }
# }
# 
# write.table(exprs_lyz, file = file.path(resdir, "LYZ_ifng_eigens_integrated.tsv"), sep = "\t", quote = F, row.names = F)

# separate

expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr12_69344099_A_G", qtl_group == "CL_0002057_IFNg_24h", approach == "separate") %>% dplyr::select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs_lyz = data.frame()
# for(i in 1:nrow(expr_meta)){
#   approach = expr_meta[i,][["approach"]]
#   cl_method = expr_meta[i,][["cl_method"]]
#   cl_id = expr_meta[i,][["cl_id"]]
#   qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#   eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#   eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results",
#                         ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                         cl_method, "eigens", eigenfile)
#   cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#   cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#   cl_eigen$sample_id = row.names(cl_eigen)
#   if(nrow(exprs_lyz) == 0){
#     print(cl_id)
#     exprs_lyz = cl_eigen
#   }else{
#     exprs_lyz = merge(exprs_lyz, cl_eigen, by = c("sample_id", "genotype_id"))
#   }
# }
# 
# write.table(exprs_lyz, file = file.path(resdir, "LYZ_ifng_eigens_separate.tsv"), sep = "\t", quote = F, row.names = F)
# 

# get LYZ gene expression 
gene_expr_lyz = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/CL_0002057_IFNg_24h.tsv", sep = "\t", stringsAsFactors = F, header = T)
lyz_gene_expr = t(gene_expr_lyz[cis_gene, , drop = F])

exprs_lyz = read.table(file = file.path(resdir, "LYZ_ifng_eigens_integrated.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs_lyz) = exprs_lyz$sample_id

# Run the analysis for LYZ
print("Start mediation")
#res1 = mediation_wrapper(gen_file = geno, exprs_data = exprs_lyz, snps = trans_lead, cis_gene_expr = lyz_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")

exprs_lyz = read.table(file = file.path(resdir, "LYZ_ifng_eigens_separate.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs_lyz) = exprs_lyz$sample_id

# Run the analysis for LYZ
print("Start mediation")
#res2 = mediation_wrapper(gen_file = geno, exprs_data = exprs_lyz, snps = trans_lead, cis_gene_expr = lyz_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")

#write.table(res1, file = file.path(resdir, paste0(cis_gene, "_integrated_mediation_ifng.tsv")), sep = "\t", quote = F, row.names = F)
#write.table(res2, file = file.path(resdir, paste0(cis_gene, "_separate_mediation_ifng.tsv")), sep = "\t", quote = F, row.names = F)

# SLC39A8 in LPS 24h
cis_gene = "ENSG00000138821" # SLC
trans_lead = "chr4_102325419_ACACT_A"
cl_ids = names(all_clusters_filt[[trans_lead]])
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
expr_meta = resgenome %>% dplyr::filter(cl_id2 %in% cl_ids, meta_id == "chr4_102325419_ACACT_A", qtl_group == "CL_0002057_LPS_24h") %>% dplyr::select(meta_id, cl_id2, cl_id, cl_method, approach, qtl_group) %>% distinct()

# extract cluster expressions
# exprs_slc = data.frame()
# for(i in 1:nrow(expr_meta)){
#    approach = expr_meta[i,][["approach"]]
#    cl_method = expr_meta[i,][["cl_method"]]
#    cl_id = expr_meta[i,][["cl_id"]]
#    qtl_group = expr_meta[i,][["qtl_group"]] # we test only in this group
#    eigenfile = paste(cl_method, qtl_group, ifelse(approach=="separate", "separate", "full"), "eigen.tsv", sep = "_")
#    eigenpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results",
#                          ifelse(approach == "separate", "separate_coexpr", "integrated_coexpr"),
#                          cl_method, "eigens", eigenfile)
#    cl_eigen = read.table(eigenpath, sep = "\t", header = T, stringsAsFactors = F)
#    cl_eigen = cl_eigen[,c(cl_id, "genotype_id")]
#    cl_eigen$sample_id = row.names(cl_eigen)
#    if(nrow(exprs_slc) == 0){
#      print(cl_id)
#      exprs_slc = cl_eigen
#    }else{
#      exprs_slc = merge(exprs_slc, cl_eigen, by = c("sample_id", "genotype_id"))
#    }
# }
# 
# write.table(exprs_slc, file = file.path(resdir, "SLC39A8_eigens.tsv"), sep = "\t", quote = F, row.names = F)

# get SLC39A8 gene expression in naive
gene_expr_slc = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/CL_0002057_LPS_24h.tsv", sep = "\t", stringsAsFactors = F, header = T)
slc_gene_expr = t(gene_expr_slc[cis_gene, , drop = F])

exprs_slc = read.table(file = file.path(resdir, "SLC39A8_eigens.tsv"), sep = "\t", stringsAsFactors = F, header = T)
row.names(exprs_slc) = exprs_slc$sample_id

# Run the analysis for SLC39A8
print("Start mediation")
res = mediation_wrapper(gen_file = geno, exprs_data = exprs_slc, snps = trans_lead, cis_gene_expr = slc_gene_expr, cov_data = T, resdir = resdir)
print("Everything is ready!")
write.table(res, file = file.path(resdir, paste0(cis_gene, "_mediation_LPS24h.tsv")), sep = "\t", quote = F, row.names = F)

### Analyse the results

# SLC
med_res1 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000138821_mediation_LPS24h.tsv", stringsAsFactors = F, header = T)
estimates = melt(med_res1[,c("snp", "cis", "trans", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                 id.vars = c("snp", "cis", "trans"))
estimates$ci = unlist(med_res1[,match(paste0(gsub("_estimate", "", as.character(estimates$variable)), "_95CI"), colnames(med_res1))])
estimates = estimates %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p1 = ggplot(estimates) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + geom_hline(yintercept = 0, linetype="dashed", 
                                                                             color = "red", size=0.5) + 
  theme_bw() + coord_flip() + xlab("") + facet_grid(~trans) + 
  ggtitle("Mediation effects for SLC39A8 locus\n(monocytes, LPS 24h)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 14))

ggsave(p1, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/SLC_mediation.jpeg", dpi = 600, width = 7)

# ARHGEF3
med_res2 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000163947_mediation.tsv", stringsAsFactors = F, header = T)
estimates2 = melt(med_res2[,c("snp", "cis", "trans", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                 id.vars = c("snp", "cis", "trans"))
estimates2$ci = unlist(med_res2[,match(paste0(gsub("_estimate", "", unique(as.character(estimates2$variable))), "_95CI"), colnames(med_res2))])
estimates2 = estimates2 %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p2 = ggplot(estimates2) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + geom_hline(yintercept = 0, linetype="dashed", 
                                                                             color = "red", size=0.5) + 
 theme_bw() + coord_flip() + xlab("") + facet_grid(~trans) + 
  ggtitle("Mediation effects for ARHGEF3 locus\n(platelets)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 12))
ggsave(p2, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ARHGEF3_mediation.jpeg", dpi = 600, width = 12)

# LYZ, IFNg 24h integrated and separate
med_res3 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_integrated_mediation_ifng.tsv", stringsAsFactors = F, header = T)
med_res3$approach = "integrated"
med_res4 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_separate_mediation_ifng.tsv", stringsAsFactors = F, header = T)
med_res4$approach = "separate"
lyz_res = rbind(med_res3, med_res4)
lyz_res$condition = "IFNg_24h"
estimates3 = melt(med_res3[,c("snp", "cis", "trans", "approach", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                  id.vars = c("snp", "cis", "trans", "approach"))
estimates3$ci = unlist(med_res3[,match(paste0(gsub("_estimate", "", unique(as.character(estimates3$variable))), "_95CI"), colnames(med_res3))])
estimates3 = estimates3 %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p3 <- ggplot(estimates3) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = 0, linetype="dashed", 
            color = "red", size=0.5) + 
  theme_bw() + coord_flip() + xlab("") + facet_wrap(trans~., ncol = 3, scales = "free_y") + 
  ggtitle("Mediation effects for LYZ locus\n(integrated approach, monocytes IFNg 24h)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 10))

ggsave(p3, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/LYZ_ifng_integrated_mediation.png", dpi = 600, height = 12, width = 14)

estimates4 = melt(med_res4[,c("snp", "cis", "trans", "approach", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                  id.vars = c("snp", "cis", "trans", "approach"))
estimates4$ci = unlist(med_res4[,match(paste0(gsub("_estimate", "", unique(as.character(estimates4$variable))), "_95CI"), colnames(med_res4))])
estimates4 = estimates4 %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p4 = ggplot(estimates4) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "red", size=0.5) + 
  theme_bw() + coord_flip() + xlab("") + facet_wrap(trans~., ncol = 3, scales = "free_y") + 
  ggtitle("Mediation effects for LYZ locus\n(separate approach, monocytes IFNg 24h)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 10))
ggsave(p4, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/LYZ_ifng_separate_mediation.png", dpi = 600, height = 12, width = 14)


# LYZ, naive integrated and separate
med_res3 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_integrated_mediation.tsv", stringsAsFactors = F, header = T)
med_res3$approach = "integrated"
med_res4 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_separate_mediation.tsv", stringsAsFactors = F, header = T)
med_res4$approach = "separate"
lyz_res_naive = rbind(med_res3, med_res4)
lyz_res_naive$condition = "naive"

estimates3 = melt(med_res3[,c("snp", "cis", "trans", "approach", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                  id.vars = c("snp", "cis", "trans", "approach"))
estimates3$ci = unlist(med_res3[,match(paste0(gsub("_estimate", "", unique(as.character(estimates3$variable))), "_95CI"), colnames(med_res3))])
estimates3 = estimates3 %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p3 <- ggplot(estimates3) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "red", size=0.5) + 
  theme_bw() + coord_flip() + xlab("") + facet_wrap(trans~., ncol = 3, scales = "free_y") + 
  ggtitle("Mediation effects for LYZ locus\n(integrated approach, monocytes naive)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 10))

ggsave(p3, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/LYZ_naive_integrated_mediation.png", dpi = 600, height = 12, width = 14)

estimates4 = melt(med_res4[,c("snp", "cis", "trans", "approach", "ACME_estimate", "ADE_estimate", "total_effect_estimate")], 
                  id.vars = c("snp", "cis", "trans", "approach"))
estimates4$ci = unlist(med_res4[,match(paste0(gsub("_estimate", "", unique(as.character(estimates4$variable))), "_95CI"), colnames(med_res4))])
estimates4 = estimates4 %>% rowwise() %>% mutate(lower = as.numeric(strsplit(ci, "...", fixed = T)[[1]][1]), upper = as.numeric(strsplit(ci, "...", fixed = T)[[1]][2]))

p4 = ggplot(estimates4) + geom_point(aes(x = variable, y = value), size= 3) + 
  geom_linerange(aes(x = variable, ymin = lower, ymax = upper)) + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "red", size=0.5) + 
  theme_bw() + coord_flip() + xlab("") + facet_wrap(trans~., ncol = 3, scales = "free_y") + 
  ggtitle("Mediation effects for LYZ locus\n(separate approach, monocytes naive)") + 
  ylab("") + theme(title = element_text(size = 16), axis.text = element_text(size = 14), 
                   strip.text = element_text(size = 10))
ggsave(p4, file = "/gpfs/hpc/home/liiskolb/transqtl_final/mediation/LYZ_naive_separate_mediation.png", dpi = 600, height = 12, width = 14)

full_lyz = rbind(lyz_res, lyz_res_naive)


# merge mediation results into a single file
med_res1 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000138821_mediation_LPS24h.tsv", stringsAsFactors = F, header = T)
med_res1$approach = "separate"
med_res1$cis_gene = "SLC39A8"
med_res1$condition = "monocytes_LPS_24h"

med_res2 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000163947_mediation.tsv", stringsAsFactors = F, header = T)
med_res2$approach = ifelse(med_res2$trans == "IC68", "integrated", "separate")
med_res2$cis_gene = "ARHGEF3"
med_res2$condition = "platelet"

med_res3 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_integrated_mediation.tsv", stringsAsFactors = F, header = T)
med_res3$approach = "integrated"
med_res3$cis_gene = "LYZ"
med_res3$condition = "monocytes_naive"

med_res4 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_separate_mediation.tsv", stringsAsFactors = F, header = T)
med_res4$approach = "separate"
med_res4$cis_gene = "LYZ"
med_res4$condition = "monocytes_naive"

med_res5 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_integrated_mediation_ifng.tsv", stringsAsFactors = F, header = T)
med_res5$approach = "integrated"
med_res5$cis_gene = "LYZ"
med_res5$condition = "monocytes_IFNg_24h"

med_res6 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/mediation/ENSG00000090382_separate_mediation_ifng.tsv", stringsAsFactors = F, header = T)
med_res6$approach = "separate"
med_res6$cis_gene = "LYZ"
med_res6$condition = "monocytes_IFNg_24h"

all_med_res = do.call(rbind, list(med_res1, med_res2, med_res3, med_res4, med_res5, med_res6))
all_med_res$ACME_p_val[all_med_res$ACME_p_val == 0] = "<2e-16"
all_med_res$ADE_p_val[all_med_res$ADE_p_val == 0] = "<2e-16"
all_med_res$total_effect_p_val[all_med_res$total_effect_p_val == 0] = "<2e-16"
all_med_res$prop_mediated_p_val[all_med_res$prop_mediated_p_val == 0] = "<2e-16"

write.table(all_med_res, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/Supplementary_file4.tsv", sep = "\t", row.names = F, quote = F)
