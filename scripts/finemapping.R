#!/usr/bin/env Rscript

# Fine-mapping the trans-QTL results
# The goal is to find the Credible Sets of variables (QTLs) for every cluster using the tool SuSiE. 

# devtools::install_github("stephenslab/susieR@0.8.0",build_vignettes = FALSE)

library("GDSArray")
library("dplyr")
library("susieR")
#library("purrr")
library("readr")
library("reshape2")
library("Matrix")
library("SNPRelate")

# Code from https://github.com/kauralasoo/RNAseq_pipeline/blob/master/analysis/finemapping/susie_finemapping.R

importVariantInformationFromGDS <- function(gdsfile){
  
  # Import individual columns
  snp_pos = GDSArray::GDSArray(gdsfile, "snp.position")
  snp_chromosome = GDSArray::GDSArray(gdsfile, "snp.chromosome")
  snp_id = GDSArray::GDSArray(gdsfile, "snp.rs.id")
  
  # Make a data frame
  snp_df = dplyr::data_frame(gds_snp_id = as.integer(names(snp_id)),
                             chromosome = as.vector(snp_chromosome),
                             pos = as.vector(snp_pos),
                             snp_id = as.vector(snp_id))
  return(snp_df)
}

extractGenotypeMatrixFromGDS <- function(chr, start, end, variant_information, gdsfile){

  if(!is.null(start)){
    #Extract variant ids from variant infromation
    var_meta = dplyr::filter(variant_information, chromosome == chr, pos > start, pos < end)
  }else{
    var_meta = variant_information
  }
  gds_ids = var_meta$gds_snp_id
  
  
  #Extract genotype from the gds file
  geno = GDSArray::GDSArray(gdsfile, "genotype")
  genotype = as.matrix(geno[gds_ids,])
  rownames(genotype) = var_meta$snp_id
  
  return(genotype)
}

# Quantile normalization
tRank <- function(x) {
  t1 <- x
  t2 <- x[!is.na(x)]
  r <- qnorm(rank(t2)/(length(rank(t2))+1))
  t1[!is.na(t1)] <- r
  t1
}

# Code from https://rdrr.io/cran/varbvs/src/R/misc.R

# Adjust variables X and continuous outcome Y so that the linear
# effects of the covariates Z are removed. This is equivalent to
# integrating out the regression coefficients corresponding to the
# covariates with respect to an improper, uniform prior; see Chipman,
# George and McCulloch, "The Practical Implementation of Bayesian
# Model Selection," 2001. It is assumed that the first column of Z is
# the intercept; that is, a column of ones.
remove.covariate.effects <- function (X, Z, y) {
  # Here I compute two quantities that are used here to remove linear
  # effects of the covariates (Z) on X and y, and later on to
  # efficiently compute estimates of the regression coefficients for
  # the covariates.
  # Add intercept.
  n <- nrow(X)
  
  if (is.null(Z)){Z <- matrix(1,n,1)}
  else{Z <- cbind(1,Z)}
  
  A   <- forceSymmetric(crossprod(Z))
  #SZy <- as.vector(solve(A,c(y %*% Z)))
  SZy <- as.matrix(solve(A,t(Z) %*% y))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  if (ncol(Z) == 1) {
    X <- scale(X,center = TRUE,scale = FALSE)
    #y <- y - mean(y)
    y <- scale(y,center = TRUE,scale = FALSE)
  } else {
    
    # The equivalent expressions in MATLAB are  
    #
    #   y = y - Z*((Z'*Z)\(Z'*y))
    #   X = X - Z*((Z'*Z)\(Z'*X))  
    #
    # This should give the same result as centering the columns of X
    # and subtracting the mean from y when we have only one
    # covariate, the intercept.
    #y <- y - c(Z %*% SZy)
    y <- y - Z %*% SZy
    X <- X - Z %*% SZX
  }
  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
}


runFineMappingQTL <- function(path, gdsfile){
  variant_info = importVariantInformationFromGDS(gdsfile)
  print("variant info is in")
  
  genotype_matrix = extractGenotypeMatrixFromGDS(gdsfile = gdsfile, start = NULL, variant_information = variant_info)
  
  dirs = list.dirs(path, recursive = F, full.names = TRUE)
  
  for (m in dirs){
    print(m)
    # to run only for funcExplorer
    #if (!grepl(pattern = "funcExplorer", x = m)){
    #  next 
    #  print(paste("Not running for ", m))
    #}
    #
    mpath = file.path(m, "eQTLres", "filtered")
    trans_resfiles = list.files(path = mpath, pattern = "*\\.tsv$", full.names = TRUE)
    
    for (f in trans_resfiles){
      cell_type = basename(f)
      cell_type = sub(pattern = "(.*)\\..*$", replacement = "\\1", cell_type)
      print(cell_type)
      
      # lead snps
      trans_filtered <- readr::read_delim(f, delim = "\t")
      
      if (nrow(trans_filtered) == 0){
        print("No filtered eQTL results")
        print(cell_type)
        next
      }
      # eigenvectors
      exprs_file = file.path(m, "eigens", paste0(sub(pattern = "_eQTLres_filtered", replacement = "\\1", cell_type), ".tsv"))
      exprs_data = read.table(exprs_file, sep="\t", header = TRUE) 
      
      # remember sample_id
      keep_ids = row.names(exprs_data)
      
      # use genotype id as row id and clusters are in columns
      rownames(exprs_data) = exprs_data$genotype_id
      
      sample_ids = intersect(rownames(exprs_data), colnames(genotype_matrix))
      
      if (length(sample_ids)==0){
        message("Sample ids don't match at all")
        next
      }
      exprs_data = exprs_data[sample_ids,]
      
      exprs_data = exprs_data[,!colnames(exprs_data) %in% c("genotype_id", "sample_id")]
      phenotype_mat = as.matrix(exprs_data)
      
      # Quantile normalize eigenvectors
      print("Quantile normalization")
      phenotype_mat = apply(phenotype_mat, 2, tRank) # genotype IDs in rows
      print("Quantile done")
      
      print("Standardize cluster profiles")
      phenotype_mat = scale(phenotype_mat)
      print("Standardization done")
      
      # Covariates
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
      gen_obj = snpgdsOpen(gdsfile)
      # Try different LD thresholds for sensitivity analysis
      snpset <- SNPRelate::snpgdsLDpruning(gen_obj, sample.id = sample_ids, ld.threshold=0.3)
      # Get all selected snp id
      snpset.id <- unlist(unname(snpset))
      print(length(snpset.id))
      pca <- SNPRelate::snpgdsPCA(gen_obj, sample.id = sample_ids, snp.id = snpset.id, num.thread=3)
      # make a data.frame of first 3 PCs
      pcs <- data.frame(sample.id = pca$sample.id,
                        EV1 = pca$eigenvect[,1],    # the first eigenvector
                        EV2 = pca$eigenvect[,2],    # the second eigenvector
                        EV3 = pca$eigenvect[,3],    # the third eigenvector
                        stringsAsFactors = FALSE)
      snpgdsClose(gen_obj)
      covs = base::merge(covs, pcs, by.x = 0, by.y = "sample.id")
      rownames(covs) = covs$Row.names
      covs$Row.names = NULL
      covs = data.matrix(covs)
      covs = covs[row.names(covs) %in% sample_ids,]
      # remove covariates with any NA-s
      covs = covs[,complete.cases(t(covs))]
      print(paste("Covariates dim: ", paste0(dim(covs), collapse=",")))
      
      # Adjust for covariates just similarly as in MatrixeQTL
      genotype_mat = t(genotype_matrix[,sample_ids])
      
      adjusted = remove.covariate.effects(X = genotype_mat, Z = covs, y = phenotype_mat)
      genotype_mat_adj = adjusted[["X"]]
      phenotype_mat_adj = adjusted[["y"]]
      
      # for every cluster, take the qtl with smallest p-value and fine-map around this snp (+-500 000bp)
      # then exclude the snps from this region and take the next snp with smallest p-value and repeat
      clusters <- unique(trans_filtered$cluster)
      fullres <- list()
      
      for (cl in clusters){
        snps <- trans_filtered[trans_filtered$cluster == cl, ]
        while (dim(snps)[1] > 0){
          snps <- snps %>% arrange(pvalue)
          lead_snp <- snps[1,][["snp"]]
          lead_snp_pos <- snps[1,][["start"]]
          lead_snp_chr <- snps[1,][["chr"]]
          
          var_meta = dplyr::filter(variant_info, chromosome == lead_snp_chr, pos > lead_snp_pos - 500000, pos < lead_snp_pos + 500000)
          gds_ids = var_meta$gds_snp_id
          
          # extract genotype matrix
          geno_mat = genotype_mat_adj[,gds_ids]
          
          phenotype_values = as.vector(phenotype_mat_adj[,cl])
            
          # fit susie
          fitted <- susie(geno_mat, phenotype_values,
                          L = 10,
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = FALSE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE, 
                          compute_univariate_zscore = TRUE)
          
          fullres[[cl]][[lead_snp]] <- fitted
          
          # exclude lead snp and its region
          snps <- snps[snps$start < lead_snp_pos - 500000 | snps$start > lead_snp_pos + 500000,]
        }
      }
    
      resdir = file.path(mpath, "finemapped")
      if(!dir.exists(resdir)){
        dir.create(resdir)
      }
      saveRDS(fullres, file.path(resdir, paste0(sub('\\..*$', '', basename(f)), "_finemapped", ".rds")))
    }
  }
  print("done!")
}


# Read in sample metadata
sample_metadata = read.table(file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data", "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)
print("sample metadata in")

# Run the Fine-mapping

path1 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"
path2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"

# Genotypes GDS file
gdsfile = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.gds"

runFineMappingQTL(path1, gdsfile)
runFineMappingQTL(path2, gdsfile)
