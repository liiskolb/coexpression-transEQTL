#!/usr/bin/env Rscript

# Add BY FDR values to filtered files


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readr))

path1 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"
path2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"

# methods = list.dirs(path1, recursive = F)
# approach = "integrated"
# for(m in methods){
#   print(m)
#   cl_method = basename(m)
#   eqtlresdir <- file.path(m, "eQTLres")
#   #dir.create(file.path(eqtlresdir, "filteredFDR"))
#   trans_resfiles <- list.files(eqtlresdir, pattern = "*gz$")
#   for (trans_resfile in trans_resfiles){
#     print(trans_resfile)
#     qtl_group = sub(x = basename(trans_resfile), pattern = "eigen_eQTLres.gz", replacement = "")
#     trans_resfile2 <- file.path(eqtlresdir, basename(trans_resfile))
#     trans_res <- readr::read_delim(trans_resfile2, delim = "\t", col_names=F)
#     colnames(trans_res) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
#     trans_res$FDRBY = p.adjust(trans_res$pvalue, method = "BY")
#     # add BY FDR values to filtered files
#     # filtered filename
#     filename = paste(sub('\\.gz$', "", trans_resfile), "filtered.tsv", sep = "_")
#     filt_trans_resfile = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), cl_method, "eQTLres", "filtered", filename)
#     
#     trans_res_filt = readr::read_delim(filt_trans_resfile, delim = "\t", col_names=T)
#     # add column with BY FDR value
#     trans_res = trans_res %>% filter(snp %in% trans_res_filt$snp, cluster %in% trans_res_filt$cluster) %>% select(snp, cluster, FDRBY)
#     trans_res_filt = merge(trans_res_filt, trans_res, by = c("snp", "cluster"), sort = FALSE, all.x = TRUE, all.y = FALSE)
#     
#     message("writing file to disk")
#     numCores <- 8
#     print(paste("number of cores:", numCores))
#     filename2 = paste(sub('\\.gz$', "", trans_resfile), "filtered2.tsv", sep = "_")
#     filt_trans_resfile2 = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), cl_method, "eQTLres", "filtered", filename2)
#     data.table::fwrite(trans_res_filt, file = filt_trans_resfile2, quote = F, sep = "\t", col.names = T, nThread = numCores, verbose = TRUE, showProgress = TRUE)
#   }
# }


methods = list.dirs(path2, recursive = F)
methods = c("/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/WGCNA", "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/PLIER")
approach = "separate"
for(m in methods){
  print(m)
  cl_method = basename(m)
  eqtlresdir <- file.path(m, "eQTLres")
  #dir.create(file.path(eqtlresdir, "filteredFDR"))
  trans_resfiles <- list.files(eqtlresdir, pattern = "*gz$")
  for (trans_resfile in trans_resfiles){
    print(trans_resfile)
    qtl_group = sub(x = basename(trans_resfile), pattern = "eigen_eQTLres.gz", replacement = "")
    trans_resfile2 <- file.path(eqtlresdir, basename(trans_resfile))
    trans_res <- readr::read_delim(trans_resfile2, delim = "\t", col_names=F)
    colnames(trans_res) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
    trans_res$FDRBY = p.adjust(trans_res$pvalue, method = "BY")
    # add BY FDR values to filtered files
    # filtered filename
    filename = paste(sub('\\.gz$', "", trans_resfile), "filtered.tsv", sep = "_")
    filt_trans_resfile = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), cl_method, "eQTLres", "filtered", filename)
    
    trans_res_filt = readr::read_delim(filt_trans_resfile, delim = "\t", col_names=T)
    # add column with BY FDR value
    trans_res = trans_res %>% filter(snp %in% trans_res_filt$snp, cluster %in% trans_res_filt$cluster) %>% select(snp, cluster, FDRBY)
    trans_res_filt = merge(trans_res_filt, trans_res, by = c("snp", "cluster"), sort = FALSE, all.x = TRUE, all.y = FALSE)
    
    message("writing file to disk")
    numCores <- 8
    print(paste("number of cores:", numCores))
    filename2 = paste(sub('\\.gz$', "", trans_resfile), "filtered2.tsv", sep = "_")
    filt_trans_resfile2 = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), cl_method, "eQTLres", "filtered", filename2)
    data.table::fwrite(trans_res_filt, file = filt_trans_resfile2, quote = F, sep = "\t", col.names = T, nThread = numCores, verbose = TRUE, showProgress = TRUE)
  }
}

