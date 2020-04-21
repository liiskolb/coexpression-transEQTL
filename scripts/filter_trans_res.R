#!/usr/bin/env Rscript

# Filter trans results

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-e", "--eqtl", help="trans-eQTL result file (bgz)")
parser$add_argument("-r", "--resdir", help="Directory of results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if (!is.null(args$eqtl)) {
  if(is.null(args$resdir)){
    args$resdir = "."
  }else{
    if (!dir.exists(args$resdir)) { dir.create(args$resdir) }
  }
  print(paste("eQTL results are in:", args$eqtl))
  print(paste("The results will be saved to:", args$resdir))
}else{
  args$print_help
  stop("At least eQTL result file should be presented", call.=FALSE)
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readr))

filter_trans <- function(trans_resfile, eqtlresdir, pv_threshold = 0.05/1000000){
  trans_res <- readr::read_delim(trans_resfile, delim = "\t", col_names=F)
  colnames(trans_res) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
  trans_res_filt <- trans_res %>% dplyr::filter(pvalue <= pv_threshold)
  filtdir <- file.path(eqtlresdir, "filtered")
  if (!dir.exists(filtdir)) { dir.create(filtdir) }
  resfile <- file.path(filtdir, paste0(sub('\\..*$', '', basename(trans_resfile)), "_filtered", ".tsv"))
  print(resfile)
  readr::write_delim(trans_res_filt, path = resfile, delim="\t")
  return(trans_res_filt)
}

trans_resfile = args$e
eqtlresdir <- file.path(args$resdir, "eQTLres")
trans_resfile <- file.path(eqtlresdir, basename(trans_resfile))

trans_res_filt = filter_trans(trans_resfile, eqtlresdir, pv_threshold = 0.05/1000000)
