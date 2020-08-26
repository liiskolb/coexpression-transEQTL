#!/usr/bin/env Rscript

library(dplyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("systemPipeR")
library(systemPipeR)
library(proxy)
library(ggplot2)
library(igraph)
library(data.table)

### Credible sets 


gatherCredibleSets <- function(path, save = T, respath = NULL){
  method = sub(pattern = "_coexpr", replacement = "\\1", basename(path))
  dirs = list.dirs(path, recursive = F, full.names = TRUE)
  if(!is.null(respath)){
    resfile = file.path(respath, paste0("credible_sets_", method, ".tsv"))
  }
  res = data.frame()
  for (m in dirs){
    print(basename(m))
    mpath = file.path(m, "eQTLres", "filtered", "finemapped")
    
    finemappingfiles = list.files(path = mpath, pattern = "*\\.rds$", full.names = TRUE)
    for (f in finemappingfiles){
      cell_type = basename(f)
      cell_type = sub(pattern = "(.*)\\..*$", replacement = "\\1", cell_type)
      if(method == "separate"){
        cell_type = sub(pattern = "_separate_eigen_eQTLres_filtered_finemapped", replacement = "\\1", cell_type)
      }else{
        cell_type = sub(pattern = "_full_eigen_eQTLres_filtered_finemapped", replacement = "\\1", cell_type)
      }
      cell_type = sub(pattern = paste0(basename(m), "_"), replacement = "\\1", cell_type)
      print(cell_type)
      fres = readRDS(f) # list of lists with susie results per cluster
      subres = data.frame()
      for (cl in names(fres)){
        print(cl)
        for (trans_lead in names(fres[[cl]])){
          print(trans_lead)
          crediblesets = fres[[cl]][[trans_lead]]$sets
          if (is.null(crediblesets$cs)){
            next
          }
          for (L in seq(1, length(crediblesets$cs))){
            i = crediblesets$cs[[L]] 
            pip = fres[[cl]][[trans_lead]]$pip[i]
            zscore = fres[[cl]][[trans_lead]]$z[i]
            snps = names(fres[[cl]][[trans_lead]]$X_column_scale_factors)[i]
            z = data.frame(cbind(i, pip, zscore, snps))
            colnames(z) = c("position", "PIP", "Z-score", "snp_id")
            z$cs = L
            z$trans_lead = trans_lead
            z$cl_id = cl
            z$cl_method = basename(m)
            z$approach = method 
            z$qtl_group = cell_type
            res = rbind(res, z)
          }
        }
      }
    }
    
  }
  if(save){
    write.table(res, file = resfile, quote = F, row.names = F, sep = "\t")
  }
  
  print("done")
  return(res)
}


path1 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"
path2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"

respath = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets"

res1 = gatherCredibleSets(path1, save = F, respath = respath) # separate
res2 = gatherCredibleSets(path2, save = F, respath = respath) # integrated

# merge all the results together
res = rbind(res1, res2)
res$cs_id = paste(res$cl_id, res$cl_method, res$approach, res$qtl_group, res$trans_lead, "cs", res$cs, sep = "_")
res$snp_id = as.character(res$snp_id)
res$chr = unlist(lapply(strsplit(res$snp_id, "_"), function(x) x[[1]][1]))
res$PIP = as.numeric(as.character(res$PIP))
res$cl_id2 = ifelse(res$approach == "separate", paste(res$cl_id, res$cl_method, res$approach, res$qtl_group, sep = "_"), paste(res$cl_id, res$cl_method, res$approach, sep = "_"))


#### run script to calculate Benjamini & Yekutieli FDR for eQTL leads (trans_FDR.R) ####


#### add p-values to trans_lead to have additional Bonferroni threshold ####

res$lead_pval = 1
res$lead_FDR = 1
res$lead_BYFDR = 1

cl_summary1 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/integrated_cl_summary.tsv", sep = "\t", header = T, stringsAsFactors = F)
cl_summary2 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/separate_cl_summary.tsv", sep = "\t", header = T, stringsAsFactors = F)

subres = res %>% select(trans_lead, cl_id, cl_method, approach, qtl_group) %>% distinct()
for (i in 1:nrow(subres)){
  approach = subres[i,][["approach"]]
  cl_method = subres[i,][["cl_method"]]
  qtl_group = subres[i,][["qtl_group"]]
  trans_lead = subres[i,][["trans_lead"]]
  #filename = ifelse(approach == "separate", paste(cl_method, qtl_group, "separate_eigen_eQTLres_filtered.tsv", sep = "_"), paste(cl_method, qtl_group, "full_eigen_eQTLres_filtered.tsv", sep = "_"))
  filename = ifelse(approach == "separate", paste(cl_method, qtl_group, "separate_eigen_eQTLres_filtered2.tsv", sep = "_"), paste(cl_method, qtl_group, "full_eigen_eQTLres_filtered2.tsv", sep = "_"))
  trans_resfile = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), cl_method, "eQTLres", "filtered", filename)
  trans_res_filt = readr::read_delim(trans_resfile, delim = "\t", col_names=T)
  sub_trans_res_filt = trans_res_filt %>% filter(snp == trans_lead)
  res[res$trans_lead==trans_lead & res$approach == approach & res$cl_method == cl_method & res$qtl_group == qtl_group,"lead_pval"] = sub_trans_res_filt$pvalue[1]
  res[res$trans_lead==trans_lead & res$approach == approach & res$cl_method == cl_method & res$qtl_group == qtl_group,"lead_FDR"] = sub_trans_res_filt$FDR[1]
  res[res$trans_lead==trans_lead & res$approach == approach & res$cl_method == cl_method & res$qtl_group == qtl_group,"lead_BYFDR"] = sub_trans_res_filt$FDRBY[1]
  
  # add pvalue threshold
  if(approach == "integrated"){
    nr_clusters = cl_summary1[cl_summary1$cl_method == cl_method,"nr_of_cl"]
  }else{
    nr_clusters = cl_summary2[cl_summary2$cl_method == cl_method & cl_summary2$qtl_group == qtl_group,"nr_of_cl"]
  }
  res[res$trans_lead==trans_lead & res$approach == approach & res$cl_method == cl_method & res$qtl_group == qtl_group,"pval_thr"] = 0.05/(1000000*nr_clusters)
}

resgenome = res
resfdr = resgenome %>% filter(lead_FDR < 0.1) # with FDR 10%
resbonf = resgenome %>% filter(lead_pval < pval_thr) # with bonferroni threshold
resby = resgenome %>% filter(lead_BYFDR < 0.1) # BY FDR 10%

#### Finding meta crediblesets ####

flattenSimMatrix <- function(cormat) {
  ut <- upper.tri(cormat, diag = T)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    sim  = (cormat)[ut]
  )
}

get_metacrediblesets <- function(res, thr = 0){
  # Find all possible pairwise overlaps between credible sets
  res = data.table(res)
  intersectmat = dcast(res[,c("snp_id","cs_id")], snp_id~cs_id, length) # snps in rows and credible sets in columns, 0-1 values 
  row.names(intersectmat) = intersectmat$snp_id
  intersectmat$snp_id = NULL
  
  ## Calculate Jaccard similarity matrix between crediblesets
  
  jaccardmat = proxy::simil(intersectmat, by_rows = F, method = "Jaccard")
  jaccardmat = as.matrix(jaccardmat)
  diag(jaccardmat) = 1 # similarity 1 with itself
  
  jaccardat = flattenSimMatrix(jaccardmat)
  
  jaccardat$row = as.character(jaccardat$row)
  jaccardat$column = as.character(jaccardat$column)
  
  graphdata <- jaccardat[jaccardat$sim > thr, ]
  g <- graph_from_data_frame(graphdata[,c("row", "column")], directed = F)
  
  ## Find connected components from this graph and define these as meta credible sets
  clu <- components(g) 
  metacs <- igraph::groups(clu)
  print(paste("Number of components:", clu$no))
  
  ## Find the combined snp lists for each metacs
  metacs_snps <- lapply(metacs, function(x) unique(res$snp_id[res$cs_id %in% x]))
  
  metacs_intersection <- lapply(metacs, function(x) Reduce(intersect, lapply(x, function(y) res$snp_id[res$cs_id %in% y])))
  
  metacs_snps_full <- lapply(metacs, function(x) names(table(res$snp_id[res$cs_id %in% x]))[which.max(table(res$snp_id[res$cs_id %in% x]))])

  ## Name meta credible sets by the intersecting snp with biggest average PIP across credible sets in the metagroup
  ## In case of metasets with no intersection, choose the snp that is present in the majority of cs-s
  ## Note that in some cases the average PIPs are equal for all the snps
  
  metacs_leads <- c()
  for (i in seq(1, length(metacs_intersection))){
    if (identical(metacs_intersection[[i]], character(0))){
      metacs_intersection[[i]] = metacs_snps_full[[i]]
    }
    if (length(metacs_intersection[[i]]) == 1){
      metacs_leads <- append(metacs_leads, metacs_intersection[[i]])
    } else{
      subres = res[res$snp_id %in% metacs_intersection[[i]] & res$cs_id %in% metacs[[i]],]
      avepips = subres %>% group_by(snp_id) %>% summarize(ave_pip = mean(PIP)) %>% arrange(-ave_pip)
      metacs_leads <- append(metacs_leads, avepips$snp_id[1])
    }
  }
  names(metacs_snps) <- metacs_leads
  names(metacs) <- metacs_leads
  return(list("meta_snps" = metacs_snps, "meta_cs" = metacs))
}


metafdr = get_metacrediblesets(resfdr, thr = 0) # FDR < 0.1
metabonf = get_metacrediblesets(resbonf, thr = 0) # Bonferroni thr
metagenome = get_metacrediblesets(resgenome, thr = 0) # Genome-wide thr
metaby = get_metacrediblesets(resby, thr = 0) # BY FDR < 0.1

### Add meta_ids to the results table
flatmetares = unlist(metagenome[["meta_snps"]])

metadf = data.frame("meta_id" = names(flatmetares), "snp_id" = flatmetares, row.names=NULL)
metadf$meta_id = gsub("^\\d+|\\d+$", "", metadf$meta_id)

resgenome = merge(resgenome, metadf, by.x = "snp_id", by.y = "snp_id")

write.table(resgenome, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/all_crediblesets.tsv", sep = "\t", quote = F, row.names = F)
saveRDS(resgenome, file = file.path(respath, "all_crediblesets.rsd"))

## Perform trans eQTL analysis over all genes for each of the genomewide meta_ids
## results are filtered in the get_trans_for_meta.R script 
## using the resgenome dataframe from all_crediblesets.tsv file

## Read in the results after adding columns for trans-filtering
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))


## Summary table & find closest genes ####

get_metacs_summary = function(metacs){
  genes_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/genes_metadata.tsv", sep= "\t", header = T, stringsAsFactors = F)
  genes_metadata2 = unique(genes_metadata[,c("gene_id", "chromosome", "gene_start", "gene_end", "gene_name", "gene_type")])
  genes_metadata2$chromosome = paste0("chr", genes_metadata2$chromosome)
  genes_metadata2$gene_start = as.numeric(genes_metadata2$gene_start)
  genes_metadata2$gene_end = as.numeric(genes_metadata2$gene_end)
  
  findist = function(x){
    subset = genes_metadata2[genes_metadata2$chromosome == x[["chr"]],]
    pos = as.integer(x[["position"]])
    genes = subset %>% filter((pos >= gene_start) & (pos <= gene_end))
    incl_genes = paste0(genes$gene_name, collapse = " ")
    # detect other closest genes
    dists = base::apply(subset[,c("gene_start", "gene_end")], 1, function(y) min(abs(y[1]-pos), abs(y[2]-pos)))
    names(dists) = subset$gene_name
    dists = sort(dists)
    close_genes = paste0(names(dists)[1:5], collapse = " ")
    return(paste(incl_genes, close_genes, sep = "; "))
  }
  
  metacs_summary = data.frame("metacs_lead" = names(metacs[["meta_cs"]]), 
                              "nr_snps" = unlist(lapply(metacs[["meta_snps"]], length)),
                              "nr_cs(nr_clusters)" = unlist(lapply(metacs[["meta_cs"]], length)))
  metacs_summary$chr = unlist(lapply(metacs_summary$metacs_lead, function(x) strsplit(as.character(x), "_")[[1]][1]))
  metacs_summary$position = unlist(lapply(metacs_summary$metacs_lead, function(x) as.numeric(strsplit(as.character(x), "_")[[1]][2])))
  metacs_summary$closest_genes = apply(metacs_summary, 1, findist)
  row.names(metacs_summary) = NULL
  metacs_summary$chr = as.numeric(gsub("chr", "", metacs_summary$chr))
  metacs_summary = metacs_summary %>% arrange(chr)
  return(metacs_summary)
}

####### gprofiler link creation (implemented before adding it to the official package) #######
gp_globals = new.env()

gp_globals$version =
  tryCatch(
    utils::packageVersion("gprofiler2"),
    error = function(e) { return("unknown_version") }
  );

# Set SSL version to TLSv1_2 with fallback to TLSv1_1
# CURL_SSLVERSION_SSLv3 is not used due to the SSLv3 vulnerability <https://access.redhat.com/articles/1232123>
# CURL_SSLVERSION_TLSv1_3 is not widespread enough to have a built-in LibreSSL support yet.
# (curl's authors may decide to change it at some point, so links to the source are provided.)
gp_globals$CURL_SSLVERSION_TLSv1_1 <- 5L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1925>
gp_globals$CURL_SSLVERSION_TLSv1_2 <- 6L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1926>

gp_globals$rcurl_opts =
  RCurl::curlOptions(useragent = paste("gprofiler2/", gp_globals$version, sep=""), sslversion = gp_globals$CURL_SSLVERSION_TLSv1_2)

get_short_link <- function(query, organism = "hsapiens", multi_query = FALSE, filter = "default", ordered_query = FALSE, significant = TRUE, 
                           exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE,
                           user_threshold = 0.05, correction_method = "g_SCS",
                           domain_scope = "annotated", custom_bg = NULL,
                           numeric_ns = "", sources = NULL){
 
  url = file.path("https://biit.cs.ut.ee", "gplink", "l")
  
  if (is.null(query)) {
    stop("Missing query")
  } else if (is.list(query)) {
    if (is.data.frame(query)){
      stop("Query can't be a data.frame. Please use a vector or list of identifiers.")
    }
    # Multiple queries
    qnames = names(query)
    if (is.null(qnames)) {
      qnames = paste("query", seq(1, length(query)), sep = "_")
      names(query) = qnames
    }
    query = lapply(query, function(x) x[!is.na(x)])
  }
  else{
    query = query[!is.na(query)]
  }
  
  if(!is.null(names(query))){
    query2 = paste0(unlist(lapply(names(query), function(x) paste(">", x, "\n", paste0(query[[x]], collapse = " ")))), collapse = "\n")
    multi_query = TRUE
  }else{
    query2 = paste0(query, collapse = " ")
  }

  body <- jsonlite::toJSON((
    list(
      url = jsonlite::unbox("https://biit.cs.ut.ee/gprofiler/gost"),
      #data_version = "e98_eg45_p14_6de5f00",
      payload = {
        list(query = jsonlite::unbox(query2), 
             organism = jsonlite::unbox(organism), 
             user_threshold = jsonlite::unbox(user_threshold),
             all_results = jsonlite::unbox(!significant),
             ordered = jsonlite::unbox(ordered_query),
             no_evidences = jsonlite::unbox(!evcodes),
             combined = jsonlite::unbox(multi_query),
             measure_underrepresentation = jsonlite::unbox(measure_underrepresentation),
             no_iea = jsonlite::unbox(exclude_iea),
             domain_scope = jsonlite::unbox(domain_scope),
             numeric_ns = jsonlite::unbox(numeric_ns),
             significance_threshold_method = jsonlite::unbox(correction_method),
             background = custom_bg,
             filter = jsonlite::unbox(filter)
             ) 
      }
    )
  ),
  auto_unbox = FALSE,
  null = "null")
  
  # Headers
  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")
  
  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2 = RCurl::getCurlHandle() # Get info about the request
  
  # Request
  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  
  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }
  
  txt <- h1$value()
  res <- jsonlite::fromJSON(txt)
  shortlink = paste0('https://biit.cs.ut.ee/gplink/l/', res$result)
  return(shortlink)
}

###### Summaries after filtering for cis clusters from full trans analysis ####
respath = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets"
resgenome =  readRDS(file = file.path(respath, "all_crediblesets.rsd"))

res_filt = resgenome %>% filter(nr_trans_in_cl > 0) # leaves 303 out of 601
# filter by significant overlap in gene-level analysis
res_filt2 = res_filt %>% filter(fisher_pval_adj < 0.05) # leaves us with 247 meta_ids out of 303

resfdr_filt = res_filt2 %>% filter(lead_FDR < 0.1) # with FDR 10%
resbonf_filt = res_filt2 %>% filter(lead_pval < pval_thr) # with bonferroni threshold
resby_filt = res_filt2 %>% filter(lead_BYFDR < 0.1)  # with BY FDR 10%

metares_filt = get_metacrediblesets(res_filt2, thr = 0) # Genome-wide thr
setdiff(names(metares_filt[["meta_snps"]]), names(metagenome[["meta_snps"]]))
# rename this meta_id back to the initial id: chr6_32543752_G_A
names(metares_filt[["meta_snps"]])[names(metares_filt[["meta_snps"]]) == "chr6_32604283_T_C"] = "chr6_32543752_G_A"
names(metares_filt[["meta_cs"]])[names(metares_filt[["meta_cs"]]) == "chr6_32604283_T_C"] = "chr6_32543752_G_A"
# rename chr9_19752862_T_C back to chr9_19754478_G_A
names(metares_filt[["meta_snps"]])[names(metares_filt[["meta_snps"]]) == "chr9_19752862_T_C"] = "chr9_19754478_G_A"
names(metares_filt[["meta_cs"]])[names(metares_filt[["meta_cs"]]) == "chr9_19752862_T_C"] = "chr9_19754478_G_A"

#metafdr_filt = get_metacrediblesets(resfdr_filt, thr = 0) # FDR thr -> not used anymore
#setdiff(names(metafdr_filt[["meta_snps"]]), names(metagenome[["meta_snps"]]))

metabonf_filt = get_metacrediblesets(resbonf_filt, thr = 0) # Bonferroni thr
setdiff(names(metabonf_filt[["meta_snps"]]), names(metagenome[["meta_snps"]]))
# rename chr9_20840034_C_T to "chr9_20818520_A_G"
names(metabonf_filt[["meta_snps"]])[names(metabonf_filt[["meta_snps"]]) == "chr9_20840034_C_T"] = "chr9_20818520_A_G"
names(metabonf_filt[["meta_cs"]])[names(metabonf_filt[["meta_cs"]]) == "chr9_20840034_C_T"] = "chr9_20818520_A_G"

metaby_filt = get_metacrediblesets(resby_filt, thr = 0) # BY FDR 10%, leaves 38
setdiff(names(metaby_filt[["meta_snps"]]), names(metagenome[["meta_snps"]]))
# merge chr6_44674288_G_T and chr6_45381772_A_G to "chr6_44900099_C_T"
names(metaby_filt[["meta_snps"]])[names(metaby_filt[["meta_snps"]]) == "chr6_44674288_G_T"] = "chr6_44900099_C_T"
names(metaby_filt[["meta_cs"]])[names(metaby_filt[["meta_cs"]]) == "chr6_44674288_G_T"] = "chr6_44900099_C_T"

metaby_filt[["meta_snps"]][["chr6_44900099_C_T"]] = unique(c(metaby_filt[["meta_snps"]][["chr6_44674288_G_T"]], metaby_filt[["meta_snps"]][["chr6_45381772_A_G"]]))
metaby_filt[["meta_cs"]][["chr6_44900099_C_T"]] = unique(c(metaby_filt[["meta_cs"]][["chr6_44674288_G_T"]], metaby_filt[["meta_cs"]][["chr6_45381772_A_G"]]))
metaby_filt[["meta_snps"]][["chr6_44674288_G_T"]] = NULL
metaby_filt[["meta_snps"]][["chr6_45381772_A_G"]] = NULL
metaby_filt[["meta_cs"]][["chr6_44674288_G_T"]] = NULL
metaby_filt[["meta_cs"]][["chr6_45381772_A_G"]] = NULL

# rename chr9_20840034_C_T to "chr9_20818520_A_G"
names(metaby_filt[["meta_snps"]])[names(metaby_filt[["meta_snps"]]) == "chr9_20840034_C_T"] = "chr9_20818520_A_G"
names(metaby_filt[["meta_cs"]])[names(metaby_filt[["meta_cs"]]) == "chr9_20840034_C_T"] = "chr9_20818520_A_G"

# rename chr4_10301761_T_A to chr4_10135501_C_G
names(metaby_filt[["meta_snps"]])[names(metaby_filt[["meta_snps"]]) == "chr4_10301761_T_A"] = "chr4_10135501_C_G"
names(metaby_filt[["meta_cs"]])[names(metaby_filt[["meta_cs"]]) == "chr4_10301761_T_A"] = "chr4_10135501_C_G"

# rename chr3_192810885_G_C to chr3_192835271_G_A
names(metaby_filt[["meta_snps"]])[names(metaby_filt[["meta_snps"]]) == "chr3_192810885_G_C"] = "chr3_192835271_G_A"
names(metaby_filt[["meta_cs"]])[names(metaby_filt[["meta_cs"]]) == "chr3_192810885_G_C"] = "chr3_192835271_G_A"

## calculate summaries
metares_summary_filt = get_metacs_summary(metares_filt)
#metafdr_summary_filt = get_metacs_summary(metafdr_filt)
metabonf_summary_filt = get_metacs_summary(metabonf_filt)
metaby_summary_filt = get_metacs_summary(metaby_filt)

# merge clusters by meta_ids (after gene level significant overlap filter)
metalist_filt <- split(as.character(res_filt2$cl_id2), as.character(res_filt2$meta_id))

# update the gprofiler queries and other statistics based on filtering results
all_clusters = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters.rds")

all_clusters_filt = list()
for (m in unique(metares_summary_filt$metacs_lead)){
  all_clusters_filt[[as.character(m)]] = full_clusters[unique(metalist_filt[[as.character(m)]])]
}
saveRDS(all_clusters_filt, "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters_filt.rds")

metares_summary_filt$gp_link = "link"
metares_summary_filt$qtl_group = ""
metares_summary_filt$qtl_label = ""
sample_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/Kolberg_2020_duplicate.tsv", sep = "\t", header = T, stringsAsFactors = F)

# include gene level trans results
gene_level_trans = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/filtered_trans_res.tsv", sep = "\t", stringsAsFactors = F)
meta_trans_genes = gene_level_trans %>% dplyr::select(snp, qtl_group, gene_id) %>% group_by(snp, qtl_group) %>% nest()
d = split(meta_trans_genes$qtl_group, meta_trans_genes$snp)
d2 = lapply(names(d), function(x) sapply(d[[x]], function(y) meta_trans_genes[meta_trans_genes$snp == x  & meta_trans_genes$qtl_group == y,][["data"]][[1]][["gene_id"]], USE.NAMES = T))
names(d2) = names(d)

for (meta in unique(metares_summary_filt$metacs_lead)){
  print(meta)
  ## Add qtl groups with significant association
  qtl_groups = res_filt %>% filter(meta_id == meta) %>% select(qtl_group) %>% distinct() %>% .$qtl_group
  metares_summary_filt[metares_summary_filt$metacs_lead == meta, "qtl_group"] = paste(qtl_groups, collapse = ", ")
  ## human readable names
  qtl_names = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) %>% .$qtl_label
  metares_summary_filt[metares_summary_filt$metacs_lead == meta, "qtl_label"] = paste(qtl_names, collapse = ", ")
  
  # get genes of significant clusters
  query_clusters = all_clusters_filt[[meta]]
  # add gene level significant genes as a query (if exists)
  if (meta %in% names(d2)){
    if (is.null(names(d2[[meta]])) ){
      qtl_groups = colnames(d2[[meta]])
      qtl_names = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 

      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label 
        query_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][,q]
      }
    }else{
      qtl_groups = names(d2[[meta]])
      qtl_names = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 
      
      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label
        query_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][[q]]
      }
    }
  }

  gp_link = get_short_link(query_clusters, multi_query = TRUE)
  metares_summary_filt[metares_summary_filt$metacs_lead == meta, "gp_link"] = gp_link
}

# add information about trans/cis associations in clusters

metares_summary_filt2 = metares_summary_filt %>% rowwise() %>% mutate(nr_trans_in_cl = paste0(unique(res_filt2[res_filt2$meta_id == metacs_lead,c("meta_id", "cl_id2", "qtl_group", "nr_trans_in_cl")])[["nr_trans_in_cl"]], collapse = ";")) %>% data.frame()

## Add info if snp remains after BY-FDR and after strict bonferroni
metares_summary_filt2$bonferroni = ifelse(metares_summary_filt2$metacs_lead %in% metabonf_summary_filt$metacs_lead, TRUE, FALSE)
#metares_summary_filt2$fdr10 = ifelse(metares_summary_filt2$metacs_lead %in% metafdr_summary_filt$metacs_lead, TRUE, FALSE)
metares_summary_filt2$byfdr10 = ifelse(metares_summary_filt2$metacs_lead %in% metaby_summary_filt$metacs_lead, TRUE, FALSE)

# for the ones that have different meta_id, match as TRUE as well
bonfdiff = setdiff(unique(metabonf_summary_filt$metacs_lead), unique(metares_summary_filt2$metacs_lead))
#fdrdiff = setdiff(unique(metafdr_summary_filt$metacs_lead), unique(metares_summary_filt2$metacs_lead))
byfdrdiff = setdiff(unique(metaby_summary_filt$metacs_lead), unique(metares_summary_filt2$metacs_lead))

### old lead matching code ####
# different lead snp: chr9_20840034_C_T
"chr9_20840034_C_T" %in% metares_filt[["meta_snps"]][["chr9_20818520_A_G"]] # bonferroni and fdr
metares_summary_filt2[metares_summary_filt2$metacs_lead == "chr9_20818520_A_G", "bonferroni"] = TRUE 
metares_summary_filt2[metares_summary_filt2$metacs_lead == "chr9_20818520_A_G", "fdr10"] = TRUE 

# different lead snps: chr6_44674288_G_T and chr6_45381772_A_G (FDR separates these)
"chr6_44674288_G_T" %in% metares_filt[["meta_snps"]][["chr6_44900099_C_T"]]
"chr6_45381772_A_G" %in% metares_filt[["meta_snps"]][["chr6_44900099_C_T"]]
metares_summary_filt2[metares_summary_filt2$metacs_lead == "chr6_44900099_C_T", "fdr10"] = TRUE 

# different lead snp for chr9_19752862_T_C from fdr
"chr9_19752862_T_C" %in% metares_filt[["meta_snps"]][["chr9_19754478_G_A"]]
metares_summary_filt2[metares_summary_filt2$metacs_lead == "chr9_19754478_G_A", "fdr10"] = TRUE 

# different lead snp for chr6_32604283_T_C from fdr
"chr6_32604283_T_C" %in% metares_filt[["meta_snps"]][["chr6_32543752_G_A"]]
metares_summary_filt2[metares_summary_filt2$metacs_lead == "chr6_32543752_G_A", "fdr10"] = TRUE 

### add coloc results ####
coloc_res_summary = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/full_coloc_res_summary.tsv", sep = "\t", stringsAsFactors = F, header= T)
metares_summary_filt3 = merge(metares_summary_filt2, coloc_res_summary, by.x = "metacs_lead", by.y = "qtl_lead", all.x = TRUE)

## Add columns about clustering approach

approach_summary = res_filt2 %>% select(approach, meta_id) %>% distinct() %>% group_by(meta_id) %>% summarise(separate = "separate" %in% approach, integrated = "integrated" %in% approach)
metares_summary_filt4 = merge(metares_summary_filt3, approach_summary, by.x = "metacs_lead", by.y = "meta_id", all.x = TRUE)

## columns about clustering methods
cl_method_summary = res_filt2 %>% select(cl_method, meta_id) %>% distinct() %>% group_by(meta_id) %>% spread(cl_method, 1) %>% data.frame()
row.names(cl_method_summary) = cl_method_summary$meta_id
cl_method_summary$meta_id = NULL
cl_method_summary[!is.na(cl_method_summary)] = TRUE
cl_method_summary[is.na(cl_method_summary)] = FALSE
cl_method_summary$meta_id = row.names(cl_method_summary)
metares_summary_filt5 = merge(metares_summary_filt4, cl_method_summary, by.x = "metacs_lead", by.y = "meta_id", all.x = TRUE)

write.table(metares_summary_filt5, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", quote = F, row.names = F)

###### For integrated clustering, show the effect size (beta, last column of trans res) or -log10(pvalue) in a cluster ####

#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
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

sample_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/Kolberg_2020_duplicate.tsv", sep = "\t", header = T, stringsAsFactors = F)

get_effectsizes <- function(res, metacs){
  metamat <- stack(metacs[["meta_cs"]])
  metamat$ind <- as.character(metamat$ind)
  colnames(metamat) <- c("cs_id", "meta_id")
  # exclude meta_id column just in case
  res = res[,!(names(res) %in% c("meta_id"))]
  res2 <- merge(res, metamat, by.x = "cs_id", by.y = "cs_id")
  
  ### Add effect-sizes for all qtl_groups
  subres = res2 %>% select(cl_id, cl_method, approach, qtl_group, chr, cl_id2, meta_id) %>% unique()
  qtl_groups = unique(subres$qtl_group)
  meta_ids = unique(subres$meta_id)
  cl_ids = unique(subres$cl_id2)
  
  effectsizedata = data.frame()
  for (lead_snp in meta_ids){
    print(lead_snp)
    lead_pos = as.numeric(strsplit(lead_snp, "_")[[1]][2])
    lead_chr = as.numeric(sub(pattern = "chr", replacement = "\\1", x = strsplit(lead_snp, "_")[[1]][1]))
    sd = subres %>% filter(meta_id == lead_snp)
    for (i in 1:nrow(sd)){
      cl_method = sd$cl_method[i]
      approach = sd$approach[i]
      cl_id = sd$cl_id[i]
      cl_id2 = sd$cl_id2[i]
      if (approach == "integrated"){
        # if integrated clustering then get pvalues for every qtl group
        for (q in qtl_groups){
          qtlfilename = paste(cl_method, q, "full_eigen_eQTLres.gz", sep = "_")
          qtlfilename = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr", cl_method, "eQTLres", qtlfilename)
          # get beta value and p_value
          param <- GRanges(c(lead_chr), IRanges(lead_pos, width=1))
          sc1 <- data.frame(scanTabixDataFrame(qtlfilename, param, col_names = FALSE))
          colnames(sc1) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
          sc1 = sc1[sc1$snp == lead_snp & sc1$cluster == cl_id,]
          if(nrow(sc1) > 0){
            beta = sc1$beta
            logpval = -log10(sc1$pvalue)
            temp = data.frame("meta_id" = lead_snp, "qtl_group" = q, "approach" = approach, "cl_method" = cl_method,"cl_id2" = cl_id2, "cl_id" = cl_id, "beta" = beta, "logpval" = logpval)
            effectsizedata = rbind(effectsizedata, temp)
          }else{
            print("Missing")
            print(cl_id2)
            
          }
        }
      }
      else{
        # if separate, then get pvalue and beta only for the specific qtl group
        q = sd$qtl_group[i]
        qtlfilename = paste(cl_method, q, "separate_eigen_eQTLres.gz", sep = "_")
        qtlfilename = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr", cl_method, "eQTLres", qtlfilename)
        # get beta value and p_value
        param <- GRanges(c(lead_chr), IRanges(lead_pos, width=1))
        sc1 <- data.frame(scanTabixDataFrame(qtlfilename, param, col_names = FALSE))
        colnames(sc1) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
        sc1 = sc1[sc1$snp == lead_snp & sc1$cluster == cl_id,]
        if(nrow(sc1) > 0){
          beta = sc1$beta
          logpval = -log10(sc1$pvalue)
          temp = data.frame("meta_id" = lead_snp, "qtl_group" = q, "approach" = approach, "cl_method" = cl_method,"cl_id2" = cl_id2, "cl_id" = cl_id, "beta" = beta, "logpval" = logpval)
          effectsizedata = rbind(effectsizedata, temp)
        }else{
          print("Missing2")
          print(cl_id2)
          print(cl_method)
          print(lead_snp)
        }
      }
    }
  }
  
  effectsizedata = unique(effectsizedata)
  effectsizedata = merge(effectsizedata, unique(sample_metadata[,c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
  effectsizedata$qtl_group_label = paste(effectsizedata$cell_type_label, effectsizedata$condition_label, sep = "_")
  return(effectsizedata)
}

effectsfdr <- get_effectsizes(resfdr_filt, metafdr_filt)
effectsbonf <- get_effectsizes(resbonf_filt, metabonf_filt)
effectsres <- get_effectsizes(res_filt, metares_filt)

write.table(effectsres, file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", quote = F)

###### describe the meta credible sets with linegraphs ####

effectsfdr$qtl_group_label = as.factor(effectsfdr$qtl_group_label)

p1 <- ggplot(effectsfdr) + geom_line(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2)) + 
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = qtl_group_label, shape = approach), size = 3) +
  xlab("QTL group") + facet_wrap(~meta_id, ncol = 3, scales = "free_y") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_shape_manual(values=c(24, 8)) + 
  scale_y_continuous(limits = c(0,NA))

ggsave(p1, filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metacsfdr_linegraphs.pdf", height = 100, width = 10, limitsize = F)

effectsbonf$qtl_group_label = as.factor(effectsbonf$qtl_group_label)

p2 <- ggplot(effectsbonf) + geom_line(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2)) + 
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = qtl_group_label, shape = approach), size = 3) +
  xlab("QTL group") + facet_wrap(~meta_id, ncol = 3, scales = "free_y") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_shape_manual(values=c(24, 8)) + 
  scale_y_continuous(limits = c(0,NA))

ggsave(p2, filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metacsbonf_linegraphs.pdf", height = 6, width = 10, limitsize = F)

effectsres$qtl_group_label = as.factor(effectsres$qtl_group_label)

p3 <- ggplot(effectsres) + geom_line(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2)) + 
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = qtl_group_label, shape = approach), size = 3) +
  xlab("QTL group") + facet_wrap(~meta_id, ncol = 3, scales = "free_y") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_shape_manual(values=c(24, 8)) + 
  scale_y_continuous(limits = c(0,NA))

ggsave(p3, filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metacsres_linegraphs.pdf", height = 160, width = 10, limitsize = F)

### Add rs-ids to snp_ids ####
# tabix /gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/dbSNP_b151_GRCh38p7.vcf.gz
dbsnp = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/dbSNP_b151_GRCh38p7.vcf.gz"

resgenome$snp_pos = unlist(lapply(resgenome$snp_id, function(x) as.numeric(strsplit(x, "_")[[1]][2])))
resgenome$rs_id = ""
for(i in 1:nrow(resgenome)){
  snp_id = resgenome[i,][["snp_id"]]
  print(snp_id)
  ch = gsub("chr", "", resgenome[i, ][["chr"]])
  snp_pos = resgenome[i, ][["snp_pos"]]
  param <- GRanges(c(ch), IRanges(snp_pos, width = 1))
  sc1 <- data.frame(scanTabixDataFrame(dbsnp, param, col_names = FALSE, col_types = readr::cols(.default = "c")))
  if(nrow(sc1) == 0){
    next
  }
  colnames(sc1) = c("chr", "pos", "rs_id", "ref", "alt", "qual", "filter", "info")
  sc1$chr = paste("chr", sc1$chr, sep = "")
  sc1$snp_id = paste(sc1$chr, sc1$pos, sc1$ref, sc1$alt, sep = "_")
  resgenome$rs_id[i] = paste(sc1[sc1$snp_id == snp_id,][["rs_id"]], collapse = "; ")
}

saveRDS(resgenome, file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
write.table(resgenome, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/all_crediblesets.tsv", sep = "\t", quote = F, row.names = F)

metares_summary_filt5 = metares_summary_filt5 %>% rowwise() %>% mutate(all_rs_ids = paste(unique(resgenome[resgenome$meta_id == metacs_lead, ][["rs_id"]]), collapse = "; ")) %>% data.frame()
metares_summary_filt5 = metares_summary_filt5 %>% rowwise() %>% mutate(rs_id = unique(resgenome[resgenome$snp_id == metacs_lead, ][["rs_id"]])) %>% data.frame()

write.table(metares_summary_filt5, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", quote = F, row.names = F)

metares_summary_filt5 %>% filter(bonferroni)
metares_summary_filt5 %>% filter(fdr10) %>% head


## Add mapping to GRCh37
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
submitdata = as.character(metares_summary_filt5$metacs_lead)
submitdata = unlist(lapply(submitdata, function(x) paste(strsplit(x, "_")[[1]][1], ":", strsplit(x, "_")[[1]][2], "-", strsplit(x, "_")[[1]][2], sep = "")))
fileConn<-file("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metacs_bed.txt")
writeLines(submitdata, fileConn)
close(fileConn)

#mapdata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metacs_hg19_map.bed", stringsAsFactors = F)
mapdata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/hglft_genome_23f1e_6d7690.bed", stringsAsFactors = F)

mapdata$GRCh38 = as.character(metares_summary_filt5$metacs_lead)
colnames(mapdata) = c("GRCh37", "GRCh38")
mapdata$GRCh37 = unlist(lapply(mapdata$GRCh37, function(x) strsplit(x, "-")[[1]][1]))

metares_summary_filt6 = merge(metares_summary_filt5, mapdata, by.x = "metacs_lead", by.y = "GRCh38")
metares_summary_filt6 = metares_summary_filt6[,c("metacs_lead", "rs_id", "GRCh37", "chr", "position", 
                                                 "nr_snps", "all_rs_ids", "nr_cs.nr_clusters.", "closest_genes",
                                                 "gp_link", "qtl_group.x", "qtl_label", "nr_trans_in_cl", 
                                                 "bonferroni", "byfdr10", "integrated", "separate",
                                                 "ICA", "PLIER", "PEER", "WGCNA", "funcExplorer", 
                                                 "qtl_group.y", "trait", "trait_name")]
write.table(metares_summary_filt6, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", quote = F, row.names = F)

# update gprofiler links
metares_summary_filt6 = read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", header = T, stringsAsFactors = F)

for (meta in unique(metares_summary_filt6$metacs_lead)){
  print(meta)
  # get genes of significant clusters
  query_clusters = all_clusters_filt[[meta]]
  # add gene level significant genes as a query (if exists)
  if (meta %in% names(d2)){
    if (is.null(names(d2[[meta]])) ){
      qtl_groups = colnames(d2[[meta]])
      qtl_names = sample_metadata %>% dplyr::select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 
      
      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label 
        query_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][,q]
      }
    }else{
      qtl_groups = names(d2[[meta]])
      qtl_names = sample_metadata %>% dplyr::select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 
      
      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label
        query_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][[q]]
      }
    }
  }
  
  gp_link = get_short_link(query_clusters, multi_query = TRUE)
  metares_summary_filt6[metares_summary_filt6$metacs_lead == meta, "gp_link"] = gp_link
}
write.table(metares_summary_filt6, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", quote = F, row.names = F)


### Find gsnpense results to detect missense variants etc ####
snpense = gsnpense(unique(res_filt$rs_id))
snp_variants = snpense[["variants"]]
snpense_res = cbind(snpense[,c("rs_id", "ensgs", "gene_names")], snp_variants)
snpense_summary = merge(snpense_res, unique(res_filt[,c("rs_id", "meta_id")]), by.x = "rs_id", by.y = "rs_id")
snpense_summary2 = melt(snpense_summary, id.vars = c("meta_id", "rs_id", "ensgs", "gene_names"))
snpense_summary3 = snpense_summary2 %>% dplyr::filter(value>0)

snpense_summary3 %>% filter(variable %like% "missense") 
snpense_summary3 = snpense_summary3 %>% arrange(meta_id)
snpense_summary4 = snpense_summary3 %>% dplyr::rowwise() %>% dplyr::mutate(ensgs = paste(ensgs, collapse = "; "), gene_names = paste(gene_names, collapse = "; "))
write.table(snpense_summary4, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_snpense.tsv", sep = "\t", quote = F, row.names = F)

meta_summ = snpense_summary3 %>% group_by(meta_id) %>% summarise(rs_ids = paste(rs_id, collapse = ","), variables = paste(variable, collapse = ","))

##### Literature overview target overlap ####
library(tabulizer)
pdffile = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fgepi.22205&file=gepi22205-sup-0007-TableS3_MultiXcan_results_2019-02-06.pdf"
out <- extract_tables(pdffile)
out <- lapply(out, function(x) x[,seq(1, 19)])
final <- do.call(rbind, out)
finaldf <- as.data.frame(final[2:nrow(final), ])
colnames(finaldf) = final[1,]

arhgeftargets = finaldf[finaldf$predname == "ARHGEF3",]
arhgeftargetgenes = as.character(arhgeftargets$`predS2 obsgene`)
arhgeftargetgenes = unlist(lapply(arhgeftargetgenes, function(x) strsplit(x, " ")[[1]][2]))

our_arhgef = unique(unlist(all_clusters_filt[["chr3_56815721_T_C"]]))
length(intersect(arhgeftargetgenes, our_arhgef))

## eQTLGen
transresfile = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"
con <- gzcon(url(transresfile))
txt <- readLines(con)
dat <- read.table(textConnection(txt), sep ="\t", header = T)
arhgeftargets2 = dat %>% filter(SNP == "rs1354034")
length(intersect(arhgeftargets2$Gene, our_arhgef))

##### Enrichment statistics for associated modules ####

metares_summary_filt6 = read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/metares_full_summary_filt_coloc.tsv", sep = "\t", stringsAsFactors = F, header = T)
all_clusters_filt = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters_filt.rds")
resgenome = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))

# for every meta_id find the enrichments of the clusters and create a table for statistics across snps and methods
# include also gene-level analysis enrichment summary
enrichment_stats = data.frame()

metadata = resgenome %>% select(meta_id, cl_method, approach, qtl_group, cl_id2, cl_size) %>% distinct()
metadata = data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)

for (meta in unique(metares_summary_filt6$metacs_lead)){
  print(meta)
  # get genes of significant clusters
  query_clusters = all_clusters_filt[[meta]]
  
  # add gene level significant genes as a query (if exists)
  gene_level_clusters = list()
  if (meta %in% names(d2)){
    if (is.null(names(d2[[meta]])) ){
      qtl_groups = colnames(d2[[meta]])
      qtl_names = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 
      
      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label 
        gene_level_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][,q]
      }
    }else{
      qtl_groups = names(d2[[meta]])
      qtl_names = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct %>% filter(qtl_group %in% qtl_groups) %>% mutate(qtl_label = paste(cell_type_label, condition_label)) 
      
      for(q in qtl_groups){
        qname = qtl_names[qtl_names$qtl_group==q,]$qtl_label
        gene_level_clusters[[paste(meta, "gene-level trans (FDR5%) in", qname)]] = d2[[meta]][[q]]
      }
    }
  }
  
  # gp enrichment
  gpres1 = gost(query_clusters, sources = c("GO", "REAC", "KEGG"))
  
  # add info to table
  subdf = data.frame(meta_id = meta, stringsAsFactors = F)
  #enrichment_stats[enrichment_stats$meta_id == meta,, drop = F] 
  subdf$cl_id = paste(names(query_clusters), collapse = ",")
  subdf = subdf %>% 
    mutate(cl_id = strsplit(cl_id, ",")) %>% 
    unnest(cl_id)
  
  subdf = merge(subdf, metadata, by.x = c("meta_id", "cl_id"), by.y = c("meta_id", "cl_id2"), sort = FALSE, all.x = T) 
  subdf = subdf %>% mutate(is_enriched = ifelse(cl_id %in% gpres1$result$query, TRUE, FALSE))
  
  # gene-level enrichment
  if (length(gene_level_clusters)>0){
    gpres2 = gost(gene_level_clusters, sources = c("GO", "REAC", "KEGG"))
    subdf2 = data.frame(cl_id = names(gene_level_clusters), stringsAsFactors = F)
    subdf2$meta_id = meta
    subdf2$qtl_group = qtl_groups
    subdf2$cl_method = "gene-level"
    subdf2$approach = "gene-level"
    subdf2$cl_size = unlist(lapply(gene_level_clusters, function(x) length(x)))
    subdf2 = subdf2 %>% mutate(is_enriched = ifelse(cl_id %in% gpres2$result$query, TRUE, FALSE))
    subdf2 = subdf2[,names(subdf)]
  }
  else{
    subdf2 = data.frame(meta_id = meta, stringsAsFactors = F)
    subdf2$cl_id = NA
    subdf2$cl_method = "gene-level"
    subdf2$approach = "gene-level"
    subdf2$cl_size = 0
    subdf2$is_enriched = FALSE
  }
  enrichment_stats = rbind(enrichment_stats, subdf)
  enrichment_stats = rbind(enrichment_stats, subdf2)
}

write.table(enrichment_stats, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/enrichment_stats.tsv", sep = "\t", quote = F, row.names = F)

# summarise the enrichment statistics over methods, approaches, etc
enrichment_stats = read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/enrichment_stats.tsv", sep = "\t", header = T, stringsAsFactors = F)

loci_stats = enrichment_stats %>% filter(cl_method!="gene-level") %>% dplyr::select(cl_id, is_enriched) %>% distinct()
mean(loci_stats$is_enriched)

loci_stats = enrichment_stats %>% filter(cl_method!="gene-level") %>% group_by(meta_id, cl_method, approach) %>% summarise(ave_enriched_modules = mean(is_enriched))

# gene level enrichment for comparison
loci_stats_genelevel = enrichment_stats %>% filter(cl_method =="gene-level") 
mean(loci_stats_genelevel$is_enriched)


p1 = ggplot(loci_stats_genelevel) + geom_histogram(aes(x=ave_enriched_modules), bins = 10) + theme_bw() + 
  facet_grid(.~cl_method) + ylab("Frequency") + xlab("Average proportion of\nenriched trans genes per locus")
p2 = ggplot(loci_stats) + geom_histogram(aes(x=ave_enriched_modules), bins = 10) + theme_bw() + 
  facet_grid(cl_method~approach) + ylab("Frequency") + xlab("Average proportion of\nenriched modules per locus")

p2 + p1

##### Replication analysis table for IFNB1, ARHGEF3 and LYZ ####

### ARHGEF3 and eQTLGen; chr3_56815721_T_C (rs1354034)

eqtlgen = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/replication/eQTLgen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
arhgef_eqtlgen = eqtlgen %>% filter(SNP == "rs1354034")
nrow(arhgef_eqtlgen) # 840 FDR 5% significant trans genes
eqtlgen_genes = arhgef_eqtlgen$Gene

# gene modules
arhgefres = data.frame(locus = "ARHGEF3", rs_id = "rs1354034", meta_id = "chr3_56815721_T_C",
                       context = "platelet",
                       cl_id = names(all_clusters_filt[["chr3_56815721_T_C"]]), 
                       cl_size = unlist(lapply(all_clusters_filt[["chr3_56815721_T_C"]], function(x) length(x))),
                       domain_size = 18383, 
                       overlap_size = unlist(lapply(all_clusters_filt[["chr3_56815721_T_C"]], function(x) length(intersect(x, eqtlgen_genes)))),
                       repl_set_size = length(eqtlgen_genes),
                       study = "eQTLGen",
                       repl_rs_id = "rs1354034",
                       replication_context = "blood",
                       stringsAsFactors = F)
row.names(arhgefres) = NULL

### ARHGEF in Nath, 2017
nath = read.xlsx("/gpfs/hpc/home/liiskolb/transqtl_final/replication/13059_2017_1279_MOESM9_ESM.xlsx", sheetName = "Trans_eQTLs_PM", header = T, startRow = 2)
nath_ensg = gconvert(nath$Gene.Name)$target
# gene modules
arhgefres2 = data.frame(locus = "ARHGEF3", rs_id = "rs1354034", meta_id = "chr3_56815721_T_C",
                       context = "platelet",
                       cl_id = names(all_clusters_filt[["chr3_56815721_T_C"]]), 
                       cl_size = unlist(lapply(all_clusters_filt[["chr3_56815721_T_C"]], function(x) length(x))),
                       domain_size = 18383, 
                       overlap_size = unlist(lapply(all_clusters_filt[["chr3_56815721_T_C"]], function(x) length(intersect(x, nath_ensg)))),
                       repl_set_size = length(nath_ensg),
                       study = "Nath_2017",
                       repl_rs_id = "rs1354034",
                       replication_context = "blood",
                       stringsAsFactors = F)
row.names(arhgefres2) = NULL

### LYZ (54 modules, compare only naive modules -> 30 modules)
lyz_cluster_names = resgenome %>% filter(meta_id == "chr12_69344099_A_G") %>% dplyr::select(cl_id2, qtl_group, fisher_pval_adj) %>% distinct() %>% filter(fisher_pval_adj<0.05, qtl_group == "CL_0002057_naive") %>% .$cl_id2
lyz_clusters = all_clusters_filt[["chr12_69344099_A_G"]]
lyz_naive_clusters = lyz_clusters[lyz_cluster_names]

# Rakitsch_2016
# extracted manually from a figure in supplementary material
# https://docs.google.com/spreadsheets/d/1Vivru6Hea3Zzod5dh_JLB0rR7Q-EajCYWA5SAhsJyI0/edit?usp=sharing
library(gsheet)
rakitsch_genes = gsheet2tbl("https://docs.google.com/spreadsheets/d/1Vivru6Hea3Zzod5dh_JLB0rR7Q-EajCYWA5SAhsJyI0/edit?usp=sharing")
rakitsch_ensg = gconvert(rakitsch_genes$gene_name)
rakitsch_ensg = unique(rakitsch_ensg$target)

lyzres = data.frame(locus = "LYZ", rs_id = "rs10784774", 
                    meta_id = "chr12_69344099_A_G",
                    context = "monocytes naive",
                       cl_id = names(lyz_naive_clusters), 
                       cl_size = unlist(lapply(lyz_naive_clusters, function(x) length(x))),
                       domain_size = 18383, 
                       overlap_size = unlist(lapply(lyz_naive_clusters, function(x) length(intersect(x, rakitsch_ensg)))),
                       repl_set_size = length(rakitsch_ensg),
                       study = "Rakitsch_2016",
                    repl_rs_id = "rs6581889",
                    replication_context = "monocytes",
                       stringsAsFactors = F)
row.names(lyzres) = NULL

# Rotival_2011, rs11177644
# https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1002367.s011&type=supplementary

out <- read.table("https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1002367.s011&type=supplementary", sep = "\t", header = T, check.names = T, skip = 1)
out = out %>% filter(Cis.Or.trans == "trans")
rotival_genes = unique(out$gene.symbol)
rotival_ensg = gconvert(rotival_genes)
rotival_ensg = unique(rotival_ensg$target)

lyzres2 = data.frame(locus = "LYZ", rs_id = "rs10784774", meta_id = "chr12_69344099_A_G",
                     context = "monocytes naive",
                    cl_id = names(lyz_naive_clusters), 
                    cl_size = unlist(lapply(lyz_naive_clusters, function(x) length(x))),
                    domain_size = 18383, 
                    overlap_size = unlist(lapply(lyz_naive_clusters, function(x) length(intersect(x, rotival_ensg)))),
                    repl_set_size = length(rotival_ensg),
                    study = "Rotival_2011",
                    repl_rs_id = "rs11177644",
                    replication_context = "monocytes",
                    stringsAsFactors = F)
row.names(lyzres2) = NULL

## FF LYZ genes 
ff_lyz = strsplit("CREB1
                  DUSP19
                  RAB27A
                  EID2B
                  ZMAT3
                  SNHG10
                  SEMA3E
                  C3ORF34
                  CDKN2AIPNL
                  LRAP
                  RAG1AP1
                  MINK1
                  KILLIN
                  CHRNA5
                  LOC100128288
                  FAM119A
                  HSPC268
                  FDXACB1
                  LOC729090
                  QRFPR
                  GRIPAP1
                  MCF2L2
                  PLEKHB2
                  PPID
                  FAM63A
                  GSDM1
                  TNFSF15
                  ZNF69
                  SS18
                  LOC100129269
                  HS.580797
                  TMEM159
                  LMOD3
                  TRIM34
                  CHST12
                  DEM1
                  NDUFC2
                  LOC643509
                  NDUFV3
                  C1ORF210
                  YRDC
                  ROCK2
                  NOD1
                  LOC401152
                  MAFF
                  PHAX
                  TBCCD1
                  ZNF765
                  CEP27
                  ZNF674
                  PRRG4
                  PILRA
                  ZNF131
                  SSSCA1
                  C11ORF63
                  RASSF6
                  BLZF1
                  LOC100128126
                  RBM3
                  TNFRSF10B
                  ANKRD44
                  CES2
                  PPM1K
                  LOC644250
                  FKBP14
                  ZNF652
                  TPM4
                  ICA1
                  HIST1H2AG
                  PRPF8
                  DENND4C
                  SYAP1
                  GNL3L
                  FKBP1P1
                  COX16
                  GPR1
                  TPM4", "\n")[[1]]
ff_lyz = gsub(" ", "", ff_lyz)
ff_ens = gconvert(ff_lyz)
ff_ensg = unique(ff_ens$target)
length(intersect(ff_ensg, rakitsch_ensg))

### IFNB1 (12 modules) and Quach LPS 6h; rs12553564
ifnb_clusters = all_clusters_filt[["chr9_20818520_A_G"]]

library(xlsx)
tmp <- tempfile(fileext = ".xlsx")
download.file(url = "http://ars.els-cdn.com/content/image/1-s2.0-S009286741631306X-mmc3.xlsx", destfile = tmp, mode="wb")
quach_ifnb = read.xlsx(tmp, sheetIndex = 1, startRow = 16)
quach_ifnb = quach_ifnb %>% filter(SNPa == "rs12553564", Condition == "LPS") %>% dplyr::select(trans.regulated.genes..up.e, trans.regulated.genes..down.e)
quach_ifnb_genes = unique(c(as.character(quach_ifnb$trans.regulated.genes..up.e), strsplit(as.character(quach_ifnb$trans.regulated.genes..down.e), " // ")[[1]]))
quach_ifnb_ensg = gconvert(quach_ifnb_genes)
quach_ifnb_ensg = unique(quach_ifnb_ensg$target)

ifnbres = data.frame( locus = "IFNB1",
                      rs_id = "rs13296842", 
                      meta_id = "chr9_20818520_A_G",
                      context = "monocytes LPS 24h",
                       cl_id = names(all_clusters_filt[["chr9_20818520_A_G"]]), 
                       cl_size = unlist(lapply(all_clusters_filt[["chr9_20818520_A_G"]], function(x) length(x))),
                       domain_size = 18383, 
                       overlap_size = unlist(lapply(all_clusters_filt[["chr9_20818520_A_G"]], function(x) length(intersect(x, quach_ifnb_ensg)))),
                       repl_set_size = length(quach_ifnb_ensg),
                       study = "Quach_2016",
                       repl_rs_id = "rs12553564",
                      replication_context = "monocytes LPS 6h",
                       stringsAsFactors = F)
row.names(ifnbres) = NULL

## merge replication results to a single table
repl_table = do.call(rbind, list(arhgefres, arhgefres2, ifnbres, lyzres, lyzres2))

## calculate overlap significance

pvals = apply(repl_table, 1, 
              function(x) {
                names(x) <- names(repl_table)
                tbl <- matrix(c(as.numeric(x["overlap_size"]), 
                                as.numeric(x["cl_size"])-as.numeric(x["overlap_size"]), 
                                as.numeric(x["repl_set_size"])-as.numeric(x["overlap_size"]), 
                                as.numeric(x["domain_size"])-as.numeric(x["cl_size"])-(as.numeric(x["repl_set_size"])-as.numeric(x["overlap_size"]))), nrow = 2)
                fisher.test(tbl, alternative="greater")$p.value
              })
repl_table$fisher_pval = pvals

repl_table = repl_table %>% group_by(meta_id) %>% mutate(adj_fisher_pval = p.adjust(fisher_pval, method = "bonferroni")) %>% data.frame()
write.table(repl_table, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/replication_overlap.tsv", sep = "\t", row.names = F, quote = F)

### Independent gene-level replication comparison ####

## using Quach, BLUEPRINT and ImmVar results as "ground truth"
# Quach_2016 -> monocytes naive, monocytes LPS 6h
# Immvar -> monocytes naive, T cell CD4
# BLUEPRINT -> monocytes, neutrophils, T-cell


respath = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets"
resgenome = readRDS(file = file.path(respath, "all_crediblesets.rsd"))

#genes_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/genes_metadata.tsv", header = T, sep = "\t")
#our_genes = unique(genes_metadata$gene_id)
  
# monocytes naive
monocyte_specific = resgenome %>% dplyr::filter(qtl_group == "CL_0002057_naive", fisher_pval_adj < 0.05, lead_BYFDR<0.1) %>% dplyr::select(meta_id, qtl_group, cl_id2) %>% distinct()
fdr_thr = 0.05
# get our gene-level results
monocytes_res = data.frame()
for (i in 1:nrow(monocyte_specific)){
  meta_id = monocyte_specific[i,][["meta_id"]]
  cl_id = monocyte_specific[i,][["cl_id2"]]
  qtl_group = "CL_0002057_naive"
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
  # Quach results
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/Quach_2016/eQTLs/Quach_2016_monocyte_naive_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% filter(adj.pval < fdr_thr) # genes in replication set
    
    monocytes_res = rbind(monocytes_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id, 
                                                    "repl_study" = "Quach_2016", "qtl_group" = qtl_group, 
                                                    cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                                    #domain_size = 18383,
                                                
                                                    domain_size = nrow(df2),
                                                    overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                                    repl_set_size = length(rep_sign$gene_id)))
  }
}

# Immvar

for (i in 1:nrow(monocyte_specific)){
  meta_id = monocyte_specific[i,][["meta_id"]]
  cl_id = monocyte_specific[i,][["cl_id2"]]
  qtl_group = "CL_0002057_naive"
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
  
  # Immvar results
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/Immvar/eQTLs/ImmVar_monocyte_CD14_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% dplyr::filter(adj.pval < fdr_thr) # genes in replication set
    
    monocytes_res = rbind(monocytes_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id,  "repl_study" = "ImmVar", "qtl_group" = qtl_group, 
                                                    cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                                    #domain_size = 18383,
                                                    domain_size = nrow(df2),
                                                    overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                                    repl_set_size = length(rep_sign$gene_id)))
  }
}

# BLUEPRINT
for (i in 1:nrow(monocyte_specific)){
  meta_id = monocyte_specific[i,][["meta_id"]]
  cl_id = monocyte_specific[i,][["cl_id2"]]
  qtl_group = "CL_0002057_naive"
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
  
  # Blueprint results
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/BLUEPRINT/eQTLs/BLUEPRINT_SE.gene_counts_cqn_norm.monocyte_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% filter(adj.pval < fdr_thr) # genes in replication set
    
    monocytes_res = rbind(monocytes_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id,  "repl_study" = "BLUEPRINT", "qtl_group" = qtl_group, 
                                                    cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                                    #domain_size = 18383,
                                                    domain_size = nrow(df2),
                                                    overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                                    repl_set_size = length(rep_sign$gene_id)))
  }
}

# T-cells
tcell_specific = resgenome %>% filter(qtl_group %in% c("CL_0000625_naive", "CL_0000624_naive"), fisher_pval_adj < 0.05, lead_BYFDR<0.1) %>% dplyr::select(meta_id, qtl_group, cl_id2) %>% distinct()
tcell_res = data.frame()
# Immvar
for (i in 1:nrow(tcell_specific)){
  meta_id = tcell_specific[i,][["meta_id"]]
  cl_id = tcell_specific[i,][["cl_id2"]]
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  qtl_group = tcell_specific[i,][["qtl_group"]]
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
  
  # Immvar
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/Immvar/eQTLs/ImmVar_T-cell_CD4_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% filter(adj.pval < fdr_thr) # genes in replication set
    
    tcell_res = rbind(tcell_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id,  "repl_study" = "ImmVar", "qtl_group" = qtl_group, 
                                                    cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                                    #domain_size = 18383,
                                                    domain_size = nrow(df2),
                                                    overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                                    repl_set_size = length(rep_sign$gene_id)))
  }
}

# Blueprint
for (i in 1:nrow(tcell_specific)){
  meta_id = tcell_specific[i,][["meta_id"]]
  cl_id = tcell_specific[i,][["cl_id2"]]
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  qtl_group = tcell_specific[i,][["qtl_group"]]
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
 
  # Blueprint
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/BLUEPRINT/eQTLs/BLUEPRINT_PE.gene_counts_cqn_norm.T-cell_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% filter(adj.pval < fdr_thr) # genes in replication set
    
    tcell_res = rbind(tcell_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id,  "repl_study" = "BLUEPRINT", "qtl_group" = qtl_group, 
                                            cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                            #domain_size = 18383,
                                            domain_size = nrow(df2),
                                            overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                            repl_set_size = length(rep_sign$gene_id)))
  }
}


# neutrophils
neutrophil_specific = resgenome %>% filter(qtl_group %in% c("CL_0000775_naive"), fisher_pval_adj < 0.05, lead_BYFDR<0.1) %>% dplyr::select(meta_id, qtl_group, cl_id2) %>% distinct()
neutrophil_res = data.frame()
# Blueprint
for (i in 1:nrow(neutrophil_specific)){
  meta_id = neutrophil_specific[i,][["meta_id"]]
  cl_id = neutrophil_specific[i,][["cl_id2"]]
  # read effectsizes for given snp
  meta_chr = as.integer(gsub("chr" ,"", strsplit(meta_id, "_")[[1]][1]))
  meta_pos = as.integer(strsplit(meta_id, "_")[[1]][2])
  qtl_group = neutrophil_specific[i,][["qtl_group"]]
  param <- GRanges(c(meta_chr), IRanges(meta_pos, width=1))
 
  # Blueprint
  df2 = try(scanTabixDataFrame("/gpfs/hpc/home/liiskolb/transqtl_final/replication/BLUEPRINT/eQTLs/BLUEPRINT_SE.gene_counts_cqn_norm.neutrophil_eQTLres.gz", param, col_names = F)[[1]], TRUE)
  if(inherits(df2, "try-error")){
    next
  }
  if(is.null(df2)){
    next
  }else{
    colnames(df2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
    df2 = df2 %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
    rep_sign = df2 %>% filter(adj.pval < fdr_thr) # genes in replication set
    
   neutrophil_res = rbind(neutrophil_res, data.frame("meta_id" = meta_id, "cl_id" = cl_id,  "repl_study" = "BLUEPRINT", "qtl_group" = qtl_group, 
                                            cl_size = length(all_clusters_filt[[meta_id]][[cl_id]]),
                                            #domain_size = 18383,
                                            domain_size = nrow(df2),
                                            overlap_size = length(intersect(all_clusters_filt[[meta_id]][[cl_id]], rep_sign$gene_id)),
                                            repl_set_size = length(rep_sign$gene_id)))
  }
}

replication_res = rbind(monocytes_res, tcell_res, neutrophil_res)

repl_pvals = apply(replication_res, 1, 
                   function(x) {
                     names(x) <- names(replication_res)
                     #print(x)
                     if(as.numeric(x["overlap_size"])>0 & as.numeric(x["repl_set_size"]) > 0){
                       tbl <- matrix(c(as.numeric(x["overlap_size"]), 
                                       as.numeric(x["cl_size"])-as.numeric(x["overlap_size"]), 
                                       as.numeric(x["repl_set_size"])-as.numeric(x["overlap_size"]), 
                                       as.numeric(x["domain_size"])-as.numeric(x["cl_size"])-(as.numeric(x["repl_set_size"])-as.numeric(x["overlap_size"]))), nrow = 2)
                       
                       fisher.test(tbl, alternative="greater")$p.value
                     }else{
                       1
                     }
                     
                   })

replication_res$fisher_pval = repl_pvals

replication_res %>% dplyr::filter(fisher_pval<0.05/nrow(replication_res))
