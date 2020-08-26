# Filter trans results
library(dplyr)
library(readr)

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
transfiles = list.files(resdir, pattern = ".gz$", full.names = T)

# Get the metacs lead ids
metares = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
meta_ids = unique(metares$meta_id)

res = data.frame()
for (f in transfiles){
  trans_res <- readr::read_delim(f, delim = "\t", col_names=F)
  colnames(trans_res) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
  qtl_group = sub(x = basename(f), pattern = "_eQTLres.gz", replacement = "")
  print(qtl_group)
  trans_res = trans_res %>% group_by(snp) %>% mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
  # check if meta_ids all in results?
  if(length(setdiff(meta_ids, trans_res$snp))>0){
    print(setdiff(meta_ids, trans_res$snp))
  }
  # filter out FDR 5% results
  trans_res_filt = trans_res %>% filter(adj.pval <= 0.05, snp %in% meta_ids) %>% data.frame()
  if (nrow(trans_res_filt)>0){
    trans_res_filt$qtl_group = qtl_group
    res = rbind(res, trans_res_filt)
  }
}

genes_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/genes_metadata.tsv", header = T, sep = "\t")
genes_metadata = unique(genes_metadata[,c("gene_id", "gene_name",  "chromosome", "gene_start", "gene_end")])

res2 = merge(res, genes_metadata, by.x = "gene_id", by.y = "gene_id")

res2 = res2 %>% mutate(is_trans = (chromosome != chr) | (chromosome == chr) & abs(gene_start - start) > 5000000)
res2 = res2 %>% mutate(is_cis = (chromosome == chr) & (start >= gene_start) & (start <= gene_end)) # is cis when snp is inside a gene

write.table(res2, "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/full_trans_res.tsv", sep = "\t", quote = F)

# summarize number of trans genes per qtl group
trans_summary = res2 %>% select(snp, gene_id, is_trans, qtl_group) %>% distinct() %>% group_by(snp, qtl_group) %>% summarise(nr_of_trans = sum(is_trans), full_nr = n())
trans_summary2 = trans_summary %>% filter(nr_of_trans > 0)

# filter out cis genes (+/- 5 mp) 
res_filt2 = res2 %>% filter(is_trans)

write.table(res_filt2, file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/filtered_trans_res.tsv", sep = "\t", quote = F)

res_filt2 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/filtered_trans_res.tsv", sep = "\t", stringsAsFactors = F)

##### Get % of genes in clusters #####
full_clusters  = list() # genes in clusters
for(cl in unique(metares$cl_id2)){
  print(cl)
  subdf = metares %>% filter(cl_id2 == cl) %>% select(approach, cl_method, qtl_group, cl_id, cl_id2) %>% distinct()
  clpath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(subdf$approach[1], "_coexpr"), subdf$cl_method[1], "clusters" )
  if (subdf$approach[1] == "integrated"){
    clfile = list.files(path = clpath, pattern = "*_clusters.tsv")
    cldata = read.table(file.path(clpath,clfile), sep = "\t", header = T, stringsAsFactors = F)
    full_clusters[[cl]] = cldata %>% filter(make.names(cl_id) == subdf$cl_id[1]) %>% .$gene_id
  }else{
    for(q in as.character(subdf$qtl_group)){
      clfile = list.files(path = clpath, pattern = "*_clusters.tsv")
      clfile = clfile[grep(x = clfile, pattern = paste0(q, "_"))]
      cldata = read.table(file.path(clpath,clfile), sep = "\t", header = T, stringsAsFactors = F)
      full_clusters[[cl]] = cldata %>% filter(make.names(cl_id) == subdf$cl_id[1]) %>% .$gene_id
    }
  }
}
saveRDS(full_clusters, "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/full_clusters.rds")

full_clusters = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/full_clusters.rds")

metalist <- split(as.character(metares$cl_id2), as.character(metares$meta_id))
all_clusters = list()
for(m in names(metalist)){
  all_clusters[[as.character(m)]] = full_clusters[unique(metalist[[as.character(m)]])]
}

saveRDS(all_clusters, "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters.rds")

all_clusters = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters.rds")

meta_trans_genes = res_filt2 %>% select(snp, qtl_group, gene_id) %>% group_by(snp, qtl_group) %>% nest()

d = split(meta_trans_genes$qtl_group, meta_trans_genes$snp)
d2 = lapply(names(d), function(x) sapply(d[[x]], function(y) meta_trans_genes[meta_trans_genes$snp == x  & meta_trans_genes$qtl_group == y,][["data"]][[1]][["gene_id"]], USE.NAMES = T))
names(d2) = names(d)

# Find clusters that don't contain any trans related genes
metares$nr_trans_in_cl = 0
for(i in 1:nrow(metares)){
  meta_id = metares[i,][["meta_id"]]
  qtl_group = metares[i,][["qtl_group"]]
  cl_id2 = metares[i,][["cl_id2"]]
  if(meta_id %in% names(d2)){
    if(qtl_group %in% names(d2[[meta_id]]) | qtl_group %in% colnames(d2[[meta_id]])){
      print(meta_id)
      print(qtl_group)
      if(!is.null(colnames(d2[[meta_id]]))){
        nr_genes = length(intersect(d2[[meta_id]][,qtl_group], all_clusters[[meta_id]][[cl_id2]]))
        nr_fdr_genes = length(d2[[meta_id]][,qtl_group])
      }else{
        nr_genes = length(intersect(d2[[meta_id]][[qtl_group]], all_clusters[[meta_id]][[cl_id2]]))
        nr_fdr_genes = length(d2[[meta_id]][[qtl_group]])
      }
      metares[i,"nr_trans_in_cl"] = nr_genes
      
      # add nr of FDR5% trans genes from gene-level analysis
      metares[i, "nr_gene_trans"] = nr_fdr_genes

      # add cluster size
      metares[i, "cl_size"] = length(all_clusters[[meta_id]][[cl_id2]])
      
      # check if snp overlapping gene is in cluster
      pos = as.numeric(strsplit(meta_id, "_")[[1]][2])
      chr = as.numeric(gsub("chr", "", strsplit(meta_id, "_")[[1]][1]))
      cis_genes = genes_metadata %>% filter(chromosome == chr, (pos >= gene_start) & (pos <= gene_end))
      if(nrow(cis_genes) == 0){
        metares[i, "cis_in_cl"] = FALSE
      }else{
        metares[i, "cis_in_cl"] = any(cis_genes$gene_id %in% all_clusters[[meta_id]][[cl_id2]])
      }
      # check if significant cis gene in the cluster
      sign_cis_genes = res2[res2$snp == meta_id & res2$qtl_group == qtl_group & res2$is_cis == TRUE,]
      metares[i, "sign_cis_in_cl"] = any(sign_cis_genes$gene_id %in% all_clusters[[meta_id]][[cl_id2]])
    }
  }
}

# exclude rows with clusters where no genes are in trans association
metares_filt = metares %>% filter(nr_trans_in_cl > 0)

# statistics about discarded relations
diff_stats = metares %>% filter(nr_trans_in_cl == 0)
diff_stats = diff_stats %>% select(cl_method, approach, cl_id2, meta_id) %>% distinct()
# how many cluster-meta_id pairs excluded as cis-relation?
summ_table = table(diff_stats$cl_method, diff_stats$approach)
# how many cluster-meta_id pairs in total?
total_table = metares %>% select(cl_method, approach, cl_id2, meta_id) %>% distinct()
total_table = table(total_table$cl_method, total_table$approach)

# How many meta_ids were excluded alltogether?
# 298 out of 601 meta_ids excluded
length(setdiff(unique(metares$meta_id), unique(metares_filt$meta_id)))

# add the filter column to the res file
saveRDS(metares, file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
write.table(metares, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/all_crediblesets.tsv", sep = "\t", quote = F)

# additional comparisons after reviews
metares = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/all_crediblesets.tsv", sep = "\t")
metares_filt = metares %>% filter(nr_trans_in_cl > 0)

# one-sided Fisher’s exact test
#fisher.test(matrix(c(74, 1000, 844, 16465), nrow = 2), alternative = "greater") ARHGEF3 example

# nr_trans_in_cl = 19. column
# cl_size = 21. column
# nr_gene_trans = 20. column
pvals = apply(metares_filt, 1, 
      function(x) {
        tbl <- matrix(c(as.numeric(x[19]), 
                        as.numeric(x[21])-as.numeric(x[19]), 
                        as.numeric(x[20])-as.numeric(x[19]), 
                                   18383-as.numeric(x[21])-(as.numeric(x[20])-as.numeric(x[19]))), nrow = 2)
        fisher.test(tbl, alternative="greater")$p.value
      })
metares_filt$fisher_pval = pvals

p = ggplot(metares_filt %>% select(cl_id2, cl_method, approach, meta_id, fisher_pval, cl_size, nr_gene_trans, nr_trans_in_cl) %>% distinct()) + 
  geom_histogram(aes(x=fisher_pval), bins = 20) + theme_bw() + facet_grid(cl_method~approach, scales = "free_y") + 
  xlab("P-value (one-sided Fisher’s exact test)") + ylab("Frequency")
ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS_fisher.png", p, height = 6)

write.table(metares_filt, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/filt_crediblesets.tsv", sep = "\t", quote = F)
metares_filt = read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/filt_crediblesets.tsv", sep = "\t", stringsAsFactors = F)

# filter with Bonferroni
adjusted_pvals = metares_filt %>% select(meta_id, cl_id2, fisher_pval) %>% distinct() %>% group_by(meta_id) %>% 
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = "bonferroni"))

metares_filt2 = base::merge(metares_filt, adjusted_pvals, by = c("meta_id", "cl_id2", "fisher_pval"), all.x = TRUE, sort = FALSE)
write.table(metares_filt2, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/filt_crediblesets.tsv", sep = "\t", quote = F)

#metares_filt3 = metares_filt2 %>% filter(fisher_pval_adj<0.05) # leaves us with 247 meta_ids out of 303

# add the new filter column to the res file
metares2 = merge(metares, metares_filt2[, c("meta_id", "cl_id2", "snp_id", "trans_lead", "cs_id", "fisher_pval", "fisher_pval_adj")], by = c("meta_id", "cl_id2", "snp_id", "trans_lead", "cs_id"), sort = FALSE, all.x = TRUE, all.y = FALSE)

saveRDS(metares2, file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
write.table(metares2, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/all_crediblesets.tsv", sep = "\t", quote = F)

