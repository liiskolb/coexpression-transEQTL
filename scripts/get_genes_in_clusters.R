# Extract cluster genes 
library(picaplot)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggthemes)

sd.outlier <- function(x) abs(x-mean(x))/sd(x) > 2

# 1. Integrated clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)

exprs_mat <- read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/Merged_ENSG_expression.tsv", sep = "\t", header = T, stringsAsFactors = F)
genenames = row.names(exprs_mat)

for (m in methods){
  if (m == "funcExplorer"){
    next
  }else {
    clfiles = list.files(file.path(resdir, m, "clusters"), pattern = "*.rds", full.names = T)
    for (cf in clfiles){
      if(grepl("PLIERpaths", cf)){
        next
      }
      print(cf)
      clres = readRDS(cf)
      if(m == "ICA"){
        #ic_covar_mx = picaplot::getCovariateMx(clres)
        ic_covar_mx = t(clres$A)
        cl_ids = colnames(ic_covar_mx)
        # extreme genes (2 std-s from the mean (distribution has zero mean))
        sds = apply(clres$S, 2, function(x)(abs(x)>2*stats::sd(x)))
        sds = sds[,cl_ids]
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "PEER"){
        weights <- clres$weights[,-1] # first column is average cluster
        row.names(weights) = genenames
        colnames(weights) = paste0("F", seq(1, dim(weights)[2]))
        sds = apply(weights, 2, sd.outlier)
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "PLIER"){
        Z = clres$Z
        colnames(Z) = rownames(clres$B)
        sds = apply(Z, 2, sd.outlier)
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "WGCNA"){
        cl_colors <- paste0("ME", clres$colors)
        cldf = data.frame("gene_id" = genenames, "cl_id" = cl_colors)
        # grey is the module of genes that are not in any of the clusters
        cldf = cldf[cldf$cl_id!="MEgrey",]
      }
      group = gsub(".rds", "", basename(cf))
      resfile = file.path(resdir, m, "clusters", paste0(group,"_full_clusters.tsv"))
      write.table(cldf, resfile, sep = "\t", quote = F, row.names = F)
    }
  }
}


# 2. Separate clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)


for (m in methods){
  if (m == "funcExplorer"){
    next
  }else {
    clfiles = list.files(file.path(resdir, m, "clusters"), pattern = "*.rds", full.names = T)
    for (cf in clfiles){
      if(grepl("PLIERpaths", cf)){
        next
      }
      print(cf)
      group = gsub(".rds", "", basename(cf))
      # read in cluster results
      clres = readRDS(cf)
      # read in the expression matrix
      exprfile = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs", paste0(gsub(paste0(m, "_"), "", group), ".tsv"))
      exprs_mat <- read.table(exprfile, sep = "\t", header = T, stringsAsFactors = F)
      genenames = row.names(exprs_mat)
      if(m == "ICA"){
        #ic_covar_mx = picaplot::getCovariateMx(clres)
        ic_covar_mx = t(clres$A)
        cl_ids = colnames(ic_covar_mx)
        # extreme genes (2 std-s from the mean (distribution has zero mean))
        sds = apply(clres$S, 2, function(x)(abs(x)>2*stats::sd(x)))
        sds = sds[,cl_ids]
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "PEER"){
        weights <- clres$weights[,-1] # first column is average cluster
        row.names(weights) = genenames
        colnames(weights) = paste0("F", seq(1, dim(weights)[2]))
        sds = apply(weights, 2, sd.outlier)
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "PLIER"){
        Z = clres$Z
        colnames(Z) = rownames(clres$B)
        sds = apply(Z, 2, sd.outlier)
        cldf = melt(sds)
        cldf = cldf[cldf$value,c("Var1", "Var2")]
        colnames(cldf) = c("gene_id", "cl_id")
      }else if (m == "WGCNA"){
        cl_colors <- paste0("ME", clres$colors)
        cldf = data.frame("gene_id" = genenames, "cl_id" = cl_colors)
        # grey is the module of genes that are not in any of the clusters
        cldf = cldf[cldf$cl_id!="MEgrey",]
      }
      resfile = file.path(resdir, m, "clusters", paste0(group,"_separate_clusters.tsv"))
      write.table(cldf, resfile, sep = "\t", quote = F, row.names = F)
    }
  }
}


## Find summaries of the cluster sizes and number of clusters
summaries = data.frame()

# Read in sample metadata
sample_metadata = read.table(file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data", "metadata", "Kolberg_2020_duplicate.tsv"), sep = "\t", header = T, stringsAsFactors = F)
# keep only the QC passed samples
sample_metadata = sample_metadata %>% filter(rna_qc_passed == T, genotype_qc_passed == T)

  
# 1. Integrated clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)
#methods = methods[methods!="funcExplorer_old"]

for(m in methods){
  clfiles = list.files(file.path(resdir, m, "clusters"), pattern = "*_clusters.tsv", full.names = T)
  for (cf in clfiles){
    print(cf)
    clusters = read.table(cf, sep = "\t", header = T, stringsAsFactors = F)
    nr_of_cl = length(unique(clusters$cl_id))
    cl_sizes = table(clusters$cl_id)
    #summaries = rbind(summaries, data.frame("method" = "integrated", "cl_method" = m, "nr_of_cl" = nr_of_cl, "min_cl_size" = min(cl_sizes), "ave_cl_size" = as.integer(mean(cl_sizes)), "max_cl_size" = max(cl_sizes)))
    summaries = rbind(summaries, data.frame("method" = "integrated", "cl_method" = m, "nr_of_cl" = nr_of_cl, "cl_size" = cl_sizes))
    }
}

summaries_stats = summaries %>% group_by(method, cl_method, nr_of_cl) %>% summarise(min_cl_size = min(cl_size.Freq), max_cl_size = max(cl_size.Freq), ave_cl_size = mean(cl_size.Freq))
write.table(summaries_stats, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/integrated_cl_summary.tsv" ,sep = "\t", quote = F, row.names = F)

summaries2 = data.frame()

# 2. Separate clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)
for(m in methods){
  clfiles = list.files(file.path(resdir, m, "clusters"), pattern = "*_clusters.tsv", full.names = T)
  for (cf in clfiles){
    print(cf)
    group = gsub("_separate_clusters.tsv", "", gsub(paste0(m, "_"), "", basename(cf)))
    clusters = read.table(cf, sep = "\t", header = T, stringsAsFactors = F)
    nr_of_cl = length(unique(clusters$cl_id))
    cl_sizes = table(clusters$cl_id)
    summaries2 = rbind(summaries2, data.frame("method" = "separate", "qtl_group" = group,  "cl_method" = m, "nr_of_cl" = nr_of_cl, "cl_size" = cl_sizes))
    #summaries2 = rbind(summaries2, data.frame("method" = "separate", "qtl_group" = group,  "cl_method" = m, "nr_of_cl" = nr_of_cl, "min_cl_size" = min(cl_sizes), "ave_cl_size" = as.integer(mean(cl_sizes)), "max_cl_size" = max(cl_sizes)))
  }
}

summaries3 = merge(summaries2, unique(sample_metadata[, c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
summaries3$qtl_group_label = paste(summaries3$cell_type_label, summaries3$condition_label, sep = "_")
summaries3 = summaries3[,!colnames(summaries3) %in% c("cell_type_label", "condition_label")]

summaries3_stats = summaries3 %>% group_by(qtl_group, qtl_group_label, method, cl_method, nr_of_cl) %>% summarise(min_cl_size = min(cl_size.Freq), max_cl_size = max(cl_size.Freq), ave_cl_size = mean(cl_size.Freq))
write.table(summaries3_stats, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/separate_cl_summary.tsv" ,sep = "\t", quote = F, row.names = F)

## Plot cluster sizes

colors = list("B cell_naive" = "#59A14F", "CD4+ T cell_naive" = "#E15759", "CD8+ T cell_naive" = "#ED9899",
              "neutrophil_naive" = "#F28E2B", "platelet_naive" = "#797979", "monocyte_naive" = "#38D4D6",
              "monocyte_LPS_2h" = "#00A2FF", "monocyte_LPS_24h" = "#417DD6",
              "monocyte_IFNg_24h" = "#4F4EA7")

intplot = ggplot(summaries) + geom_boxplot(aes(x = cl_method, y = cl_size.Freq), fill = "white", color = "maroon") + xlab("") + 
  ylab("Cluster size") + scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), limits = c(0,3200)) + 
  theme_pubr(base_size = 20) + 
  geom_text(data = summaries %>% select(cl_method, nr_of_cl) %>% distinct(), aes(x = cl_method, y = 3100, label = nr_of_cl), size = 6) + 
  #labs(subtitle = "Integrated clustering") + 
  theme(axis.ticks.x = element_blank()) 
intplot

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_cluster_stats1.png", intplot, width = 7, height = 5)

### Plot number of clusters

intplot2 = ggplot(summaries %>% select(cl_method, nr_of_cl) %>% distinct()) + 
  geom_point(aes(x = cl_method, y = nr_of_cl), size = 6, color = "salmon4") + 
  geom_segment( aes(x = cl_method, y = nr_of_cl-5, xend=cl_method, yend=0)) + 
  xlab("") + ylab("Number of clusters") + 
  theme_pubr(base_size = 20) +
  theme(axis.ticks.x = element_blank())
  
intplot2
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_cluster_stats1_1.png", intplot2, width = 7, height = 5)


summaries3$qtl_group_label = factor(summaries3$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
sepplot = ggplot(summaries3, aes(x = qtl_group_label, y = as.numeric(cl_size.Freq), 
                                 fill = qtl_group_label)) + 
  geom_boxplot(position = position_dodge(width = 1)) + xlab("") + ylab("Cluster size") + 
  theme_pubr(base_size = 20, legend = "right") + 
  geom_text(data = summaries3 %>% select(cl_method, qtl_group_label, nr_of_cl) %>% distinct(), color = "black", aes(label = nr_of_cl, y = 3100), size = 5, position = position_dodge(width = 1)) + 
  facet_grid(.~cl_method, switch="both") +
  theme(axis.text.x = element_blank(), 
        #axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 16), 
        strip.placement = "outside") + 
  #labs(subtitle = "Separate clustering") +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), limits = c(0,3200)) + 
  scale_fill_manual(values = unlist(colors), name = "") 
#+ scale_fill_manual(values = unlist(colors), name = "")
sepplot

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Suppl_cluster_stats2.png", sepplot, width = 15, height = 6.5)

### Merge plots into one ###

patch =  (wrap_elements(sepplot + ggtitle("Cell types separately"))) / (wrap_elements(intplot + ggtitle("Integrated")) + plot_spacer()) + 
  plot_layout(heights = c(4,3)) 
patch
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS1.png", patch, width = 24, height = 10)

### Separate clustering, number of clusters distribution
summaries4 = summaries3 %>% select(qtl_group_label, cl_method, nr_of_cl) %>% distinct()
sepplot2 = ggplot(summaries4, aes(x = cl_method, y = as.numeric(nr_of_cl), color = qtl_group_label)) + 
  geom_point(size = 6) + theme_pubr(base_size = 20) + ylab("Number of clusters") + xlab("") + 
  scale_color_manual(values = unlist(colors), guide = F) + labs(color = "Group") + 
  theme(axis.ticks.x = element_blank())

sepplot2

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_cluster_stats2.png", sepplot2, width = 7, height = 5)



### Sample sizes barplots ####
sample_metadata$qtl_group_label = paste(sample_metadata$cell_type_label, sample_metadata$condition_label, sep = "_")
sample_metadata$study2 = unlist(lapply(sample_metadata$sample_id, function(x) {
  if(any(grep("IPC", x))){return("CEDAR")}
  if(any(grep("SAMPLE_", x))){return("Kasela_2017")}
  if(any(grep("b_", x))){return("Fairfax_2012")}
  if(any(grep("CD16_", x))){return("Naranbhai_2015")}
  else{return("Fairfax_2014")}}))

sample_summary = sample_metadata %>% dplyr::group_by(qtl_group, qtl_group_label) %>% dplyr::summarize(nr_samples = dplyr::n())
sample_summary = sample_summary %>% arrange(-nr_samples)
sample_summary$qtl_group_label = factor(sample_summary$qtl_group_label, levels = sample_summary$qtl_group_label)

sample_summary2 = sample_metadata %>% group_by(qtl_group, qtl_group_label, study2) %>% dplyr::summarize(nr_samples = dplyr::n())
sample_summary2 = sample_summary2 %>% arrange(-nr_samples)
sample_summary2$qtl_group_label = factor(sample_summary2$qtl_group_label, levels = sample_summary$qtl_group_label)


#samplesizep = ggplot(sample_summary) + geom_bar(aes(x = qtl_group_label, y = nr_samples, fill = qtl_group_label), stat = "identity", width=0.5) + 
#  theme_pubr(base_size = 14) + xlab("") + ylab("Number of samples") + 
#  labs(fill = "Condition") + scale_fill_manual(values = unlist(colors))
#samplesizep

#ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_sample_sizes1.pdf", samplesizep, width = 15, height = 6)
#ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_sample_sizes1.png", samplesizep, width = 15, height = 6)

library(plyr)

samplesizep2 = ggplot(sample_summary2, aes(x = factor(study2), y = nr_samples, fill = qtl_group_label)) + 
  geom_bar(stat = "identity", width=0.5) + 
  theme_pubr(base_size = 20) + xlab("") + ylab("Number of samples") + 
  labs(fill = "Condition") + scale_fill_manual(values = unlist(colors), guide = F) + 
  geom_text(aes(label=nr_samples, y = nr_samples), position = position_stack(vjust = 0.5, reverse = FALSE), color="white", size=8)
samplesizep2

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_sample_sizes2.pdf", samplesizep2, width = 12, height = 6)
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_sample_sizes2.png", samplesizep2, width = 12, height = 6)

## MDS plot

mat_norm = read.table(file.path("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/", "Merged_ENSG_expression.tsv"), sep = "\t", header = T, stringsAsFactors = F)
expr_matrix = as.matrix(mat_norm)

# Perform MDS
dist = cor(expr_matrix, method = "pearson")
fit <- MASS::isoMDS(1 - dist, k=2)

mds_matrix = as.data.frame(fit$points) %>%
  as_tibble() %>%
  dplyr::mutate(sample_id = rownames(fit$points)) %>%
  dplyr::left_join(as.data.frame(sample_metadata), by = "sample_id")

mds_plot = ggplot2::ggplot(mds_matrix, aes(x = V1, y = V2)) + 
  geom_point(aes(color = qtl_group_label), size = 3, alpha = 0.6) + 
  scale_color_manual(values = unlist(colors), guide = F) + 
  xlab("MDS Coordinate 1") + ylab("MDS Coordinate 2") + 
  theme_pubr(base_size = 20) 
  #labs(color = "Condition")  
  #xlim(-0.3, 0.6) + 
  #ylim(-0.3, 0.6) 

mds_plot

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_MDS.pdf", mds_plot, width = 7, height = 6)
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_MDS.png", mds_plot, width = 7, height = 6)


## Number of clusters with trans-eQTLs pval<5.0e-8
library(readr)

# 1. Integrated clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)
filtsummary1 = data.frame()

for (m in methods){
  filtdir = file.path(resdir, m, "eQTLres", "filtered")
  filtfiles = list.files(filtdir, full.names= T, pattern = "*.tsv")
  subdf = data.frame()
  for (f in filtfiles){
    print(f)
    group = gsub("_full_eigen_eQTLres_filtered.tsv", "", gsub(paste0(m, "_"), "", basename(f)))
    trans_filtered = readr::read_delim(f, delim = "\t")
    nr_cl = length(unique(trans_filtered$cluster)) # nr of clusters with at least on genome-wide trans relation
    filtsummary1 = rbind(filtsummary1, data.frame("method" = "integrated", "qtl_group" = group,  "cl_method" = m, "nr_cl_with_trans" = nr_cl))
  }
}

filtsummary1 = merge(filtsummary1, unique(sample_metadata[, c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
filtsummary1$qtl_group_label = paste(filtsummary1$cell_type_label, filtsummary1$condition_label, sep = "_")
filtsummary1 = filtsummary1[,!colnames(filtsummary1) %in% c("cell_type_label", "condition_label")]
filtsummary4 = merge(filtsummary1[,c("qtl_group", "qtl_group_label", "method", "cl_method", "nr_cl_with_trans")], summaries, by.x = c("method", "cl_method"), by.y = c("method", "cl_method"))

write.table(filtsummary4, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/integrated_cl_trans_summary.tsv" ,sep = "\t", quote = F, row.names = F)

# 2. Separate clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)
filtsummary2 = data.frame()

for (m in methods){
  filtdir = file.path(resdir, m, "eQTLres", "filtered")
  filtfiles = list.files(filtdir, full.names= T, pattern = "*.tsv")
  subdf = data.frame()
  for (f in filtfiles){
    print(f)
    group = gsub("_separate_eigen_eQTLres_filtered.tsv", "", gsub(paste0(m, "_"), "", basename(f)))
    trans_filtered = readr::read_delim(f, delim = "\t")
    nr_cl = length(unique(trans_filtered$cluster)) # nr of clusters with at least on genome-wide trans relation
    filtsummary2 = rbind(filtsummary2, data.frame("method" = "separate", "qtl_group" = group,  "cl_method" = m, "nr_cl_with_trans" = nr_cl))
  }
}
filtsummary2 = merge(filtsummary2, unique(sample_metadata[, c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
filtsummary2$qtl_group_label = paste(filtsummary2$cell_type_label, filtsummary2$condition_label, sep = "_")
filtsummary2 = filtsummary2[,!colnames(filtsummary2) %in% c("cell_type_label", "condition_label")]

filtsummary3 = merge(summaries3, filtsummary2[,c("qtl_group", "method", "cl_method", "nr_cl_with_trans")], by.x = c("qtl_group", "method", "cl_method"), by.y = c("qtl_group", "method", "cl_method"))
write.table(filtsummary3, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/separate_cl_trans_summary.tsv" ,sep = "\t", quote = F, row.names = F)

## Plot the statistics
library(ggplot2)

plotdata1 = filtsummary4 %>% mutate(prop = nr_cl_with_trans/nr_of_cl)
plotdata2 = filtsummary3 %>% mutate(prop = nr_cl_with_trans/nr_of_cl)
plotdata = rbind(plotdata1, plotdata2[,colnames(plotdata1)])
p <- ggplot(plotdata) + geom_bar(aes(x = qtl_group_label, y = prop, fill = cl_method), stat="identity", position=position_dodge()) + theme_bw() +
  theme(text = element_text(size=14), axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  xlab("") + ylim(0,1) + ylab("Proportion") + facet_grid(method~.) + ggtitle("Proportion of clusters with at least one lead SNP\nwith genome-wide threshold <5x10^-8") + 
  labs("fill" = "Clustering method")
ggsave(p, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eqtl_proportions.pdf", width = 8, height = 6)
ggsave(p, file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eqtl_proportions.png", width = 8, height = 6)
