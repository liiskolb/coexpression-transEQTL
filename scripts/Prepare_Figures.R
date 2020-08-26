### Helpers ####

library(dplyr)
library(ggplot2)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(gprofiler2)
library(Rsamtools)
library(scales)
library(tidyr)
library(data.table)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library("wesanderson")
library(grid)

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

transManhattan <- function(trans_resfile1, cl_id1, param, bigres, meta_snp1, cis_lead_snp){
  transres <- data.frame(scanTabixDataFrame(trans_resfile1, param, col_names = FALSE))
  colnames(transres) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
  transres1 <- transres %>% filter(cluster == cl_id1)
  # credible set
  csres1 = bigres %>% filter(cl_id == cl_id1, meta_id == meta_snp1) 
  plottitle = paste0("trans-eQTL for ", csres1$cl_id[1], " (", csres1$cl_method[1], ", ", csres1$approach[1], ")")
  trans_lead1 = csres1 %>% select(trans_lead) %>% distinct() %>% .$trans_lead
  cs_snps1 = csres1 %>% select(snp_id) %>% distinct() %>% .$snp_id
  
  highlight_snps = c(cis_lead_snp, meta_snp1, trans_lead1)
  transres1$show_label = transres1$snp %in% highlight_snps
  transres1$highlight = transres1$snp %in% cs_snps1
  
  # create plots
  p1 <- ggplot(transres1) + geom_point(aes(x = start, y = -log10(pvalue)), 
                                       alpha = 0.6, size = 3, color = "#8c8c8c") +
    theme_pubr(base_size = 30) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  # highlight credible set
  p1 <- p1 + geom_point(data = transres1 %>% filter(highlight == T), aes(x = start, y = -log10(pvalue)), size = 3, colour = "#00b159")
  # show text labels for leads
  p1 <- p1 + geom_point(data = transres1 %>% filter(show_label == T), 
                        aes(x = start, y = -log10(pvalue)), 
                        size = 6, shape = 21, colour = "black", fill = "#d11141") + 
    geom_text(data = transres1 %>% filter(show_label == T), 
              aes(x = start, y = -log10(pvalue), label = snp), 
              size = 10, nudge_y = 2, nudge_x = 10000)
  p1 <- p1 + labs(subtitle = plottitle)
  return(p1)
}

publish_gosttable2 <- function(gostres, highlight_terms = NULL, use_colors = TRUE, rot = 90, show_columns = c("source", "term_name", "term_size", "intersection_size"), filename = NULL, fontsize = 32){
  # gostres is the GOSt response list (contains results and metadata) or a data frame
  term_id <- p_values <- query <- p_value <- NULL
  
  if (class(gostres) == "list"){
    if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
    df <- gostres$result
  } else if (class(gostres) == "data.frame"){
    df <- gostres
  } else {
    stop("The input 'gostres' should be a data frame or a list from the gost() function.")
  }
  
  # make sure all the essential column names are there
  if (!"term_id" %in% colnames(df)) stop("The required column 'term_id' is missing from the input")
  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the input")
  
  # selected terms
  if (is.null(highlight_terms)) {
    # show full table if no terms given
    highlight_terms = df
  }
  
  if (is.data.frame(highlight_terms)){
    message("The input 'highlight_terms' is a data.frame. The column 'term_id' will be used.")
    if ("term_id" %in% colnames(highlight_terms)){
      highlight_terms <- highlight_terms$term_id
    }
    else{
      stop("No column named 'term_id'.")
    }
  }
  
  subdf <- base::subset(df, term_id %in% highlight_terms)
  
  if (nrow(subdf) == 0){
    stop("None of the term IDs in the 'highlight_terms' were found from the results.")
  }
  
  subdf$id <- match(subdf$term_id, highlight_terms)
  
  # default column names to show
  show_columns <- unique(append(show_columns, c("term_id", "p_value")))
  gp_colnames <- c("source", "term_id", "term_name", "term_size", "query_size", "intersection_size", "p_value", "intersection_sizes", "query_sizes")
  
  colnames <- gp_colnames[which(gp_colnames %in% show_columns)]
  
  # include non gprofiler columns
  if (length(setdiff(show_columns, gp_colnames)) > 0){
    colnames <- append(colnames, setdiff(show_columns, gp_colnames))
  }
  
  # If multiquery
  if ("p_values" %in% colnames(subdf)){
    if ("meta" %in% names(gostres)){
      meta <- gostres$meta
      subdf$query <- list(names(meta$query_metadata$queries))
    } else {
      qnames = paste("query", seq(1, length(subdf$p_values[[1]])), sep = "_")
      subdf$query <- list(qnames)
    }
    spread_col = c("p_values")
    if ("query_sizes" %in% show_columns){
      spread_col = append(spread_col, "query_sizes")
    }
    if ("intersection_sizes" %in% show_columns){
      spread_col = append(spread_col, "intersection_sizes")
    }
    # spread the data frame to correct form
    subdf <- tidyr::unnest(data = subdf, cols = c(spread_col, query))
    subdf <- dplyr::rename(subdf, p_value = p_values)
    subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
    showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
    #showdf <- subdf[,stats::na.omit(match(c(colnames, "query", spread_col), names(subdf)))]
    showdf <- tidyr::pivot_wider(showdf, names_from = query, values_from = c(p_value, spread_col[spread_col!="p_values"]))
  } else {
    if ("query" %in% names(subdf) & length(unique(subdf$query)) > 1){
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
      showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
      spread_col <- c("p_value", "intersection_size", "query_size")
      spread_col <- intersect(colnames(showdf), spread_col)
      showdf <- tidyr::pivot_wider(showdf, names_from = query, values_from = spread_col, names_prefix = ifelse(length(spread_col) == 1,"pval ", ""))
      
    } else {
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
      showdf <- subdf[,stats::na.omit(match(colnames, names(subdf)))]
    }
  }
  # find the column to color
  idx <- which(grepl(pattern = "val", x = names(showdf)))
  # remove pval from the name
  names(showdf)[idx] = unlist(lapply(names(showdf)[idx], function(x) gsub(pattern = "pval ", replacement= "", x = x)))
  
  # Prepare table
  if (use_colors){
    colours <- matrix("white", nrow(showdf), ncol(showdf))
    cols <- sapply(showdf[,idx], function(x) mapViridis(-log10(as.numeric(x))))
    colours[,idx] <- cols
    fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
    fontcolours[,idx] <- "white"
  } else {
    order_of_cl = names(showdf)[idx]
    # add empty columns next to pvalue columns that show the color scale (to the right)
    showdf[paste0(1, order_of_cl)] <- NA
    # update the order with extra pvalue color code column
    order_of_cl2 = c(rbind(paste0(1, order_of_cl), order_of_cl))
    showdf = showdf[,c(names(showdf)[1:min(idx)-1], order_of_cl2)]
    colours <- matrix("white", nrow(showdf), ncol(showdf))
    # add colors to empty columns
    temp_df = showdf[,order_of_cl2]
    temp_cols <- sapply(temp_df, function(x) ifelse(!is.na(x), mapViridis(-log10(as.numeric(x))), "white"))
    # switch values
    temp_cols[,seq(1,ncol(temp_cols),2)] = temp_cols[,seq(2,ncol(temp_cols),2)]
    temp_cols[,seq(2,ncol(temp_cols),2)] = "white"
    colours[,which(names(showdf) %in% order_of_cl2)] <- temp_cols
    fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
    # remove column names from color scale column
    names(showdf)[startsWith(names(showdf), "1")] = ""
  }
  
  fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
  
  showdf[is.na(showdf)] <- ""
  # change the column names of term_id and term_name
  names(showdf)[names(showdf) == "term_id"] = "Term ID"
  names(showdf)[names(showdf) == "term_name"] = "Term name"
  th <- gridExtra::ttheme_default(base_size = fontsize,
                                  padding = grid::unit(c(5, 5), "mm"),
                                  core=list(
                                    padding.h = grid::unit(c(5,5), "mm"),
                                    padding.v = grid::unit(c(5,5), "mm"),
                                    bg_params = list(fill = colours, col="black", lwd = 0.5),
                                    fg_params=list(hjust = 0, x = 0.01, col=fontcolours, fontface=fontfaces)),
                                  colhead=list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"),
                                               fg_params=list(col="black", fontface="bold", rot = rot)),
                                  rowhead=list(fg_params=list(col="black", fontface="bold")))
  
  #vp <- viewport(width = 0.9) 
  tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
  #tb$widths[-1] <- rep(unit(1, "null"), 2)
  #tb$heights <- rep(unit(1/nrow(tb), "null"), nrow(tb))
  
  #h <- grid::unit.c(sum(tb$heights))
  #w <- grid::unit.c(sum(tb$widths))
  tg <- gridExtra::grid.arrange(tb, ncol = 1)
  
  #p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
  
  if (is.null(filename)){
    return(tg)
  } else{
    imgtype <- strsplit(basename(filename), split="\\.")[[1]][-1]
    
    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }
    
    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", "bmp")){
      width = grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.2
      height = grid::convertHeight(sum(tg$heights), "in", TRUE) + 0.2
      ggplot2::ggsave(filename = filename, plot = p, height = height, width = width)
      message("The image is saved to ", filename)
      return(p)
    } else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}

### Load data and define scales ####
colors2 = list("CL_0000236_naive" = "#59A14F", "CL_0000624_naive" = "#E15759", "CL_0000625_naive" = "#ED9899",
               "CL_0000775_naive" = "#F28E2B", "CL_0000233_naive" = "#797979", "CL_0002057_naive" = "#38D4D6",
               "CL_0002057_LPS_2h" = "#00A2FF", "CL_0002057_LPS_24h" = "#417DD6",
               "CL_0002057_IFNg_24h" = "#4F4EA7")

colors = list("B cell_naive" = "#59A14F", "CD4+ T cell_naive" = "#E15759", "CD8+ T cell_naive" = "#ED9899",
              "neutrophil_naive" = "#F28E2B", "platelet_naive" = "#797979", "monocyte_naive" = "#38D4D6",
              "monocyte_LPS_2h" = "#00A2FF", "monocyte_LPS_24h" = "#417DD6",
              "monocyte_IFNg_24h" = "#4F4EA7")

cs_colors = list("lead SNP" = "#d11141", "cis cs" = "chocolate1", "trans cs" = "#00b159")

all_clusters = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_clusters.rds")
respath = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets"
res = readRDS(file.path(respath, "all_crediblesets.rsd"))
# exclude cis clusters
res_filt = res %>% filter(nr_trans_in_cl > 0)

# cis credible sets
ciscs = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/cis_cs.tsv", sep = "\t", header = T, stringsAsFactors = F)

sample_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/Kolberg_2020_duplicate.tsv", sep = "\t", header = T, stringsAsFactors = F)
genes_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/genes_metadata.tsv", header = T, sep = "\t")
genes_metadata = unique(genes_metadata[,c("gene_id", "gene_name",  "chromosome", "gene_start", "gene_end")])
genes_metadata = genes_metadata %>% filter(!chromosome %in% c("MT"))

exprs = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/Merged_ENSG_expression_long.tsv", sep = "\t", stringsAsFactors = F, header = T)

eigens_sep = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/full_eigens_separate.tsv", sep = "\t", stringsAsFactors = F)
eigens_sep = merge(eigens_sep, unique(sample_metadata[,c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
eigens_sep$qtl_group_label = paste(eigens_sep$cell_type_label, eigens_sep$condition_label, sep = "_")
eigens_sep$qtl_group_label = factor(eigens_sep$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
eigens_sep$variable = make.names(eigens_sep$variable)

eigens = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/full_eigens_integrated.tsv", sep = "\t", stringsAsFactors = F)
eigens = merge(eigens, unique(sample_metadata[,c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
eigens$qtl_group_label = paste(eigens$cell_type_label, eigens$condition_label, sep = "_")
eigens$qtl_group_label = factor(eigens$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))

### Fig 1 ####

# Genotypes PCA - data from Nurlan Kerimov
pcares = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/pop_assign/projections_comb.tsv", sep = "\t", header = T, stringsAsFactors = F)

pcp = ggplot(pcares %>% filter(superpopulation_code!="AMR")) + 
  geom_point(aes(x = PC1, y = PC2, color = superpopulation_code), 
      size = 4, alpha = 0.6) + 
  theme_pubr(base_size = 30, legend = "right") + xlab("PC1") + ylab("PC2") + 
  coord_fixed()+
  #theme(aspect.ratio=1) + 
  scale_color_brewer("Population", palette = "Dark2", type = "qual") + 
  guides(colour = guide_legend(override.aes = list(alpha=1))) 

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_PCA.png", pcp, width = 10, height = 10)

# MDS plot

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

mds_plot
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_MDS.png", mds_plot, width = 7, height = 6)

# Meta credible sets Manhattan plot
# nominal results 
res = readRDS(file = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets", "all_crediblesets.rsd"))
res_filt = res %>% filter(nr_trans_in_cl > 0)
sample_metadata = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/metadata/Kolberg_2020_duplicate.tsv", sep = "\t", header = T, stringsAsFactors = F)

group_data = sample_metadata %>% select(qtl_group, cell_type_label, condition_label) %>% distinct() %>% mutate(qtl_group_label = paste(cell_type_label, condition_label, sep = "_")) %>% select(qtl_group, qtl_group_label)

effectsres <- read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", header = T, stringsAsFactors = F)
resgenome = res_filt %>% select(cl_id, cl_method, approach, qtl_group, lead_pval, meta_id) %>% distinct()

resgenome$meta_chr = unlist(lapply(resgenome$meta_id, function(x) gsub("chr", "", strsplit(x, "_")[[1]][1])))
resgenome$meta_pos = unlist(lapply(resgenome$meta_id, function(x) as.numeric(strsplit(x, "_")[[1]][2])))
resgenome = merge(resgenome, group_data, by = "qtl_group")

resgenome$meta_chr = factor(resgenome$meta_chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
resgenome = merge(resgenome, effectsres, by.x = c("meta_id", "cl_id", "cl_method", "approach", "qtl_group"), by.y = c("meta_id", "cl_id", "cl_method", "approach", "qtl_group"))

resfdr2.tmp = resgenome %>% dplyr::group_by(meta_chr) %>%
  dplyr::summarise(chr_len = max(meta_pos)) %>%
  # cumulative position
  dplyr::mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # add this info to the initial dataset
  dplyr::left_join(resgenome, ., by = c("meta_chr" = "meta_chr")) %>%
  # cumulative position to each gene
  dplyr::arrange(meta_chr, meta_pos) %>%
  dplyr::mutate(meta_poscum = meta_pos + tot)

axis.set = resfdr2.tmp %>% dplyr::group_by(meta_chr) %>% 
  dplyr::summarize(center = (max(meta_poscum) + min(meta_poscum)) / 2, start = min(meta_poscum), end = max(meta_poscum))
axis.set$meta_chr = factor(axis.set$meta_chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

resfdr2.tmp = resfdr2.tmp %>% arrange(meta_poscum, lead_pval)

metacsp = ggplot(resfdr2.tmp, aes(x = meta_poscum, y = -log10(lead_pval), color = qtl_group_label.y)) + 
  geom_point(size = 15, alpha = 0.8) +
  theme_pubr(base_size = 55) + scale_color_manual(values = unlist(colors), guide = F) + 
  scale_x_continuous(label = axis.set$meta_chr, breaks = axis.set$center, 
                     limits = c(1993370, NA), expand = c(0.01, 0.01)) +
  geom_vline(xintercept = axis.set$start, linetype=2, size = 2, color = "grey") + 
  scale_y_continuous(breaks = c(0, 8, 20, 30, 40, 50, 60, 70), limits = c(0,75),
                     expand = c(0.01,0.01)) +  
  xlab("CHR") + theme(axis.text.x = element_text(size = 45)) 

metacsp2 <- metacsp + ggrepel::geom_text_repel(data = resfdr2.tmp %>% filter(meta_chr == 12, meta_id == "chr12_69344099_A_G") %>% arrange(lead_pval, beta) %>% head(1), label = "LYZ", segment.size = 1.5, segment.color = "black", min.segment.length = 0, color = "black", nudge_y = 2, nudge_x = 100000000, size = 30)
metacsp2 <- metacsp2 + ggrepel::geom_text_repel(data = resfdr2.tmp %>% filter(meta_chr == 9, meta_id == "chr9_20818520_A_G") %>% arrange(lead_pval, beta) %>% head(1), label = "IFNB1", segment.size = 1.5, segment.color = "black", min.segment.length = 0, color = "black", nudge_y = 3, nudge_x = 100000000, size = 30)
metacsp2 <- metacsp2 + ggrepel::geom_text_repel(data = resfdr2.tmp %>% filter(meta_chr == 4, meta_id == "chr4_102325419_ACACT_A") %>% arrange(lead_pval, beta) %>% head(1), label = "SLC39A8", segment.size = 1.5, segment.color = "black", min.segment.length = 0, color = "black", nudge_y = 8, nudge_x = 100000000, size = 30)
metacsp2 <- metacsp2 + ggrepel::geom_text_repel(data = resfdr2.tmp %>% filter(meta_chr == 3, meta_id == "chr3_56815721_T_C") %>% arrange(lead_pval, beta) %>% head(1), label = "ARHGEF3" , segment.size = 1.5, segment.color = "black", min.segment.length = 0, color = "black", nudge_y = 5, nudge_x = -100000000, size = 30)
metacsp2

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_metacs_manhattan.png", metacsp2, width = 45, height = 20, limitsize = F)

# Number of modules

# 1. Integrated clustering
resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr"
methods = list.dirs(resdir, recursive = F, full.names = F)
summaries = data.frame()
for(m in methods){
  clfiles = list.files(file.path(resdir, m, "clusters"), pattern = "*_clusters.tsv", full.names = T)
  for (cf in clfiles){
    print(cf)
    clusters = read.table(cf, sep = "\t", header = T, stringsAsFactors = F)
    nr_of_cl = length(unique(clusters$cl_id))
    cl_sizes = table(clusters$cl_id)
    summaries = rbind(summaries, data.frame("method" = "integrated", "cl_method" = m, "nr_of_cl" = nr_of_cl, "cl_size" = cl_sizes))
  }
}
summaries_stats = summaries %>% group_by(method, cl_method, nr_of_cl) %>% summarise(min_cl_size = min(cl_size.Freq), max_cl_size = max(cl_size.Freq), ave_cl_size = mean(cl_size.Freq))
write.table(summaries_stats, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/integrated_cl_summary.tsv" ,sep = "\t", quote = F, row.names = F)

intplot2 = ggplot(summaries %>% select(cl_method, nr_of_cl) %>% distinct()) +
geom_point(aes(x = cl_method, y = nr_of_cl), size = 6, color = "salmon4") +
geom_segment( aes(x = cl_method, y = nr_of_cl-5, xend=cl_method, yend=0)) +
xlab("") + ylab("Number of modules") +
theme_pubr(base_size = 20) +
theme(axis.ticks.x = element_blank())
intplot2
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_cluster_stats1_1.png", intplot2, width = 7, height = 5)

# 2. Separate clustering
summaries2 = data.frame()
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
  }
}

summaries3 = merge(summaries2, unique(sample_metadata[, c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
summaries3$qtl_group_label = paste(summaries3$cell_type_label, summaries3$condition_label, sep = "_")
summaries3 = summaries3[,!colnames(summaries3) %in% c("cell_type_label", "condition_label")]

summaries3_stats = summaries3 %>% group_by(qtl_group, qtl_group_label, method, cl_method, nr_of_cl) %>% summarise(min_cl_size = min(cl_size.Freq), max_cl_size = max(cl_size.Freq), ave_cl_size = mean(cl_size.Freq))
write.table(summaries3_stats, "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/separate_cl_summary.tsv" ,sep = "\t", quote = F, row.names = F)

summaries4 = summaries3 %>% select(qtl_group_label, cl_method, nr_of_cl) %>% distinct()
sepplot2 = ggplot(summaries4, aes(x = cl_method, y = as.numeric(nr_of_cl), 
                                  color = qtl_group_label)) +
geom_point(size = 6) + theme_pubr(base_size = 20) + ylab("Number of modules") + xlab("") +
scale_color_manual(values = unlist(colors), guide = F) +
theme(axis.ticks.x = element_blank())
sepplot2
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig1_cluster_stats2.png", sepplot2, width = 7, height = 5)

### Fig S1 - module sizes ####

intplot = ggplot(summaries) + geom_boxplot(aes(x = cl_method, y = cl_size.Freq), fill = "white", color = "maroon") + xlab("") + 
  ylab("Cluster size") + scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), limits = c(0,3200)) + 
  theme_pubr(base_size = 20) + 
  geom_text(data = summaries %>% select(cl_method, nr_of_cl) %>% distinct(), aes(x = cl_method, y = 3100, label = nr_of_cl), size = 6) + 
  #labs(subtitle = "Integrated clustering") + 
  theme(axis.ticks.x = element_blank()) 
intplot

summaries3$qtl_group_label = factor(summaries3$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
sepplot = ggplot(summaries3, aes(x = qtl_group_label, y = as.numeric(cl_size.Freq), 
                                 fill = qtl_group_label)) + 
  geom_boxplot(position = position_dodge(width = 1)) + xlab("") + ylab("Cluster size") + 
  theme_pubr(base_size = 20, legend = "right") + 
  geom_text(data = summaries3 %>% select(cl_method, qtl_group_label, nr_of_cl) %>% distinct(), color = "black", aes(label = nr_of_cl, y = 3100), size = 5, position = position_dodge(width = 1)) + 
  facet_grid(.~cl_method, switch="both") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 16), 
        strip.placement = "outside") + 
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), limits = c(0,3200)) + 
  scale_fill_manual(values = unlist(colors), name = "") 

# Merge plots into one 

patch = (wrap_elements(sepplot + ggtitle("Cell types separately"))) / (wrap_elements(intplot + ggtitle("Integrated")) + plot_spacer()) + 
  plot_layout(heights = c(10,8)) 
patch
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS1.png", patch, width = 24, height = 15)


### Fig 2 - ARHGEF3 results ####

meta_snp2 = "chr3_56815721_T_C"
arh_res = res_filt %>% filter(meta_id == meta_snp2) # 3 clusters, only in platelet
table(arh_res$cl_id2)
table(arh_res$qtl_group)
table(arh_res$trans_lead) # 1 lead for all the clusters

# Part 1: Manhattan plots
chr = 3
position = 56815721
pheno_id2 = "ILMN_1781010" # ARHGEF3, ENSG00000163947
width = 100000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))

# cis plot
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0000233_naive.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE))
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% filter(phenotype_id == pheno_id2)
cisres = cisres %>% arrange(pvalue)
head(cisres)
cis_lead_snp2 = cisres %>% filter(is_top == 1) %>% .$snp
cis_p2 <- ggplot(cisres) + geom_point(aes(x = start, y = -log10(pvalue)), size = 3,
               alpha = 0.6, color = "#8c8c8c") +
  theme_pubr(base_size = 30) + theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
cis_p2 <- cis_p2 + labs(subtitle = "cis-eQTL for ARHGEF3 (platelets)")

# cis credible set
arhgef_cis_cs = ciscs %>% filter(phenotype_id == pheno_id2, qtl_group == "CL_0000233_naive")

# trans plot
cl_id = "X6.WIERENGA_STAT5A_TARGETS_DN"
trans_resfile2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/PLIER/eQTLres/PLIER_CL_0000233_naive_separate_eigen_eQTLres.gz"
p2 <- transManhattan(trans_resfile2, cl_id, param, arh_res, meta_snp2, NULL)

p2 <- p2 + xlab("CHR 3 position (bp)") + theme(axis.title.x=element_text(),
               axis.text.x=element_text(),
               axis.ticks.x=element_line())
p2 <- p2 + labs(subtitle = "trans-eQTL for\nX6.WIERENGA_STAT5A_TARGETS_DN \n(PLIER, separate, platelets)")
# add rsid
p2[["layers"]][[4]] = NULL
p2 = p2 + ylim(0,10) + geom_text_repel(data = p2$data %>% filter(snp == meta_snp2),
            aes(x = start, y = -log10(pvalue)), label = "rs1354034",
                 size = 12, nudge_y = 1, nudge_x = 10000)

# highlight lead snp on cis plot
cis_p2 <- cis_p2 + ylim(0,20) + geom_point(data = cis_p2$data %>% filter(snp == meta_snp2),
                              aes(x = start, y = -log10(pvalue)), size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p2$data %>% filter(snp == meta_snp2),
                  aes(x = start, y = -log10(pvalue)), label = "rs1354034",
                  size = 12, nudge_y = 1, nudge_x = 10000)

# highlight cis credible set in trans and cis plot
cis_p2 <- cis_p2 + geom_point(data = cisres %>% filter(snp %in% arhgef_cis_cs$variant_id),
                              aes(x = start, y = -log10(pvalue)), size = 3,
                              colour = "chocolate1")
# highlight trans credible set in cis plot
trans_cs = p2$data %>% filter(highlight == T)
cis_p2 <- cis_p2 + geom_point(data = cisres %>% filter(snp %in% trans_cs$snp),
                              aes(x = start, y = -log10(pvalue)), size = 3, colour = "#00b159")


# Part 2: line plot
effectsres <- read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", header = T, stringsAsFactors = F)
effectsres$qtl_group_label = factor(effectsres$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
effectsres_sub = effectsres %>% filter(meta_id == meta_snp2)
linep2 <- ggplot(effectsres_sub) + geom_line(aes(x = as.factor(qtl_group_label), 
                                                 y = logpval, group = cl_id2)) +
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = cl_method, shape = approach), 
             size = 5, stroke = 2) +
  xlab("") + theme_pubr(base_size = 30, legend = "right") +
  scale_x_discrete(drop=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values=c(24, 8)) +
  scale_y_continuous(limits = c(0,NA)) + ylab("-log10(p-value)")

# Part 3: Manhattan plot

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
f <- file.path(resdir, "CL_0000233_naive_eQTLres.gz")
transres <- readr::read_delim(f, delim = "\t", col_names=F)
colnames(transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
transres = transres %>% group_by(snp) %>% mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
genes_in_arh_cl = all_clusters[[meta_snp2]][["X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive"]]
transres2 = transres %>% filter(snp == meta_snp2)
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_arh_cl
transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))
transres2.tmp = transres2 %>% group_by(chromosome) %>%
  summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  arrange(chromosome, gene_start) %>%
  mutate(gene_startcum = gene_start + tot)

axis.set = transres2.tmp %>% group_by(chromosome) %>%
summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)
axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp2 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: rs1354034, Condition: platelet, Module: X6.WIERENGA_STAT5A_TARGETS_DN")

fp2 <- fp2 + geom_point(data = transres2.tmp %>% filter(highlight == T),
             aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
             color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.001,
    gene_name!="ARHGEF3"),
    aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
    size = 10, point.padding = 0.2, force = 10)

arhdat = transres2.tmp %>% filter(gene_name == "ARHGEF3")
fp2 <- fp2 + geom_point(data = arhdat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = arhdat,
        aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
        color = "red", fontface = "bold", point.padding = 0.5, nudge_x = 10000, nudge_y = 3, size = 10) + 
  ylim(-14,14)
fp2
  
# Part 4: g:Profiler
by_genes = transres2.tmp %>% filter(adj.pval < 0.05) %>% select(gene_id) %>% distinct()
arh_clusters = all_clusters[[meta_snp2]]
# significant by gene eQTLs
arh_clusters[["Gene-level"]] = by_genes$gene_id
names(arh_clusters) = unlist(lapply(names(arh_clusters), function(x)
  ifelse(x %like% "PLIER", paste(strsplit(x,"_")[[1]][1], "\n", length(arh_clusters[[x]]), " genes" , "\n", "PLIER", sep = "") ,
                       ifelse(x %like% "funcExplorer", paste(paste(strsplit(x,"_")[[1]][1:2], collapse = "_"), "\n", length(arh_clusters[[x]]), " genes" , "\n", "funcExplorer", sep =""),
                                            paste(strsplit(x, "_")[[1]][1], "\n", length(arh_clusters[[x]]), " genes" , "\n", strsplit(x, "_")[[1]][2], sep ="")))))
gostres2 <- gost(arh_clusters, significant = F)
highlight_ids2 <- c("GO:0030168", "KEGG:04611", "REAC:R-HSA-76002", "GO:0005925", "GO:0071281")
level_idx = c(which(names(arh_clusters) %like% "X6"), which(!names(arh_clusters) %like% "X6"))
gostres2$result$query = factor(gostres2$result$query, levels = names(arh_clusters)[level_idx])
gostres2$result = gostres2$result %>% arrange(query)
gpt2 <- publish_gosttable2(gostres2, rot = 0, use_colors = F, highlight_terms = highlight_ids2, show_columns = c("term_name"))

# Part 5: association plots
# GWAS plot: filtered the file to region in command line
# gwasfile = "/gpfs/hpc/home/a72094/datasets/summary_stats/GWASCatalog/Astle_2016/AstleWJ_27863252_GCST004599/harmonised/27863252-GCST004599-EFO_0004584.h.tsv.gz"
gwasfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/platelet_volume_arhgef3_gwas.tsv"

gwasdat = read.table(gwasfile, sep = "\t", header = T, stringsAsFactors = F)
gwasdat$hm_variant_id = paste0("chr", gwasdat$hm_variant_id)
# filter region further to 100000bp
gwasdat = gwasdat %>% filter(hm_pos >= 56815721-100000 & hm_pos <= 56815721+100000)

# create plots
gwasp <- ggplot(gwasdat) + geom_point(aes(x = hm_pos, y = mlog10p),
           alpha = 0.6, size = 3, color = "#8c8c8c") +
  theme_pubr(base_size = 30)

gwasp <- gwasp + labs(subtitle = "Mean platelet volume GWAS") +
  xlab("") + ylab("-log10(p_value)") + theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

# Add colocalisation PP4 values
gwasp <- gwasp + annotate("text", x = 56720000, y = 750, hjust = 0, label = "cis PP4: 0.99\ntrans PP4: 0.996", size = 10)

# highlight lead SNP
gwasp <- gwasp + geom_point(data = gwasdat %>% filter( hm_variant_id == meta_snp2),
                            aes(x = hm_pos, y = mlog10p),
                            shape = 21, size = 6, fill = "#d11141") +
  geom_text_repel(data = gwasdat %>% filter( hm_variant_id == meta_snp2),
                  aes(x = hm_pos, y = mlog10p), label = "rs1354034",
                  size = 12, nudge_y = 50, nudge_x = 20000)

# highlight trans cs
gwasp <- gwasp + geom_point(data = gwasp$data %>% filter(hm_variant_id %in% trans_cs$snp),
                            aes(x = hm_pos, y = mlog10p), size = 3, colour = "#00b159")

# Cluster expression 

x6_eigens = eigens_sep %>% filter(variable == "X6.WIERENGA_STAT5A_TARGETS_DN", qtl_group == "CL_0000233_naive")
params = GRanges(c(3), IRanges(start = 56815721, end = 56815721))
tabix_header = Rsamtools::headerTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz")
tabix_header2 = tabix_header$header[229]
tabix_header2 = unlist(lapply(tabix_header2, function(x) strsplit(x, "\t")[[1]]))
tabix_list = Rsamtools::scanTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz", param = params)
lead_genotype = unlist(lapply(tabix_list, function(x) strsplit(x, "\t")[[1]]))
gen_df = data.frame(lead_genotype)
gen_df$genotype_id = tabix_header2
gen_df = gen_df[!gen_df$genotype_id %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),, drop = F]
row.names(gen_df) = NULL
colnames(gen_df) = c(meta_snp2, "genotype_id")
gen_df[[meta_snp2]] = as.character(gen_df[[meta_snp2]])
gen_df[[meta_snp2]] = unlist(lapply(gen_df[[meta_snp2]], function(x) strsplit(x, ":")[[1]][1]))
gen_df[[meta_snp2]] = unlist(lapply(gen_df[[meta_snp2]], function(x) ifelse(x == "0|0", "TT", ifelse(x == "1|1", "CC", "TC")) ))
gen_df[[meta_snp2]] = factor(gen_df[[meta_snp2]], levels = c("TT", "TC", "CC"))
x6_eigens = merge(x6_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp2 <- ggplot(x6_eigens) + geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) +
  facet_wrap(.~qtl_group_label, nrow = 1) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  ylab("Eigengene expression") +
  labs(subtitle = "X6.WIERENGA_STAT5A_TARGETS_DN\n(PLIER, separate, platelets)") + 
  xlab("rs1354034") + 
  theme(plot.margin = margin(t = 1,r = 10,b = 1, l = 1, "cm"))
genp2

# ARHGEF3 expression
arhgef_expr = exprs %>% filter(gene_id == "ENSG00000163947")
arhgef_expr = merge(arhgef_expr, sample_metadata[,c("sample_id", "genotype_id", "qtl_group", "cell_type_label", "condition_label")], by.x = c("variable", "qtl_group"), by.y = c("sample_id", "qtl_group"), all.y = F)
arhgef_expr$qtl_group_label = paste(arhgef_expr$cell_type_label, arhgef_expr$condition_label, sep = "_")
arhgef_expr$qtl_group_label = factor(arhgef_expr$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))

arhgef_expr = merge(arhgef_expr, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp3 <- ggplot(arhgef_expr %>% filter(qtl_group == "CL_0000233_naive")) +
  geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) +
  facet_wrap(.~qtl_group_label, nrow = 1) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  ylab("Gene expression (log2)") + xlab("rs1354034") +
  labs(subtitle = "ARHGEF3")
genp3

genp4 <- ggplot(arhgef_expr) +
  geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) +
  facet_wrap(.~qtl_group_label, nrow = 1) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  ylab("Gene expression (log2)") + xlab("rs1354034") +
  labs(subtitle = "ARHGEF3")
genp4

# create legend
testdata = cis_p2$data[1:3,]
testdata$cs = "cis cs"
testdata$cs[1] = "trans cs"
testdata$cs[3] = "lead SNP"
testdata$cs = factor(testdata$cs)
levels(testdata$cs) = c("cis cs", "trans cs", "lead SNP")

legend_p <- ggplot(testdata) + geom_point(aes(x = 0, y = 0, color = cs), size = 8) +
  scale_color_manual(name = "", values = unlist(cs_colors)) +
  theme_pubr() +
  theme(legend.text = element_text(size = 30), legend.position = "top",
  axis.text.y = element_text(), axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.y = element_blank(), axis.line.y = element_blank(),
  strip.text.x = element_blank(),
  strip.background = element_blank(),
  panel.spacing.x=unit(0, "lines"),
  panel.spacing.y=unit(0, "lines"),
  axis.line = element_blank(),
  legend.box.spacing = unit(-1, 'cm'),
  legend.margin = margin(0, 0, 0, 0, "cm"),
  plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0))
legend_p

leg <- get_legend(legend_p)
th <- sum(leg$heights)
leg2 = gridExtra::grid.arrange(leg, heights = unit.c(unit(1, "null"), th))

# Merge all the plots
part1 <- wrap_elements(plot = wrap_elements(leg2) /gwasp / (cis_p2 + ylim(0,24)) / p2 + plot_layout(heights = c(0.5,5,5,5)))
part2 <- wrap_elements(linep2) / wrap_plots(genp2, genp3, ncol = 2, widths = c(3,3)) + plot_layout(heights = c(4,2))
part3 <- wrap_elements(fp2) / (wrap_elements(gpt2)) + plot_layout(heights = c(10, 6))
patchfigg <- (part1 | part2) / part3 + plot_layout(heights = c(14, 12), widths = c(3, 6))
patchfigg <- patchfigg + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 40, face = "bold"))
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig2.png", patchfigg, 
       width = 32, height = 28)

### Fig 3 - SLC39A8 results ####

meta_snp3 = "chr4_102325419_ACACT_A"
slc_res = res_filt %>% filter(meta_id == meta_snp3) # 1 cluster, only in LPS_24h
table(slc_res$cl_id2)
table(slc_res$qtl_group)
table(slc_res$trans_lead) # 1 lead

# Part 1: Manhattan plots
chr = 4
position = 102325419
pheno_id1 = "ILMN_1695316" # SLC39A8, ENSG00000138821
#pheno_id2 = "ILMN_2233539" # SLC39A8
#pheno_id3 = "ILMN_1714965" # NFKB1
#missense = "chr4_102267552_C_T"
#phenos = c("ILMN_1695316", "ILMN_2233539", "ILMN_1714965")
width = 100000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))


# cis plot naive
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_naive.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE)[[1]])
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% dplyr::filter(phenotype_id == pheno_id1)
cis_lead_snp3 = cisres %>% filter(is_top == 1) 

cis_p3 <- ggplot(cisres) + geom_point(aes(x = start, y = -log10(pvalue)), size = 3,
                      alpha = 0.6, color = "#8c8c8c") +
theme_pubr(base_size = 30) + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
cis_p3 <- cis_p3 + labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)\n", "(monocytes naive)"))

# highlight trans lead
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == meta_snp3),
                           aes(x = start, y = -log10(pvalue)), size = 6, shape = 21, fill = "#d11141") +
  geom_text_repel(data = cis_p3$data %>% filter(snp == meta_snp3),
                             aes(x = start, y = -log10(pvalue)),
                            size = 12, nudge_x = 1000, nudge_y = 1, label = "rs75562818")


# Kim-Hellmuth cis plot LPS 90min

kim90 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/Kim-Hellmuth/nominal_lps90.txt", sep = " ", stringsAsFactors = F)
colnames(kim90) = c("probe_id", "snp", "distance", "pval", "beta")

cisp1 = ggplot(kim90 %>% filter(probe_id == pheno_id1, distance <= 64014+width & distance >= -width + 64014)) +
  geom_point(aes(x = distance, y = -log10(pval)), size = 3, alpha = 0.6, color = "#8c8c8c") +
  theme_pubr(base_size = 30) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)", "\nKim-Hellmuth et al. (monocytes LPS 90min)"))

slc_sub = res_filt %>% filter(meta_id == "chr4_102325419_ACACT_A") %>% select(snp_id, rs_id) %>% distinct()
slc_sub$rs_id[slc_sub$snp_id == "chr4_102325419_ACACT_A"] = "rs75562818"
slc_sub$rs_id[slc_sub$snp_id == "chr4_102340060_T_A"] = "rs62329461"

kim90[kim90$snp == "rs145021182","snp"] = "rs75562818"
kim90 = merge(kim90, slc_sub, by.x = "snp", by.y = "rs_id", all.x = T)

cisp1 = cisp1 + geom_point(data = kim90 %>% filter(probe_id == pheno_id1, snp_id == meta_snp3, 
                                                   distance <= 64014+width & distance >= -width + 64014),
             aes(x = distance, y = -log10(pval)), size = 6, shape = 21, fill = "#d11141") +
  ylim(0, 10) +
  geom_text_repel(data = kim90 %>% filter(probe_id == pheno_id1, snp_id == meta_snp3, 
                                          distance <= 64014+width & distance >= -width + 64014),
        aes(x = distance, y = -log10(pval), label = snp), size = 12, nudge_x = 1000, nudge_y = 1) + 
  ylab("-log10(pvalue")

# trans plot
cl_id = "Cluster_10413"
trans_resfile3 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/funcExplorer/eQTLres/funcExplorer_CL_0002057_LPS_24h_separate_eigen_eQTLres.gz"

p <- transManhattan(trans_resfile3, cl_id, param, slc_res, meta_snp3, NULL)

p <- p + xlab("CHR 4 position (bp)") + theme(axis.title.x=element_text(),
                                               axis.text.x=element_text(),
                                               axis.ticks.x=element_line())
p <- p + labs(subtitle = "trans-eQTL for Cluster_10413\n(funcExplorer, separate, monocytes LPS 24h)")

# add rs-id
p[["layers"]][[4]] = NULL
p = p + ylim(0,10) + geom_text_repel(data = p$data %>% filter(snp == meta_snp3),
                         aes(x = start, y = -log10(pvalue)),
                      label = "rs75562818",
                      size = 12, nudge_y = 1, nudge_x = 10000)

# highlight trans credible set
trans_cs = p$data %>% filter(highlight == T)
cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% trans_cs$snp),
  aes(x = start, y = -log10(pvalue)), 
  size = 3, colour = "#00b159")

cisp1 <- cisp1 + geom_point(data = cisp1$data %>% filter(snp %in% slc_sub$rs_id),
                            aes(x = distance, y = -log10(pval)), 
                            size = 3, colour = "#00b159" )

# highlight cis lead
cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == cis_lead_snp3$snp),
                              aes(x = start, y = -log10(pvalue)),
                              size = 6, shape = 21, fill = "#d11141") + 
        geom_text_repel(data = cis_p3$data %>% filter(snp == cis_lead_snp3$snp),
                        aes(x = start, y = -log10(pvalue)),
                        size = 12, nudge_x = 1000, nudge_y = 1, label = "rs11097779")

cisp1 <- cisp1 + geom_point(data = cisp1$data %>% filter(snp == "rs11097779"), 
                            aes(x = distance, y = -log10(pval)),
                            size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cisp1$data %>% filter(snp == "rs11097779"),
                  aes(x = distance, y = -log10(pval)),
                  size = 12, nudge_x = 1000, nudge_y = 1, label = "rs11097779")

p <- p + geom_point(data = p$data %>% filter(snp == cis_lead_snp3$snp),
                    aes(x = start, y = -log10(pvalue)),
                    size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p$data %>% filter(snp == cis_lead_snp3$snp),
                  aes(x = start, y = -log10(pvalue)),
                  size = 12, nudge_x = 1000, nudge_y = 1, label = "rs11097779")


# highlight cis credible set

slc_cis_cs = ciscs %>% filter(phenotype_id == pheno_id1, qtl_group == "CL_0002057_naive")

# get rs-ids
dbsnp = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/dbSNP_b151_GRCh38p7.vcf.gz"
rs_ids <- data.frame(scanTabixDataFrame(dbsnp, param, col_names = FALSE, col_types = readr::cols(.default = "c")))
colnames(rs_ids) = c("chr", "pos", "rs_id", "ref", "alt", "qual", "filter", "info")
rs_ids$chr = paste("chr", rs_ids$chr, sep = "")
rs_ids$snp_id = paste(rs_ids$chr, rs_ids$pos, rs_ids$ref, rs_ids$alt, sep = "_")

rs_ids = rs_ids %>% select(snp_id, rs_id, ref, alt, pos)

slc_cis_cs = merge(slc_cis_cs, rs_ids, by.x = "variant_id", by.y = "snp_id", all.x = T) # map to rsids

slc_cis_cs[slc_cis_cs$variant_id == "chr4_102362400_G_A", "rs_id"] = "rs56202913"
slc_cis_cs[slc_cis_cs$variant_id == "chr4_102376963_G_T", "rs_id"] = "rs185667882"
slc_cis_cs[slc_cis_cs$variant_id == "chr4_102377181_G_A", "rs_id"] = "rs62329499"
slc_cis_cs[slc_cis_cs$variant_id == "chr4_102380473_C_T", "rs_id"] = "rs4437253"
slc_cis_cs[slc_cis_cs$variant_id == "chr4_102381096_G_A", "rs_id"] = "rs4594728"
slc_cis_cs[slc_cis_cs$variant_id == "chr4_102381689_C_T", "rs_id"] = "rs7695249"

cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% slc_cis_cs$variant_id),
                     aes(x = start, y = -log10(pvalue)), size = 3, colour = "chocolate1")

cisp1 <- cisp1 + geom_point(data = cisp1$data %>% filter(snp %in% slc_cis_cs$rs_id),
                            aes(x = distance, y = -log10(pval)), size = 3, colour = "chocolate1")

p <- p + geom_point(data = p$data %>% filter(snp %in% slc_cis_cs$variant_id),
                  aes(x = start, y = -log10(pvalue)), size = 3, colour = "chocolate1")

# Part 2: line plot
effectsres <- read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", header = T, stringsAsFactors = F)
effectsres$qtl_group_label = factor(effectsres$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
effectsres_sub = effectsres %>% filter(meta_id == meta_snp3)

linep3 <- ggplot(effectsres_sub) + geom_line(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2)) + 
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = cl_method, shape = approach), 
             size = 5, stroke = 2) +
  xlab("") + theme_pubr(base_size = 30, legend = "right") + 
  scale_x_discrete(drop=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_shape_manual(values=c(24, 8)) + 
  scale_y_continuous(limits = c(0,NA)) + ylab("-log10(p-value)")

# Part 3: Manhattan plot

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
f <- file.path(resdir, "CL_0002057_LPS_24h_eQTLres.gz")
transres <- readr::read_delim(f, delim = "\t", col_names=F)
colnames(transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
transres = transres %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))

genes_in_slc_cl = all_clusters[[meta_snp3]]
transres2 = transres %>% filter(snp == meta_snp3)

# add gene metadata
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_slc_cl[[1]]

transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))

transres2.tmp = transres2 %>% dplyr::group_by(chromosome) %>%
  dplyr::summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  dplyr::mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  dplyr::left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  dplyr::arrange(chromosome, gene_start) %>%
  dplyr::mutate(gene_startcum = gene_start + tot) 

axis.set = transres2.tmp %>% dplyr::group_by(chromosome) %>% 
  dplyr::summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)

axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp3 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  ylim(-10,10)+
  theme(legend.position="none") +
  labs(subtitle = "SNP: rs75562818, Condition: monocytes LPS 24h, Module: Cluster_10413")

fp3 <- fp3 + geom_point(data = transres2.tmp %>% filter(highlight == T),
            aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
           color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.05,
                 gene_name!="SLC39A8"),
                aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                size = 10, point.padding = 0.5, force = 10)
# highlight cis eqtl
slcdat = transres2.tmp %>% filter(gene_name == "SLC39A8")
fp3 <- fp3 + geom_point(data = slcdat, aes(x = gene_startcum,
                 y = -log10(pvalue)*sign(beta)),
                 size = 6, stroke = 3, color = "red", shape = 21) +
  geom_text_repel(data = slcdat,
                   aes(x = gene_startcum, y = -log10(pvalue)*sign(beta),
                       label = gene_name), 
                  color = "red", fontface = "bold", point.padding = 0.5, size = 10)

# Part 4: g:Profiler
by_genes = transres2.tmp %>% filter(adj.pval < 0.05) %>% select(gene_id) %>% distinct()
slc_clusters = all_clusters[[meta_snp3]]

# significant by gene eQTLs
slc_clusters[["Gene-level"]] = by_genes$gene_id

names(slc_clusters) = unlist(lapply(names(slc_clusters), function(x)
 ifelse(x %like% "PLIER", paste(strsplit(x,"_")[[1]][1], "\n", length(slc_clusters[[x]]), " genes" , "\n", "PLIER", sep = "") ,
     ifelse(x %like% "funcExplorer", paste(paste(strsplit(x,"_")[[1]][1:2], collapse = "_"), "\n", length(slc_clusters[[x]]), " genes" , "\n", "funcExplorer", sep =""),
     paste(strsplit(x, "_")[[1]][1], "\n", length(slc_clusters[[x]]), " genes" , "\n", strsplit(x, "_")[[1]][2], sep ="")))))

gostres3 <- gost(slc_clusters)
highlight_ids3 <- c("GO:0008270", "REAC:R-HSA-5661231", "GO:0010273", "GO:0006882", "WP:WP3529", "TF:M00650_0")
gpt3 <- publish_gosttable2(gostres3, rot = 0, use_colors = F, highlight_terms = highlight_ids3, show_columns = c("term_name"))

# Genotype plot

slc_eigens = eigens_sep %>% filter(variable == "Cluster_10413", qtl_group == "CL_0002057_LPS_24h")
params = GRanges(c(4), IRanges(start = 102325419, end = 102325419))
tabix_header = Rsamtools::headerTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz")
tabix_header2 = tabix_header$header[229]
tabix_header2 = unlist(lapply(tabix_header2, function(x) strsplit(x, "\t")[[1]]))
tabix_list = Rsamtools::scanTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz", param = params)
lead_genotype = unlist(lapply(tabix_list, function(x) strsplit(x, "\t")[[1]]))
gen_df = data.frame(lead_genotype)
gen_df$genotype_id = tabix_header2
gen_df = gen_df[!gen_df$genotype_id %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),, drop = F]
row.names(gen_df) = NULL
colnames(gen_df) = c(meta_snp3, "genotype_id")
gen_df[[meta_snp3]] = as.character(gen_df[[meta_snp3]])
gen_df[[meta_snp3]] = unlist(lapply(gen_df[[meta_snp3]], function(x) strsplit(x, ":")[[1]][1]))
gen_df[[meta_snp3]] = unlist(lapply(gen_df[[meta_snp3]], function(x) ifelse(x == "0|0", "0", ifelse(x == "1|1", "2", "1")) ))
gen_df[[meta_snp3]] = factor(gen_df[[meta_snp3]], levels = c("0", "1", "2"))
slc_eigens = merge(slc_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")
genp <- ggplot(slc_eigens) +
  geom_boxplot(aes(x = chr4_102325419_ACACT_A, y = value, fill = qtl_group_label),
         width = 0.5, show.legend = FALSE) +
    facet_wrap(.~qtl_group_label, nrow = 1) +
    scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
    ylab("Eigengene expression") +
    labs(subtitle = "Cluster_10413\n(funcExplorer, separate, monocytes LPS 24h)") +
    xlab("rs75562818") + 
  theme(plot.margin = margin(1,10,1,1,"cm"))

# Merge plots
part23 <- wrap_elements(linep3 + ylim(0,10)) / wrap_plots(genp + plot_spacer() + plot_layout(widths = c(3,2))) + plot_layout(heights = c(4,2))
part13 <- wrap_elements(plot = wrap_elements(leg2) / (cis_p3 + ylim(0,19)) / cisp1 / p + ylim(0,12) + plot_layout(heights = c(0.5,5,5,5)))
part33 <- wrap_elements(fp3) / (wrap_elements(gpt3)) + plot_layout(heights = c(10, 10))
patchfigg3 <- (part13 | part23) / part33 + plot_layout(heights = c(14, 12), widths = c(3, 6))
patchfigg3 <- patchfigg3 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 40, face = "bold"))

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig3.png", patchfigg3, width = 30, height = 25)


### Fig 4.1 - LD heatmap ####
library(LDlinkR)
library(gaston)
library(LDheatmap)
require(snpStats)
require(grid)

rs_ids = c("rs75562818", "rs7692921", "rs6850893", "rs7664683", "rs147218183", "rs62329460", "rs2054394",
"rs62329462", "rs4146610", "rs62329461", "rs13107325", "rs11097779")
rs_pos = gsnpense(rs_ids) # Ghr38

ldmat = LDmatrix(rs_ids, pop = "EUR", token = "cdc6c8e6308d")
row.names(ldmat) = ldmat$RS_number
ldmat$RS_number = NULL

snp_pos = rs_pos[rs_pos$rs_id %in% row.names(ldmat),c("rs_id", "start")]
row.names(snp_pos) = snp_pos$rs_id
snp_pos = snp_pos[row.names(ldmat),]

png("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/SLC_LD_plot.png", height = 650, units = "px")
LDheatmap(gdat = as.matrix(ldmat),
  color = rev(brewer.pal(8, "Reds")), title = "", flip = T,
  geneMapLocation = 0,
  SNP.name = c("rs13107325", "rs75562818", "rs11097779", "rs7692921"))

# remove physical length text
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp=gpar(cex=0))
grid.edit(gPath("ldheatmap", "Key", "title"), label = "")
# remove the gene map lines
grid.edit("diagonal", gp=gpar(cex=0, col = "white"))
grid.edit("segments", gp=gpar(cex=0, col = "white"))
# change symbols
grid.edit("symbols", pch = 20, gp = gpar(cex = 1, col = "black"))
# edit font size
grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(col = "black", hjust = 20))
# add grid
grid.edit(gPath("ldheatmap", "heatMap", "heatmap"), gp = gpar(col = "white", lwd = 1.5))
dev.off()

### Fig 4.2 - motif logos ####
library("seqLogo")
# MTF1
m = read.table("http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/PWMs/Files/M02922_2.00.txt")
m$Pos = NULL
pw <- makePWM(t(m))
png("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig4_2.png", width = 700, height = 300)
seqLogo(pw, yaxis = F)
dev.off()

# NFKB
m2 = read.table("http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/PWMs/Files/M08006_2.00.txt", header = T)
m2$Pos = NULL
pw2 <- makePWM(t(m2))
png("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig4_3.png", width = 700, height = 300)
seqLogo(pw2, yaxis = F)
dev.off()

### Fig 4.3 - SLC expression ####
#exprs = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/data/exprs/Merged_ENSG_expression_long.tsv")
slc_expr = exprs %>% filter(gene_id == "ENSG00000138821")
### mean foldchange
slc_expr %>% group_by(qtl_group_label) %>% summarise(avg = mean(value))

slc_p = ggplot(slc_expr %>% filter(qtl_group %like% "CL_0002057")) +
  geom_boxplot(aes(x = qtl_group_label, y = value, fill = qtl_group_label), show.legend = F) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  labs(subtitle = "SLC39A8") +
  ylab("Gene expression (log2)") + xlab("monocytes") +
  scale_x_discrete(labels = c("naive", "LPS 2h", "LPS 24h", "IFNg 24h")) +
  ylim(0, 16)
slc_p
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/Fig5_2.png", slc_p, width = 9, height = 8)

### Fig S2 - ARHGEF3 example about credible sets ####
meta_snp2 = "chr3_56815721_T_C"

cl_ids = c("IC68_ICA_integrated", "X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive", "Cluster_12953_funcExplorer_separate_CL_0000233_naive")
crediblesets = readRDS("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/all_crediblesets.rsd")
arhgef_crediblesets = crediblesets %>% filter(meta_id == meta_snp2)

arhgef_crediblesets$pos = as.factor(arhgef_crediblesets$cl_id)
arhgef_crediblesets$highlight = arhgef_crediblesets$rs_id == "rs1354034"

csp = ggplot(arhgef_crediblesets) + geom_point(aes(x = snp_pos, y = as.factor(cl_method), color = cl_id), size = 8) +
  geom_point(aes(x = snp_pos, y =  as.factor(cl_method)), size = 8, pch = 21) + 
  geom_text_repel(data = arhgef_crediblesets %>% filter(!highlight), 
                  aes(x = snp_pos, y = as.factor(cl_method), label = rs_id),
                  nudge_y = 0.3, segment.colour = "grey50") +
  geom_text_repel(data = arhgef_crediblesets %>% filter(highlight), 
                  aes(x = snp_pos, y = as.factor(cl_method), label = rs_id),
                  nudge_y = 0.5, nudge_x = 1000,segment.colour = "grey50", color = "red", fontface = "bold") +
  theme_pubr(base_size = 20, legend = "top") + 
  
  theme(axis.text.y = element_text(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.title = element_blank()) + 
  scale_x_continuous(limits = c(56814971-100, 56846416+100)) + 
  xlab("CHR 3 position (bp)") + 
  scale_color_brewer(palette="Set1") + labs(subtitle = "ARHGEF3 credible sets")
csp
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS2.png", csp, width = 10, height = 5)  

### Fig S3 - IFNB1 results ####

meta_snp1 = "chr9_20818520_A_G"
chr = 9
position = 20818520
pheno_id1 = "ILMN_1682245" # IFNB1
width = 300000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))

ifnb_res = res_filt %>% filter(meta_id == meta_snp1) # 12 clusters, only in LPS_24h

# cis plot
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_LPS_2h.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE))
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% filter(phenotype_id == pheno_id1)
cisres %>% arrange(pvalue) %>% head
cis_lead_snp = cisres %>% filter(is_top == 1) %>% .$snp

cis_p <- ggplot(cisres) + geom_point(aes(x = start, y = -log10(pvalue)),
                              alpha = 0.6, size = 3, color = "#8c8c8c") +
  theme_pubr(base_size = 30) + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(), axis.ticks.x=element_blank())
cis_p <- cis_p + labs(subtitle = "cis-eQTL for IFNB1 (monocytes LPS 2h)")

# cis credible set
ifnb_cis_cs = ciscs %>% filter(phenotype_id == pheno_id1, qtl_group == "CL_0002057_LPS_2h")

# example trans plot
cl_id4 = "IC140"
trans_resfile4 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0002057_LPS_24h_full_eigen_eQTLres.gz"
p4 <- transManhattan(trans_resfile4, cl_id4, param, ifnb_res, meta_snp1, NULL)
p4 <- p4 + xlab("CHR 9 position (bp)") + theme(axis.title.x=element_text(),
                                               axis.text.x=element_text(),
                                               axis.ticks.x=element_line()) + 
  labs(subtitle = "trans-eQTL for IC140\n(ICA, integrated, monocytes LPS 24h)")

p4$layers[[4]] = NULL
p4 <- p4 + geom_text_repel(data = p4$data %>% filter(snp == meta_snp1), 
                           aes(x = start, y = -log10(pvalue), label = snp),
                           size = 12, nudge_y = 2, nudge_x = 10000)

# highlight cis lead on cis plot
cis_p <- cis_p + geom_point(data = cis_p$data %>% filter(snp == cis_lead_snp),
                            aes(x = start, y = -log10(pvalue)), size = 6, 
                            shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p$data %>% filter(snp == cis_lead_snp),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 1, nudge_x = 10000)


# highlight trans lead on cis plot
cis_p <- cis_p + geom_point(data = cis_p$data %>% filter(snp == meta_snp1),
                                           aes(x = start, y = -log10(pvalue)), size = 6, 
                            shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p$data %>% filter(snp == meta_snp1),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 2, nudge_x = 10000)

# highlight trans credible set in cis plot
trans_cs = p4$data %>% filter(highlight == T)
cis_p <- cis_p + geom_point(data = cisres %>% filter(snp %in% trans_cs$snp),
                              aes(x = start, y = -log10(pvalue)), size = 3, colour = "#00b159")

# highlight cis lead on trans plot

p4 <- p4 + geom_point(data = p4$data %>% filter(snp == cis_lead_snp),
                            aes(x = start, y = -log10(pvalue)), size = 6, 
                            shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p4$data %>% filter(snp == cis_lead_snp),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 1, nudge_x = 10000)

# Part 2: line plot
effectsres <- read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", header = T, stringsAsFactors = F)
effectsres$qtl_group_label = factor(effectsres$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
effectsres_sub = effectsres %>% filter(meta_id == meta_snp1)
linep <- ggplot(effectsres_sub) + geom_line(aes(x = as.factor(qtl_group_label), 
                                                 y = logpval, group = cl_id2)) +
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, color = cl_method, shape = approach), 
             size = 5, stroke = 2) +
  xlab("") + theme_pubr(base_size = 30, legend = "right") +
  scale_x_discrete(drop=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values=c(24, 8)) +
  scale_y_continuous(limits = c(0,NA)) + ylab("-log10(p-value)")

# Part 3: Manhattan plot

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
f <- file.path(resdir, "CL_0002057_LPS_24h_eQTLres.gz")
transres <- readr::read_delim(f, delim = "\t", col_names=F)
colnames(transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
transres = transres %>% group_by(snp) %>% mutate(adj.pval = p.adjust(pvalue, method = "fdr"))

genes_in_ifnb_cl = all_clusters[[meta_snp1]][["IC140_ICA_integrated"]]
transres2 = transres %>% filter(snp == meta_snp1)
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_ifnb_cl
transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))
transres2.tmp = transres2 %>% group_by(chromosome) %>%
  summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  arrange(chromosome, gene_start) %>%
  mutate(gene_startcum = gene_start + tot)

axis.set = transres2.tmp %>% group_by(chromosome) %>%
  summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)
axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: chr9_20818520_A_G, Condition: monocytes LPS 24h, Module: IC140")

fp <- fp + geom_point(data = transres2.tmp %>% filter(highlight == T),
                        aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.00001,
                                                  gene_name!="IFNB1"),
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                  size = 10, point.padding = 0.1, force = 10)

ifndat = transres2.tmp %>% filter(gene_name == "IFNB1")
fp <- fp + geom_point(data = ifndat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = ifndat,
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
                  color = "red", fontface = "bold", point.padding = 0.5, size = 10)
fp

# Part 4: g:Profiler enrichment
by_genes = transres2.tmp %>% filter(adj.pval < 0.05) %>% select(gene_id) %>% distinct()
ifnb_clusters = all_clusters[[meta_snp1]]

# significant by gene eQTLs
ifnb_clusters[["Gene-level"]] = by_genes$gene_id

names(ifnb_clusters) = unlist(lapply(names(ifnb_clusters), function(x) 
  ifelse(x %like% "PLIER", paste(strsplit(x,"_")[[1]][1], "\n", length(ifnb_clusters[[x]]), " genes" , "\n", "PLIER", sep = "") , 
         ifelse(x %like% "funcExplorer", paste(paste(strsplit(x,"_")[[1]][1:2], collapse = "_"), "\n", length(ifnb_clusters[[x]]), " genes" , "\n", "funcExplorer", sep =""),
                paste(strsplit(x, "_")[[1]][1], "\n", length(ifnb_clusters[[x]]), " genes" , "\n", strsplit(x, "_")[[1]][2], sep ="")))))

gostres <- gost(ifnb_clusters)

highlight_ids <- c("GO:0051607", "GO:0034340", "TF:M11685_0", "REAC:R-HSA-913531", "WP:WP619")

level_idx = c(which(names(ifnb_clusters) == "IC140\n445 genes\nICA"), which(names(ifnb_clusters) != "IC140\n445 genes\nICA"))
gostres$result$query = factor(gostres$result$query, 
                              levels = names(ifnb_clusters)[level_idx])


gostres$result = gostres$result %>% arrange(query)
gpt <- publish_gosttable2(gostres, highlight_terms = highlight_ids, use_colors = F, rot = 90, show_columns = c("term_name"), fontsize = 26)

# Part 5: eigengene expression
eigens = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/full_eigens_integrated.tsv", sep = "\t", stringsAsFactors = F)
eigens = merge(eigens, unique(sample_metadata[,c("qtl_group", "cell_type_label", "condition_label")]), by.x = "qtl_group", by.y = "qtl_group")
eigens$qtl_group_label = paste(eigens$cell_type_label, eigens$condition_label, sep = "_")
eigens$qtl_group_label = factor(eigens$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))

ic140_eigens = eigens %>% filter(variable == "IC140")

params = GRanges(c(9), IRanges(start = 20818520, end = 20818520))
tabix_header = Rsamtools::headerTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz")
tabix_header2 = tabix_header$header[229]
tabix_header2 = unlist(lapply(tabix_header2, function(x) strsplit(x, "\t")[[1]]))
tabix_list = Rsamtools::scanTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz", param = params)
lead_genotype = unlist(lapply(tabix_list, function(x) strsplit(x, "\t")[[1]]))
gen_df = data.frame(lead_genotype)
gen_df$genotype_id = tabix_header2
gen_df = gen_df[!gen_df$genotype_id %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),, drop = F]
row.names(gen_df) = NULL
colnames(gen_df) = c(meta_snp1, "genotype_id")
gen_df[[meta_snp1]] = as.character(gen_df[[meta_snp1]])
gen_df[[meta_snp1]] = unlist(lapply(gen_df[[meta_snp1]], function(x) strsplit(x, ":")[[1]][1]))
gen_df[[meta_snp1]] = unlist(lapply(gen_df[[meta_snp1]], function(x) ifelse(x == "0|0", "AA", ifelse(x == "1|1", "GG", "AG")) ))
gen_df[[meta_snp1]] = factor(gen_df[[meta_snp1]], levels = c("AA", "AG", "GG"))

ic140_eigens = merge(ic140_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp <- ggplot(ic140_eigens) + geom_boxplot(aes(x = chr9_20818520_A_G, y = value, fill = qtl_group_label), 
                                            show.legend = FALSE) + 
  facet_wrap(.~qtl_group_label, nrow = 1) + 
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) + 
  ylab("Eigengene expression") + labs(subtitle = "IC140 (ICA, integrated)")

genp

# Merge all the plots together
part1 <- wrap_elements(plot = (wrap_elements(leg2) / (cis_p + ylim(0, 10)) / (p4 + ylim(0,16))) + plot_layout(heights = c(0.5,5,5)))
part2 <- wrap_elements(linep)
part3 <- wrap_elements(fp) / wrap_elements(genp) / wrap_elements(gpt)
patchfigg2 <- ((part1 | part2) + plot_layout(widths = c(2,3))) / wrap_plots(part3) + plot_layout(heights = c(6,14))
patchfigg2 <- patchfigg2 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 40, face = "bold"))
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS3.png", patchfigg2, 
       width = 38, height = 28)

### Fig S4 - LYZ results ####

meta_snp4 = "chr12_69344099_A_G"

lyz_res = res_filt %>% filter(meta_id == meta_snp4)  # monocytes, 16 different leads, 58 clusters
table(lyz_res$cl_id2)
table(lyz_res$qtl_group)
table(lyz_res$trans_lead) 

# Part 1: Manhattan plots, monocytes naive
chr = 12
position = 69344099
pheno_id4 = "ILMN_1815205" # LYZ
#width = 100000
width = 200000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))

# cis plot
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_naive.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE))
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% filter(phenotype_id == pheno_id4)
cisres = cisres %>% arrange(pvalue)
cis_lead_snp4 = cisres %>% filter(is_top == 1) %>% .$snp

cis_p4 <- ggplot(cisres) + geom_point(aes(x = start, y = -log10(pvalue)), 
                                      alpha = 0.6, size = 3, color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p4 <- cis_p4 + labs(subtitle = "cis-eQTL for LYZ (monocytes naive)")

# cis credible set
lyz_cis_cs = ciscs %>% filter(phenotype_id == pheno_id4, qtl_group == "CL_0002057_naive")

# trans plots
cl_id = "IC86"

# IFNg 24h
trans_resfile4 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0002057_IFNg_24h_full_eigen_eQTLres.gz"
p41 <- transManhattan(trans_resfile4, cl_id, param, lyz_res %>% filter(qtl_group == "CL_0002057_IFNg_24h"), meta_snp4, NULL)
p41 <- p41 + labs(subtitle = "trans-eQTL for IC86 (ICA, integrated, monocytes IFNg 24h)") 
p41[["layers"]][[4]] = NULL

# LPS 2h
#trans_resfile42 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0002057_LPS_2h_full_eigen_eQTLres.gz"
#p42 <- transManhattan(trans_resfile42, cl_id, param, lyz_res %>% filter(qtl_group == "CL_0002057_LPS_2h"), meta_snp4, NULL)

# LPS 24h
#trans_resfile43 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0002057_LPS_24h_full_eigen_eQTLres.gz"
#p43 <- transManhattan(trans_resfile43, cl_id, param, lyz_res %>% filter(qtl_group == "CL_0002057_LPS_24h"), meta_snp4, NULL)

# naive
trans_resfile44 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0002057_naive_full_eigen_eQTLres.gz"
p44 <- transManhattan(trans_resfile44, cl_id, param, lyz_res %>% filter(qtl_group == "CL_0002057_naive"), meta_snp4, NULL)

# add x-axis labels and title
p44 <- p44 + xlab("CHR 12 position (bp)") + labs(subtitle = "trans-eQTL for IC86 (ICA, integrated, monocytes naive)") + 
  theme(axis.title.x=element_text(), axis.text.x=element_text(), axis.ticks.x=element_line())
p44[["layers"]][[4]] = NULL

# highlight cis credible set
cis_p4 <- cis_p4 + geom_point(data = cisres %>% filter(snp %in% lyz_cis_cs$variant_id), 
                              aes(x = start, y = -log10(pvalue)), size = 3, colour = "chocolate1")

# highlight cis lead
cis_p4 <- cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == cis_lead_snp4),
                              aes(x = start, y = -log10(pvalue)), size = 6, 
                              shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == cis_lead_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 3, nudge_x = 70000)

# highlight trans lead
cis_p4 <- cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == meta_snp4),
                              aes(x = start, y = -log10(pvalue)), size = 6, 
                              shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == meta_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 5, nudge_x = -70000) + ylim(0,50)

p41 <- p41 + geom_point(data = p41$data %>% filter(snp == meta_snp4),
                    aes(x = start, y = -log10(pvalue)), size = 6, 
                    shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p41$data %>% filter(snp == meta_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 5, nudge_x = -70000) 

p44 <- p44 + geom_point(data = p44$data %>% filter(snp == meta_snp4),
                        aes(x = start, y = -log10(pvalue)), size = 6, 
                        shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p44$data %>% filter(snp == meta_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 5, nudge_x = -70000) 

# highlight cis lead on trans plot

p41 <- p41 + geom_point(data = p41$data %>% filter(snp == cis_lead_snp4),
                      aes(x = start, y = -log10(pvalue)), size = 6, 
                      shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p41$data %>% filter(snp == cis_lead_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 7, nudge_x = 70000) + ylim(0,70)

p44 <- p44 + geom_point(data = p44$data %>% filter(snp == cis_lead_snp4),
                        aes(x = start, y = -log10(pvalue)), size = 6, 
                        shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p44$data %>% filter(snp == cis_lead_snp4),
                  aes(x = start, y = -log10(pvalue), label = snp),
                  size = 12, nudge_y = 7, nudge_x = 70000) + ylim(0,60)

# highlight trans cs on cis plot
trans_cs = p41$data %>% filter(highlight == T)
cis_p4 <- cis_p4 + geom_point(data = cisres %>% filter(snp %in% trans_cs$snp),
                            aes(x = start, y = -log10(pvalue)), size = 3, colour = "#00b159")


  
# Part 2: Line plots

effectsres <- read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/all_crediblesets/effectsres_filt.tsv", sep = "\t", header = T, stringsAsFactors = F)
effectsres$qtl_group_label = factor(effectsres$qtl_group_label, levels = c("B cell_naive", "CD4+ T cell_naive", "CD8+ T cell_naive", "monocyte_naive", "monocyte_LPS_2h", "monocyte_LPS_24h", "monocyte_IFNg_24h", "neutrophil_naive", "platelet_naive"))
effectsres_sub = effectsres %>% filter(meta_id == meta_snp4)

linep4 <- ggplot(effectsres_sub) + geom_line(aes(x = as.factor(qtl_group_label), 
                                                 y = logpval, group = cl_id2)) + 
  geom_point(aes(x = as.factor(qtl_group_label), y = logpval, group = cl_id2, 
                 color = cl_method, shape = approach), size = 5, stroke = 2) +
  xlab("") + theme_pubr(base_size = 30, legend = "right") + 
  scale_x_discrete(drop=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_shape_manual(values=c(24, 8)) + 
  scale_y_continuous(limits = c(0,NA)) + ylab("-log10(p-value)")

linep4

# Part 3: Manhattan plot 

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
f <- file.path(resdir, "CL_0002057_IFNg_24h_eQTLres.gz")
transres <- readr::read_delim(f, delim = "\t", col_names=F)
colnames(transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
transres = transres %>% dplyr::group_by(snp) %>% dplyr::mutate(adj.pval = p.adjust(pvalue, method = "fdr"))

genes_in_lyz_cl = all_clusters[[meta_snp4]][["IC86_ICA_integrated"]]

transres2 = transres %>% filter(snp == meta_snp4)

# add gene metadata
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_lyz_cl 

transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))

transres2.tmp = transres2 %>% dplyr::group_by(chromosome) %>%
  dplyr::summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  dplyr::mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  dplyr::left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  dplyr::arrange(chromosome, gene_start) %>%
  dplyr::mutate(gene_startcum = gene_start + tot) 

axis.set = transres2.tmp %>% dplyr::group_by(chromosome) %>% 
  dplyr::summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)

axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp4 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: chr12_69344099_A_G, Condition: monocytes IFNg 24h, Module: IC86")

fp4 <- fp4 + geom_point(data = transres2.tmp %>% filter(highlight == T),
                      aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                      color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 1.9*10**(-40),
                                                  gene_name!="LYZ"),
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                  size = 10, point.padding = 0.1, force = 10)
# highlight cis eqtl
lyzdat = transres2.tmp %>% filter(gene_name == "LYZ")

fp4 <- fp4 + geom_point(data = lyzdat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                      size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = lyzdat,
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
                  color = "red", fontface = "bold", point.padding = 0.5, size = 10)
fp4

# Part 4: eigenvector expression

ic86_eigens = eigens %>% filter(variable == "IC86")

params = GRanges(c(12), IRanges(start = 69344099, end = 69344099))
tabix_header = Rsamtools::headerTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz")
tabix_header2 = tabix_header$header[229]
tabix_header2 = unlist(lapply(tabix_header2, function(x) strsplit(x, "\t")[[1]]))
tabix_list = Rsamtools::scanTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz", param = params)
lead_genotype = unlist(lapply(tabix_list, function(x) strsplit(x, "\t")[[1]]))
gen_df = data.frame(lead_genotype)
gen_df$genotype_id = tabix_header2
gen_df = gen_df[!gen_df$genotype_id %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),, drop = F]
row.names(gen_df) = NULL
colnames(gen_df) = c(meta_snp4, "genotype_id")
gen_df[[meta_snp4]] = as.character(gen_df[[meta_snp4]])
gen_df[[meta_snp4]] = unlist(lapply(gen_df[[meta_snp4]], function(x) strsplit(x, ":")[[1]][1]))
gen_df[[meta_snp4]] = unlist(lapply(gen_df[[meta_snp4]], function(x) ifelse(x == "0|0", "AA", ifelse(x == "1|1", "GG", "AG")) ))
gen_df[[meta_snp4]] = factor(gen_df[[meta_snp4]], levels = c("AA", "AG", "GG"))

ic86_eigens = merge(ic86_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp4 <- ggplot(ic86_eigens) + geom_boxplot(aes(x = chr12_69344099_A_G, y = value, 
                                                fill = qtl_group_label), show.legend = FALSE) + 
  facet_wrap(.~qtl_group_label, nrow = 1) + 
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) + 
  ylab("Eigengene expression") + labs(subtitle = "IC86 (ICA, integrated)")

genp4

# Merge all the plots together
part1 <- wrap_elements(plot = wrap_elements(leg2) / (cis_p4) / p41 / p44  + plot_layout(heights = c(0.5,5,5,5)))
part2 <- wrap_elements(linep4) 
part3 <- wrap_elements(fp4) / wrap_elements(genp4) 
patchfigg3 <- (part1 | part2) / part3 + plot_layout(heights = c(12, 16), widths = c(10, 5))
patchfigg3 <- patchfigg3 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 40, face = "bold"))
ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS4.png", patchfigg3, 
       width = 38, height = 30)

### Fig S5 - Effect sizes heatmap ####
# all metagroups
integrated_cl = res_filt %>% filter(approach == "integrated") %>%
  select(meta_id, cl_id, cl_id2, lead_pval, cl_method, approach, chr) %>% distinct()

heatmapdata = data.frame()

for(i in 1:nrow(integrated_cl)){
  approach = "integrated"
  approach2 = "full"
  m = integrated_cl[i,][["cl_method"]]
  meta = integrated_cl[i,][["meta_id"]]
  chr = as.numeric(gsub("chr", "", integrated_cl[i,][["chr"]]))
  pos = as.numeric(strsplit(meta, "_")[[1]][2])
  cl_id = integrated_cl[i,][["cl_id"]]
  param = GRanges(c(chr), IRanges(start = pos, end = pos))
  respath = file.path("/gpfs/hpc/home/liiskolb/transqtl_final/results", paste0(approach, "_coexpr"), m, "eQTLres")
  transfiles = list.files(path = respath, pattern = "*.gz$", full.names = T)
  print(cl_id)
  for (f in transfiles){
    transres <- data.frame(scanTabixDataFrame(f, param, col_names = FALSE))
    colnames(transres) <- c("chr", "start", "end", "snp", "cluster", "statistic", "pvalue", "FDR", "beta")
    transres = transres %>% filter(snp == meta & cluster == cl_id)
    if (nrow(transres)>0){
      qtl_group = gsub(paste("_", approach2, "_eigen_eQTLres.gz", sep = ""), "", gsub(paste(m, "_", sep = ""), "", basename(f)))
      transres$qtl_group = qtl_group
      transres$approach = approach
      transres$cl_method = m
      heatmapdata  = rbind(heatmapdata, transres)
    }
  }
}

heatmapdata = unique(heatmapdata)

heatmapdata$row_idx = paste(heatmapdata$cluster, heatmapdata$snp, sep = "_")

## for every snp keep the cluster with lowest p_value
keep_filt = heatmapdata %>% group_by(snp) %>% filter(pvalue == min(pvalue)) %>% data.frame
heatmapdata_filt = heatmapdata %>% filter(cluster %in% keep_filt$cluster)

# korruta kõik ühe rea beetad läbi absoluutväärtuselt kõige suurema beeta märgiga, 
# nii et kõige tugevam eQTL effect oleks alati positiivne ning teised siis selles suhtes

heatmapdata_scale = heatmapdata_filt %>% filter(row_idx %in% keep_filt$row_idx) %>% group_by(snp, cluster) %>% mutate(beta2 = beta*sign(beta[abs(beta) == max(abs(beta))]))
heatmapdata_scale = merge(heatmapdata_scale, unique(sample_metadata[,c("qtl_group", "cell_type", "condition", "cell_type_label", "condition_label")]), by = "qtl_group")
heatmapdata_scale$qtl_group_label = paste(heatmapdata_scale$cell_type_label, ifelse(heatmapdata_scale$condition_label == "naive", "", heatmapdata_scale$condition_label))
heatmapdata_scale$qtl_group_label = factor(heatmapdata_scale$qtl_group_label, levels = c("B cell ", "CD4+ T cell ", "CD8+ T cell ", "neutrophil ", "platelet ", "monocyte ", "monocyte LPS_2h", "monocyte LPS_24h", "monocyte IFNg_24h"))

colordf = data.frame(t(data.frame(colors2)))
colnames(colordf) = c("color")
colordf$qtl_group = rownames(colordf)

heatmapdata_scale = merge(heatmapdata_scale, colordf, by = "qtl_group")

annotcol = unique(heatmapdata_scale[,c("qtl_group_label", "color")])
row.names(annotcol) = annotcol$qtl_group_label
annotcol$color = as.character(annotcol$color)

fdr_ids = res_filt %>% filter(lead_FDR<0.1) %>% select(meta_id) %>% distinct()
bonf_ids = res_filt %>% filter(lead_pval<pval_thr) %>% select(meta_id) %>% distinct()

annotrow = data.frame("snp" = unique(heatmapdata_scale$snp)) 
annotrow[["FDR10%"]] = ifelse(annotrow$snp %in% fdr_ids$meta_id, "FDR", "FALSE")
annotrow$Bonferroni = ifelse(annotrow$snp %in% bonf_ids$meta_id, "Bonferroni", "FALSE")
row.names(annotrow) = annotrow$snp
annotrow$snp = NULL

annotcolors = list("qtl_group_label" = annotcol$color)
names(annotcolors[[1]]) = annotcol$qtl_group_label
annotcolors[["FDR10%"]] = c("FDR" = "red", "FALSE" = "grey")
annotcolors[["Bonferroni"]] = c("Bonferroni" = "blue", "FALSE" = "grey")

beta_wide4 = reshape2::dcast(heatmapdata_scale, snp~qtl_group_label, value.var="beta2")
row.names(beta_wide4) = beta_wide4$snp
beta_wide4$snp = NULL

sample_snps = list("chr3_56815721_T_C" = "ARHGEF3", "chr12_69344099_A_G" = "LYZ",
                   "chr9_20818520_A_G" = "IFNB1", "chr4_102325419_ACACT_A" = "SLC39A8",
                   "chr2_133579759_C_T" = "NCKAP5", "chr15_36752984_G_A" = "C15orf41"
)

row_labels = ifelse(row.names(beta_wide4) %in% names(sample_snps), sample_snps[row.names(beta_wide4)], "")
annotation_row = annotrow
breaksList = seq(-0.8, 0.8, by = 0.02)
hp = pheatmap(beta_wide4, fontsize_col = 15, 
              cutree_rows = 9, angle_col = "45",
              fontsize_row = 8, clustering_method = "complete", 
              clustering_distance_rows = "correlation", 
              cluster_cols = F, cluster_rows = T, 
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
              breaks = breaksList,
              annotation_col = annotcol[,c("qtl_group_label"), drop = F], 
              annotation_row = annotrow,
              annotation_colors = annotcolors,
              annotation_legend = FALSE, 
              annotation_names_col = FALSE,
              cellheight = 2,
              labels_row = row_labels)

# arrows added manually later
png("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS5.png", 
    height = 800, width = 800, res = 150)
hp
dev.off()

### Fig S6 - ARHGEF3 all modules ####

meta_snp2 = "chr3_56815721_T_C"
arh_res = res_filt %>% filter(meta_id == meta_snp2) # 3 clusters, only in platelet

# Part 1: Manhattan plots
chr = 3
position = 56815721
pheno_id2 = "ILMN_1781010" # ARHGEF3, ENSG00000163947
width = 100000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))

# trans plot 1
cl_id = "X6.WIERENGA_STAT5A_TARGETS_DN"
trans_resfile2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/PLIER/eQTLres/PLIER_CL_0000233_naive_separate_eigen_eQTLres.gz"
p2 <- transManhattan(trans_resfile2, cl_id, param, arh_res, meta_snp2, NULL)

p2 <- p2 + xlab("CHR 3 position (bp)") + theme(axis.title.x=element_text(),
                                               axis.text.x=element_text(),
                                               axis.ticks.x=element_line())
p2 <- p2 + labs(subtitle = "trans-eQTL for\nX6.WIERENGA_STAT5A_TARGETS_DN \n(PLIER, separate, platelets)")
# add rsid
p2[["layers"]][[4]] = NULL
p2 = p2 + ylim(0,10) + geom_text_repel(data = p2$data %>% filter(snp == meta_snp2),
                                       aes(x = start, y = -log10(pvalue)), label = "rs1354034",
                                       size = 10, nudge_y = 1, nudge_x = 10000)

# trans plot 2
cl_id = "IC68"
trans_resfile2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/integrated_coexpr/ICA/eQTLres/ICA_CL_0000233_naive_full_eigen_eQTLres.gz"
p3 <- transManhattan(trans_resfile2, cl_id, param, arh_res, meta_snp2, NULL)
# add x-axis labels and title
p3 <- p3 + xlab("CHR 3 position (bp)") + theme(axis.title.x=element_text(),
                                               axis.text.x=element_text(),
                                               axis.ticks.x=element_line())
p3 <- p3 + labs(subtitle = "trans-eQTL for IC68 \n(ICA, integrated)")
# add rsid
p3[["layers"]][[4]] = NULL
p3 = p3 + ylim(0,10) + geom_text_repel(data = p3$data %>% filter(snp == meta_snp2),
                                       aes(x = start, y = -log10(pvalue)), label = "rs1354034",
                                       size = 10, nudge_y = 1, nudge_x = 10000)

# trans plot 3
cl_id = "Cluster_12953"
trans_resfile2 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/funcExplorer/eQTLres/funcExplorer_CL_0000233_naive_separate_eigen_eQTLres.gz"
p4 <- transManhattan(trans_resfile2, cl_id, param, arh_res, meta_snp2, NULL)
# add x-axis labels and title
# add x-axis labels and title
p4 <- p4 + xlab("CHR 3 position (bp)") + theme(axis.title.x=element_text(),
                                               axis.text.x=element_text(),
                                               axis.ticks.x=element_line())
p4 <- p4 + labs(subtitle = "trans-eQTL for Cluster_12953\n(funcExplorer, separate, platelets)")
# add rsid
p4[["layers"]][[4]] = NULL
p4 = p4 + ylim(0,10) + geom_text_repel(data = p4$data %>% filter(snp == meta_snp2),
                          aes(x = start, y = -log10(pvalue)), label = "rs1354034",
                          size = 10, nudge_y = 1, nudge_x = 10000)

# Genotype plots
x6_eigens = eigens_sep %>% filter(variable == "X6.WIERENGA_STAT5A_TARGETS_DN", qtl_group == "CL_0000233_naive")
params = GRanges(c(3), IRanges(start = 56815721, end = 56815721))
tabix_header = Rsamtools::headerTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz")
tabix_header2 = tabix_header$header[229]
tabix_header2 = unlist(lapply(tabix_header2, function(x) strsplit(x, "\t")[[1]]))
tabix_list = Rsamtools::scanTabix(file = "/gpfs/hpc/home/liiskolb/transqtl_final/data/genotypes/Kolberg_2020.MAF005.norm.vcf.gz", param = params)
lead_genotype = unlist(lapply(tabix_list, function(x) strsplit(x, "\t")[[1]]))
gen_df = data.frame(lead_genotype)
gen_df$genotype_id = tabix_header2
gen_df = gen_df[!gen_df$genotype_id %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),, drop = F]
row.names(gen_df) = NULL
colnames(gen_df) = c(meta_snp2, "genotype_id")
gen_df[[meta_snp2]] = as.character(gen_df[[meta_snp2]])
gen_df[[meta_snp2]] = unlist(lapply(gen_df[[meta_snp2]], function(x) strsplit(x, ":")[[1]][1]))
gen_df[[meta_snp2]] = unlist(lapply(gen_df[[meta_snp2]], function(x) ifelse(x == "0|0", "TT", ifelse(x == "1|1", "CC", "TC")) ))
gen_df[[meta_snp2]] = factor(gen_df[[meta_snp2]], levels = c("TT", "TC", "CC"))
x6_eigens = merge(x6_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp2 <- ggplot(x6_eigens) + geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) +
  facet_wrap(.~qtl_group_label, nrow = 1) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  ylab("Eigengene expression") +
  labs(subtitle = "X6.WIERENGA_STAT5A_TARGETS_DN\n(PLIER, separate, platelets)") + xlab("rs1354034")
genp2

# IC68
ic68_eigens = eigens %>% filter(variable == "IC68")
ic68_eigens = merge(ic68_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp22 <- ggplot(ic68_eigens) + geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) + 
  facet_wrap(.~qtl_group_label, nrow = 1) + 
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) + 
  ylab("Eigengene expression") + labs(subtitle = "IC68 (ICA, integrated)") + xlab("rs1354034")

genp22

# Cluster_12953
cl12953_eigens = eigens_sep %>% filter(variable == "Cluster_12953", qtl_group == "CL_0000233_naive")
cl12953_eigens = merge(cl12953_eigens, gen_df, by.x = "genotype_id", by.y = "genotype_id")

genp23 <- ggplot(cl12953_eigens) + geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) + 
  facet_wrap(.~qtl_group_label, nrow = 1) + 
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) + 
  ylab("Eigengene expression") + labs(subtitle = "Cluster_12953\n(funcExplorer, separate, platelets)") +
  xlab("rs1354034")

genp23

# Manhattan plots

resdir = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/"
f <- file.path(resdir, "CL_0000233_naive_eQTLres.gz")
transres <- readr::read_delim(f, delim = "\t", col_names=F)
colnames(transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")
transres = transres %>% group_by(snp) %>% mutate(adj.pval = p.adjust(pvalue, method = "fdr"))
genes_in_arh_cl = all_clusters[[meta_snp2]][["X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive"]]
transres2 = transres %>% filter(snp == meta_snp2)
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_arh_cl
transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))
transres2.tmp = transres2 %>% group_by(chromosome) %>%
  summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  arrange(chromosome, gene_start) %>%
  mutate(gene_startcum = gene_start + tot)

axis.set = transres2.tmp %>% group_by(chromosome) %>%
  summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)
axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp2 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: rs1354034, Condition: platelet, Module: X6.WIERENGA_STAT5A_TARGETS_DN")

fp2 <- fp2 + geom_point(data = transres2.tmp %>% filter(highlight == T),
                        aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.001,
                                                  gene_name!="ARHGEF3"),
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                  size = 10, point.padding = 0.1, force = 10)

arhdat = transres2.tmp %>% filter(gene_name == "ARHGEF3")
fp2 <- fp2 + geom_point(data = arhdat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = arhdat,
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
                  color = "red", fontface = "bold", point.padding = 0.2, size = 10)
fp2

# IC68
genes_in_arh_cl = all_clusters[[meta_snp2]][["IC68_ICA_integrated"]]
transres2 = transres %>% filter(snp == meta_snp2)
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_arh_cl
transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))
transres2.tmp = transres2 %>% group_by(chromosome) %>%
  summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  arrange(chromosome, gene_start) %>%
  mutate(gene_startcum = gene_start + tot)

axis.set = transres2.tmp %>% group_by(chromosome) %>%
  summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)
axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp22 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: rs1354034, Condition: platelet, Module: IC68")

fp22 <- fp22 + geom_point(data = transres2.tmp %>% filter(highlight == T),
                        aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.001,
                                                  gene_name!="ARHGEF3"),
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                  size = 10, point.padding = 0.1, force = 10)

arhdat = transres2.tmp %>% filter(gene_name == "ARHGEF3")
fp22 <- fp22 + geom_point(data = arhdat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                        size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = arhdat,
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
                  color = "red", fontface = "bold", point.padding = 0.2, size = 10)
fp22

# Cluster_12953
genes_in_arh_cl = all_clusters[[meta_snp2]][["Cluster_12953_funcExplorer_separate_CL_0000233_naive"]]
transres2 = transres %>% filter(snp == meta_snp2)
transres2 = merge(transres2, genes_metadata, by.x = "gene_id", by.y = "gene_id")

# highlight trans relations
transres2$highlight = transres2$gene_id %in% genes_in_arh_cl
transres2$chromosome = factor(transres2$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))
transres2$highlight = factor(transres2$highlight, levels = c(TRUE, FALSE))
transres2.tmp = transres2 %>% group_by(chromosome) %>%
  summarise(chr_len = max(gene_start)) %>%
  # cumulative position
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(transres2, ., by = c("chromosome" = "chromosome")) %>%
  # cumulative position to each gene
  arrange(chromosome, gene_start) %>%
  mutate(gene_startcum = gene_start + tot)

axis.set = transres2.tmp %>% group_by(chromosome) %>%
  summarize(center = (max(gene_startcum) + min(gene_startcum)) / 2)
axis.set$chromosome = factor(axis.set$chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"))

fp23 <- ggplot(transres2.tmp, aes(x=gene_startcum, y=-log10(pvalue)*sign(beta))) +
  # Show all points
  geom_point(aes(color=adj.pval < 0.05), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  xlab("CHR") +
  theme_pubr(base_size = 30) +
  theme(legend.position="none") +
  labs(subtitle = "SNP: rs1354034, Condition: platelet, Module: Cluster_12953")

fp23 <- fp23 + geom_point(data = transres2.tmp %>% filter(highlight == T),
                          aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                          color = "darkblue", size = 3) +
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) +
  geom_text_repel(data = transres2.tmp %>% filter(highlight == T, adj.pval < 0.05,
                                                  gene_name!="ARHGEF3"),
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name),
                  size = 10, point.padding = 0.5, force = 10)

arhdat = transres2.tmp %>% filter(gene_name == "ARHGEF3")
fp23 <- fp23 + geom_point(data = arhdat, aes(x = gene_startcum, y = -log10(pvalue)*sign(beta)),
                          size = 8, stroke = 2, color = "red", shape = 21) + 
  geom_text_repel(data = arhdat,
                  aes(x = gene_startcum, y = -log10(pvalue)*sign(beta), label = gene_name), 
                  color = "red", fontface = "bold", nudge_y = 3, point.padding = 1, size = 10)
fp23

# arhgef expression across cell types
genp4 <- ggplot(arhgef_expr) +
  geom_boxplot(aes(x = chr3_56815721_T_C, y = value, fill = qtl_group_label), show.legend = FALSE) +
  facet_wrap(.~qtl_group_label, nrow = 1) +
  scale_fill_manual(values = unlist(colors)) + theme_pubr(base_size = 30) +
  ylab("Gene expression (log2)") + xlab("rs1354034") +
  labs(subtitle = "ARHGEF3")
genp4

# Merge all plots together

part1 = p2 / p3 / p4
part2 = (fp2 + ylim(-15,15)) / (fp22 + ylim(-15,15))/ (fp23 + ylim(-15,15))
part2 = part2 
part3 = wrap_plots(genp2 + theme(plot.margin = margin(1,10,1,1,"cm"))) | wrap_plots(genp22 + theme(strip.text = element_text(size = 18))) | wrap_plots(genp23 + theme(plot.margin = margin(1,10,1,1,"cm")))
part3 = part3 + plot_layout(widths = c(3, 14, 3))
patchfigg6 = part1 | part2
patchfigg6 = patchfigg6 + plot_layout(widths = c(5,14))
patchfigg7 = patchfigg6 / part3 / genp4
patchfigg7 = patchfigg7 + plot_layout(heights = c(16,4,4))
patchfigg7 = patchfigg7 + plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(size = 40, face = "bold"))

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS6.png", patchfigg7, 
       width = 42, height = 32, limitsize = F)

### Fig S7 - eQTLGen replication ####
# only trans associations
my_transres = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/full_trans_res.tsv", sep = "\t", header = T, stringsAsFactors = F)
eqtlgen = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/replication/eQTLgen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
meta_id = "chr3_56815721_T_C" # rs1354034

arhgef_my_trans = my_transres %>% filter(snp == meta_id, is_trans == TRUE, qtl_group == "CL_0000233_naive") %>% select(chr, start, end, snp, gene_id, statistic, pvalue, FDR, beta) %>% mutate(group = "Kolberg_2020")

f = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/CL_0000233_naive_eQTLres.gz"
param = GRanges(c(3), IRanges(start = 56815721, end = 56815721))
platelet_transres <- data.frame(scanTabixDataFrame(f, param, col_names = F))
colnames(platelet_transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")

arhgef_trans2 = eqtlgen %>% filter(SNP == "rs1354034") %>% mutate(group = "eQTLGen")
arhgef_trans2 = arhgef_trans2 %>% select(Gene, Zscore, group)
colnames(arhgef_trans2) = c("gene_id", "beta", "group")

arhgef_my_trans2 = platelet_transres %>% filter(snp == meta_id, gene_id %in% arhgef_trans2$gene_id) %>% 
  select(chr, start, end, snp, gene_id, statistic, pvalue, FDR, beta) %>% mutate(group = "Kolberg_2020")
arhgef_my_trans2 = arhgef_my_trans2 %>% select(gene_id, beta, group)

effectres2 = rbind(arhgef_my_trans2, arhgef_trans2)
effectres_wide2 = spread(effectres2[,c("gene_id", "group", "beta")], group, beta)

effectres_wide2$trans_in_Kolberg_2020 = effectres_wide2$gene_id %in% arhgef_my_trans$gene_id

##### color by cluster genes ####
aList = all_clusters[["chr3_56815721_T_C"]]
effectres_wide2$trans_in_cluster = unlist(lapply(effectres_wide2$gene_id, function(y) paste(names(aList)[sapply(1:3,function(x){ y %in% aList[[x]]})], collapse = "; ")))
effectres_wide2$trans_in_cluster[effectres_wide2$trans_in_cluster == ""] = "Not in any module"
effectres_wide2$trans_in_cluster[effectres_wide2$trans_in_cluster %like% "Cluster_"] = "Cluster_12953"
effectres_wide2$trans_in_cluster[effectres_wide2$trans_in_cluster == "X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive"] = "X6.WIERENGA_STAT5A_TARGETS_DN"
effectres_wide2$trans_in_cluster[effectres_wide2$trans_in_cluster == "IC68_ICA_integrated"] = "IC68"
effectres_wide2$trans_in_cluster[effectres_wide2$trans_in_cluster == "IC68_ICA_integrated; X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive"] = "IC68 and X6.WIERENGA_STAT5A_TARGETS_DN"
effectres_wide2$trans_in_cluster = factor(effectres_wide2$trans_in_cluster)
effectres_wide2$trans_in_cluster = factor(effectres_wide2$trans_in_cluster, levels = c("Not in any module", "X6.WIERENGA_STAT5A_TARGETS_DN", "IC68", "IC68 and X6.WIERENGA_STAT5A_TARGETS_DN", "Cluster_12953"))
#levels(effectres_wide2$trans_in_cluster) = c("Not in any module", "X6.WIERENGA_STAT5A_TARGETS_DN", "IC68", "IC68 and X6.WIERENGA_STAT5A_TARGETS_DN", "Cluster_12953")
effectres_wide2 = effectres_wide2 %>% arrange(trans_in_cluster)

test = ggscatter(data = effectres_wide2[complete.cases(effectres_wide2),],
          x = "eQTLGen", y = "Kolberg_2020",
          add = "none", alpha = 0.8, fill = "darkgrey",
          cor.coef = FALSE,
          color = "trans_in_cluster",
          group = "trans_in_cluster",
          palette = unlist(list("Not in any module" = "darkgrey", "Cluster_12953" = "darkgreen", "X6.WIERENGA_STAT5A_TARGETS_DN" = "blue", "IC68" = "chocolate1", "IC68 and X6.WIERENGA_STAT5A_TARGETS_DN" = "maroon")),
          size = 2,
          fullrange = TRUE,
          cor.coeff.args = list(method = "pearson", size = 5,  label.sep = "\n"),
          add.params = list(color = "blue", size = 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylim(-1,1) +
  xlim(-60,60) +
  theme_pubr(base_size = 20, legend = "right") + labs(subtitle = "rs1354034\neQTLGen FDR < 0.05 trans-eQTLs ")
test
test + stat_cor(data = effectres_wide2[complete.cases(effectres_wide2),] %>% filter(trans_in_cluster %like% "X6.WIERENGA_STAT5A_TARGETS_DN"), method = "pearson", size = 5,  label.sep = "\n")
test + stat_cor(data = effectres_wide2[complete.cases(effectres_wide2),] %>% filter(trans_in_cluster %like% "IC68"), method = "pearson", size = 5,  label.sep = "\n")


#####
effectres_wide2 = effectres_wide2 %>% arrange(trans_in_Kolberg_2020)
arhgef2 = ggscatter(data = effectres_wide2[complete.cases(effectres_wide2),],
                    x = "eQTLGen", y = "Kolberg_2020",
                    add = "none", alpha = 0.8, fill = "darkgrey",
                    cor.coef = FALSE,
                    color = "trans_in_Kolberg_2020",
                    group = "trans_in_Kolberg_2020",
                    palette = c("darkgrey", "darkgreen"),
                    size = 2,
                    fullrange = TRUE,
                    cor.coeff.args = list(method = "pearson", size = 5,  label.sep = "\n"),
                    add.params = list(color = "blue", size = 1)) +
 geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 ylim(-1,1) +
 xlim(-60,60) +
  xlab("eQTLGen (Z-score)") + 
  ylab("Kolberg_2020 (beta)")
theme_pubr(base_size = 20, legend = "bottom") + labs(subtitle = "rs1354034\neQTLGen FDR < 0.05 trans-eQTLs ")

arhgef2 = arhgef2 + stat_cor(data = effectres_wide2[complete.cases(effectres_wide2),] %>% filter(trans_in_Kolberg_2020 == T), method = "pearson", color = "darkgreen", size = 5,  label.sep = "\n")
ggsave(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS7.png", arhgef2, height = 4, width = 5)

### Fig S8 - SLC all conditions ####
meta_snp3 = "chr4_102325419_ACACT_A"

slc_res = res_filt %>% filter(meta_id == meta_snp3) # 1 cluster, only in LPS_24h
chr = 4
position = 102325419
pheno_id1 = "ILMN_1695316" # SLC39A8, ENSG00000138821
pheno_id2 = "ILMN_2233539" # SLC39A8
missense = "chr4_102267552_C_T"
width = 200000

phenos = c("ILMN_1695316", "ILMN_2233539")

param <- GRanges(c(chr), IRanges(position - width, end = position + width))

## trans plot
cl_id = "Cluster_10413"
trans_resfile3 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/funcExplorer/eQTLres/funcExplorer_CL_0002057_LPS_24h_separate_eigen_eQTLres.gz"
p <- transManhattan(trans_resfile3, cl_id, param, slc_res, meta_snp3, NULL)
# add x-axis labels and title
p <- p + xlab("CHR 4 position (bp)") +
  labs(subtitle = "trans-eQTL for Cluster_10413\n(funcExplorer, separate, monocytes LPS 24h)") + 
  theme(axis.title.x=element_text(), axis.text.x=element_text(), axis.ticks.x=element_line())

p[["layers"]][[4]] = NULL

p <- p + geom_text_repel(data = p$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 1) + ylim(0,10)
# add missense
p_miss = p + geom_point(data = p$data %>% filter(snp == missense), aes(x = start, y = -log10(pvalue)),
             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = p$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), nudge_y = 3, label = "rs13107325", size = 10)

# LPS 24h
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_LPS_24h.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE)[[1]])
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% dplyr::filter(phenotype_id %in% phenos)

cis_lead_snp3 = cisres %>% filter(is_top == 1) 

cis_p3 <- ggplot(cisres %>% filter(phenotype_id == pheno_id1)) + 
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, 
             color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p3 <- cis_p3 + labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)\n", "(monocytes LPS 24h)"))
#cis_p3 <- cis_p3 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id1) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, linetype="dotted", color = "black")

cis_p4 <- ggplot(cisres %>% filter(phenotype_id == pheno_id2)) +
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p4 <- cis_p4 + labs(subtitle = paste("cis-eQTL for",  pheno_id2 ,"(SLC39A8)\n" , "(monocytes LPS 24h)"))
#cis_p4 <- cis_p4 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id2) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, 
#                              linetype="dotted", color = "black")

# add missense
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", nudge_x = 10000, nudge_y = 3, size = 10)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", 
                  nudge_x = -10000, nudge_y = 15, size = 10)

# highlight trans lead 
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 1)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 10)

# highlight trans credible set
trans_cs = p_miss$data %>% filter(highlight == T)

cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

# highlight cis leads
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), size = 10, nudge_y = 2) + ylim(0,18)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 5) + ylim(0,85)

p_miss2 = p_miss + geom_point(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                              aes(x = start, y = -log10(pvalue)),
                              size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 1) + ylim(0,10)

# cis credible sets
slc_cis_cs = ciscs %>% filter(phenotype_id == pheno_id1, qtl_group == "CL_0002057_LPS_24h")
slc_cis_cs$cs_id = paste(slc_cis_cs$cs_id, pheno_id1)
slc_cis_cs2 = ciscs %>% filter(phenotype_id == pheno_id2, qtl_group == "CL_0002057_LPS_24h")
slc_cis_cs2$cs_id = paste(slc_cis_cs2$cs_id, pheno_id2)
slc_cis_cs_full = rbind(slc_cis_cs, slc_cis_cs2)
slc_cis_cs_full$cs_id = factor(slc_cis_cs_full$cs_id)
levels(slc_cis_cs_full$cs_id) = c("CS1 ILMN_2233539", "CS2 ILMN_2233539", "CS1 ILMN_1695316")

colorvalues = c("CS1 ILMN_2233539" = "maroon" , "CS2 ILMN_2233539" = "maroon1",  "CS1 ILMN_1695316" = "chocolate1")

# highlight cis credible set
cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id1), 
                              aes(x = start, y = -log10(pvalue), color = cs_id), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)

cis_p4 <- cis_p4 + geom_point(data = cis_p4$data %>% filter(phenotype_id == pheno_id2) %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id2), 
                              aes(x = start, y = -log10(pvalue), color = as.factor(cs_id)), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)


# merge
suppl1 = cis_p3/cis_p4/p_miss2

# naive
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_naive.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE)[[1]])
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% dplyr::filter(phenotype_id %in% phenos)

cis_lead_snp3 = cisres %>% filter(is_top == 1) 

cis_p3 <- ggplot(cisres %>% filter(phenotype_id == pheno_id1)) + 
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, 
             color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p3 <- cis_p3 + labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)\n", "(monocytes naive)"))
#cis_p3 <- cis_p3 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id1) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, linetype="dotted", color = "black")

cis_p4 <- ggplot(cisres %>% filter(phenotype_id == pheno_id2)) +
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p4 <- cis_p4 + labs(subtitle = paste("cis-eQTL for",  pheno_id2 ,"(SLC39A8)\n" , "(monocytes naive)"))
#cis_p4 <- cis_p4 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id2) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, 
#                              linetype="dotted", color = "black")

# add missense
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", nudge_x = 10000, nudge_y = 3, size = 10)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", 
                  nudge_x = -10000, nudge_y = 15, size = 10)

# highlight trans lead 
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 1)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 10)

# highlight trans credible set
trans_cs = p_miss$data %>% filter(highlight == T)

cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

# highlight cis leads
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), size = 10, nudge_y = 2) + ylim(0,18)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 8) + ylim(0,85)

p_miss2 = p_miss + geom_point(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                              aes(x = start, y = -log10(pvalue)),
                              size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 1) + ylim(0,10)

# cis credible sets
slc_cis_cs = ciscs %>% filter(phenotype_id == pheno_id1, qtl_group == "CL_0002057_naive")
slc_cis_cs$cs_id = paste(slc_cis_cs$cs_id, pheno_id1)
slc_cis_cs2 = ciscs %>% filter(phenotype_id == pheno_id2, qtl_group == "CL_0002057_naive")
slc_cis_cs2$cs_id = paste(slc_cis_cs2$cs_id, pheno_id2)
slc_cis_cs_full = rbind(slc_cis_cs, slc_cis_cs2)
slc_cis_cs_full$cs_id = factor(slc_cis_cs_full$cs_id)
levels(slc_cis_cs_full$cs_id) = c("CS1 ILMN_2233539", "CS2 ILMN_2233539", "CS1 ILMN_1695316")

colorvalues = c("CS1 ILMN_2233539" = "maroon" , "CS2 ILMN_2233539" = "maroon1",  "CS1 ILMN_1695316" = "chocolate1")

# highlight cis credible set
cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id1), 
                              aes(x = start, y = -log10(pvalue), color = cs_id), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)

cis_p4 <- cis_p4 + geom_point(data = cis_p4$data %>% filter(phenotype_id == pheno_id2) %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id2), 
                              aes(x = start, y = -log10(pvalue), color = as.factor(cs_id)), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)


suppl2 <- cis_p3/cis_p4/p_miss2

# LPS 2h
cis_resfile = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/Kolberg_2020_QC2/Kolberg_2020_CL_0002057_LPS_2h.nominal.sorted.txt.gz"
cisres <- data.frame(scanTabixDataFrame(cis_resfile, param, col_names = FALSE)[[1]])
colnames(cisres) <- c("phenotype_id",  "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "nr_variants_tested", "pheno_variant_dist", "snp", "chr", "start", "snp_end", "pvalue", "beta", "is_top")
cisres <- cisres %>% dplyr::filter(phenotype_id %in% phenos)

cis_lead_snp3 = cisres %>% filter(is_top == 1) 

cis_p3 <- ggplot(cisres %>% filter(phenotype_id == pheno_id1)) + 
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, 
             color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p3 <- cis_p3 + labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)\n", "(monocytes LPS 2h)"))
#cis_p3 <- cis_p3 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id1) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, linetype="dotted", color = "black")

cis_p4 <- ggplot(cisres %>% filter(phenotype_id == pheno_id2)) +
  geom_point(aes(x = start, y = -log10(pvalue)), size = 3, alpha = 0.6, color = "#8c8c8c") + 
  theme_pubr(base_size = 30) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
cis_p4 <- cis_p4 + labs(subtitle = paste("cis-eQTL for",  pheno_id2 ,"(SLC39A8)\n" , "(monocytes LPS 2h)"))
#cis_p4 <- cis_p4 + geom_vline(data = cisres %>% filter(phenotype_id == pheno_id2) %>% head(1), 
#                              aes(xintercept = phenotype_start), size = 1, 
#                              linetype="dotted", color = "black")

# add missense
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", nudge_x = 10000, nudge_y = 3, size = 10)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == missense), 
                             aes(x = start, y = -log10(pvalue)),
                             shape = 18, size = 8, color = "blue") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == missense), 
                  aes(x = start, y = -log10(pvalue)), label = "rs13107325", 
                  nudge_x = -10000, nudge_y = 15, size = 10)

# highlight trans lead 
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 1)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp == meta_snp3), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp == meta_snp3), aes(x = start, y = -log10(pvalue)), 
                  label = "rs75562818", size = 10, nudge_y = 10)

# highlight trans credible set
trans_cs = p_miss$data %>% filter(highlight == T)

cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% trans_cs$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 3, color = "#00b159") 

# highlight cis leads
cis_p3 = cis_p3 + geom_point(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p3$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), size = 10, nudge_y = 2) + ylim(0,18)

cis_p4 = cis_p4 + geom_point(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                             aes(x = start, y = -log10(pvalue)),
                             size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = cis_p4$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 8) + ylim(0,85)

p_miss2 = p_miss + geom_point(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                              aes(x = start, y = -log10(pvalue)),
                              size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = p_miss$data %>% filter(snp %in% cis_lead_snp3$snp), 
                  aes(x = start, y = -log10(pvalue), label = snp), 
                  size = 10, nudge_y = 1) + ylim(0,10)

# cis credible sets
slc_cis_cs = ciscs %>% filter(phenotype_id == pheno_id1, qtl_group == "CL_0002057_LPS_2h")
slc_cis_cs$cs_id = paste(slc_cis_cs$cs_id, pheno_id1)
slc_cis_cs2 = ciscs %>% filter(phenotype_id == pheno_id2, qtl_group == "CL_0002057_LPS_2h")
slc_cis_cs2$cs_id = paste(slc_cis_cs2$cs_id, pheno_id2)
slc_cis_cs_full = rbind(slc_cis_cs, slc_cis_cs2)
slc_cis_cs_full$cs_id = factor(slc_cis_cs_full$cs_id)
levels(slc_cis_cs_full$cs_id) = c("CS1 ILMN_2233539", "CS2 ILMN_2233539", "CS1 ILMN_1695316")

colorvalues = c("CS1 ILMN_2233539" = "maroon" , "CS2 ILMN_2233539" = "maroon1",  "CS1 ILMN_1695316" = "chocolate1")

# highlight cis credible set
cis_p3 <- cis_p3 + geom_point(data = cis_p3$data %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id1), 
                              aes(x = start, y = -log10(pvalue), color = cs_id), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)

cis_p4 <- cis_p4 + geom_point(data = cis_p4$data %>% filter(phenotype_id == pheno_id2) %>% left_join(slc_cis_cs_full, by = c("snp" = "variant_id")) %>% filter(snp %in% slc_cis_cs_full$variant_id, cs_id %like% pheno_id2), 
                              aes(x = start, y = -log10(pvalue), color = as.factor(cs_id)), 
                              size = 3) + 
  scale_color_manual(name = "", values = colorvalues, drop = FALSE)

suppl3 <- cis_p3/cis_p4/p_miss2


# Save as supplementary figure

supplall = (suppl2)| (suppl3) | (suppl1)
supplall = supplall + plot_annotation(tag_levels = c("A"))

ggsave(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS8.png", supplall, 
       height = 25, width = 40)


### Fig S9 - SLC Kim-Hellmuth all conditions ####

chr = 4
position = 102325419
pheno_id1 = "ILMN_1695316" # SLC39A8, ENSG00000138821
#pheno_id2 = "ILMN_2233539" # SLC39A8
width = 200000
param <- GRanges(c(chr), IRanges(position - width, end = position + width))

slc_sub = res_filt %>% filter(meta_id == "chr4_102325419_ACACT_A") %>% select(snp_id, rs_id) %>% distinct()
slc_sub$rs_id[slc_sub$snp_id == "chr4_102325419_ACACT_A"] = "rs75562818"
slc_sub$rs_id[slc_sub$snp_id == "chr4_102340060_T_A"] = "rs62329461"

kimnaive = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/Kim-Hellmuth/nominal_ctrl.txt", sep = " ", stringsAsFactors = F)
kim6h = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/Kim-Hellmuth/nominal_lps6h.txt", sep = " ", stringsAsFactors = F)
kim90 = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/Kim-Hellmuth/nominal_lps90.txt", sep = " ", stringsAsFactors = F)

colnames(kimnaive) = c("probe_id", "snp", "distance", "pval", "beta")
colnames(kim6h) = c("probe_id", "snp", "distance", "pval", "beta")
colnames(kim90) = c("probe_id", "snp", "distance", "pval", "beta")

cispnaive = ggplot(kimnaive %>% filter(probe_id == pheno_id1, distance <= 64014+width & distance >= -width + 64014)) + 
  geom_point(aes(x = distance, y = -log10(pval)), size = 3, alpha = 0.6, color = "#8c8c8c") +
  theme_pubr(base_size = 20) + 
  ylim(0,5) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) + 
  labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)", "\nKim-Hellmuth et al. (monocytes naive)"))

kimnaive[kimnaive$snp == "rs145021182","snp"] = "rs75562818"
kim6h[kim6h$snp == "rs145021182","snp"] = "rs75562818"
kim90[kim90$snp == "rs145021182","snp"] = "rs75562818"

kimnaive = merge(kimnaive, slc_sub, by.x = "snp", by.y = "rs_id", all.x = T)

kimnaive$show_label = kimnaive$snp_id %in% c(meta_snp3)
## add trans credible set 
cispnaive <- cispnaive + geom_point(data = kimnaive %>% filter(probe_id == pheno_id1, snp %in% slc_sub$rs_id,
                                                               distance <= 64014+width & distance >= -width + 64014), 
                                    aes(x = distance, y = -log10(pval)), size = 3, colour = "#00b159") 
# highlight lead SNP
cispnaive = cispnaive + geom_point(data = kimnaive %>% filter(probe_id == pheno_id1, 
                                  distance <= 64014+width & distance >= -width + 64014, snp_id == meta_snp3), 
                                   aes(x = distance, y = -log10(pval)), size = 6, shape = 21, fill = "#d11141") +
  geom_text_repel(data = kimnaive %>% filter(probe_id == pheno_id1, 
                                             show_label == T, distance <= 64014+width & distance >= -width + 64014), 
                  aes(x = distance, y = -log10(pval), label = snp), size = 10, nudge_x = 1000, nudge_y = 1)  


cisp6h = ggplot(kim6h %>% filter(probe_id == pheno_id1, distance <= 64014+width & distance >= -width + 64014)) + 
  geom_point(aes(x = distance, y = -log10(pval)), size = 3, alpha = 0.6, color = "#8c8c8c") +
  theme_pubr(base_size = 20) + 
  ylim(0,5) +
  labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)", "\nKim-Hellmuth et al. (monocytes LPS 6h)")) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

kim6h = merge(kim6h, slc_sub, by.x = "snp", by.y = "rs_id", all.x = T)
kim6h$show_label = kim6h$snp_id %in% c(meta_snp3)
## add trans credible set juurde
cisp6h <- cisp6h + geom_point(data = kim6h %>% filter(probe_id == pheno_id1, snp %in% slc_sub$rs_id,
                                                      distance <= 64014+width & distance >= -width + 64014), 
                              aes(x = distance, y = -log10(pval)), size = 3, colour = "#00b159") 
# highlight lead snp
cisp6h = cisp6h + geom_point(data = kim6h %>% filter(probe_id == pheno_id1, 
                                      distance <= 64014+width & distance >= -width + 64014, snp_id == meta_snp3), 
                             aes(x = distance, y = -log10(pval)), size = 6, shape = 21, fill = "#d11141") +
  geom_text_repel(data = kim6h %>% filter(probe_id == pheno_id1, snp_id == meta_snp3, distance <= 64014+width & distance >= -width + 64014), 
                  aes(x = distance, y = -log10(pval), label = snp), size = 10, nudge_x = 1000, nudge_y = 1) 

# LPS 90min
cisp90 = ggplot(kim90 %>% filter(probe_id == pheno_id1, distance <= 64014+width & distance >= -width + 64014)) + 
  geom_point(aes(x = distance, y = -log10(pval)), size = 3, alpha = 0.6, color = "#8c8c8c") +
  theme_pubr(base_size = 20) + 
  labs(subtitle = paste("cis-eQTL for",  pheno_id1, "(SLC39A8)", "\nKim-Hellmuth et al. (monocytes LPS 90min)")) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

kim90 = merge(kim90, slc_sub, by.x = "snp", by.y = "rs_id", all.x = T)
kim90$show_label = kim90$snp_id %in% c(meta_snp3)
## add trans credible set juurde
cisp90 <- cisp90 + geom_point(data = kim90 %>% filter(probe_id == pheno_id1, snp %in% slc_sub$rs_id,
                                                      distance <= 64014+width & distance >= -width + 64014), 
                              aes(x = distance, y = -log10(pval)), size = 3, colour = "#00b159") 
# highlight lead snp
cisp90 = cisp90 + geom_point(data = kim90 %>% filter(probe_id == pheno_id1, 
                                                     distance <= 64014+width & distance >= -width + 64014, snp_id == meta_snp3), 
                             aes(x = distance, y = -log10(pval)), size = 6, shape = 21, fill = "#d11141") +
  geom_text_repel(data = kim90 %>% filter(probe_id == pheno_id1, snp_id == meta_snp3, distance <= 64014+width & distance >= -width + 64014), 
                  aes(x = distance, y = -log10(pval), label = snp), size = 10, nudge_x = 1000, nudge_y = 1) 

cl_id = "Cluster_10413"
trans_resfile3 = "/gpfs/hpc/home/liiskolb/transqtl_final/results/separate_coexpr/funcExplorer/eQTLres/funcExplorer_CL_0002057_LPS_24h_separate_eigen_eQTLres.gz"
ps <- transManhattan(trans_resfile3, cl_id, param, slc_res, meta_snp3, NULL)
# add x-axis labels and title
ps <- ps + xlab("CHR 4 position (bp)") +
  labs(subtitle = "trans-eQTL for Cluster_10413\n(funcExplorer, separate, monocytes LPS 24h)") + 
  theme(axis.title.x=element_text(), axis.text.x=element_text(), axis.ticks.x=element_line())

ps[["layers"]][[4]] = NULL
ps[["layers"]][[3]] = NULL

# add rs ids
rs_data = ps$data %>% filter(snp %in% c(meta_snp3))
rs_data = merge(rs_data, slc_sub, by.x = "snp", by.y = "snp_id")
ps <- ps + ylim(0,12) + geom_point(data = rs_data, 
                                   aes(x = start, y = -log10(pvalue)), 
                                   size = 6, shape = 21, fill = "#d11141") + 
  geom_text_repel(data = rs_data, 
                  aes(x = start, y = -log10(pvalue), label = rs_id),
                  size = 10, force = 15, nudge_y = 2) + theme_pubr(base_size = 20)

mpatch1 = wrap_elements(cispnaive/(ps))
mpatch2 = wrap_elements((cisp90 + ylim(0,8)) /ps) 
mpatch3 = wrap_elements(cisp6h/ps) 

mpatch = mpatch1 + mpatch2 + mpatch3 + plot_annotation(tag_levels = "A") + plot_layout(ncol = 2) & theme(plot.tag = element_text(size = 24, face = "bold"))

ggsave("/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS9.png", mpatch, width = 17, height = 12)

### Fig S10 - SLC replication in FF2018 ####
meta_id = "chr4_102325419_ACACT_A"
cl_id = "Cluster_10413"

# only trans associations
my_transres = read.table("/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/full_trans_res.tsv", sep = "\t", header = T, stringsAsFactors = F)

# effect sizes
width = 0
param = GRanges(c(4), IRanges(start = 102325419-width, end = 102325419+width))

f = "/gpfs/hpc/home/liiskolb/transqtl_final/replication/Fairfax_2018/eQTLs/CL_0002057_LPS_24h_eQTLres.gz"
ff_transres <- data.frame(scanTabixDataFrame(f, param, col_names = F))
colnames(ff_transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")

slc_my_trans = my_transres %>% filter(snp == meta_id, is_trans == TRUE) %>% select(chr, start, end, snp, gene_id, statistic, pvalue, FDR, beta) %>% mutate(group = "Kolberg_2020")
slc_ff_trans = ff_transres %>% filter(gene_id %in% slc_my_trans$gene_id, snp == meta_id) %>% mutate(group = "Fairfax_2018")

effectres = rbind(slc_my_trans, slc_ff_trans)
effectres_wide = spread(effectres[,c("gene_id", "group", "beta")], group, beta)

## naive

f = "/gpfs/hpc/home/liiskolb/transqtl_final/replication/Fairfax_2018/eQTLs/CL_0002057_naive_eQTLres.gz"
ff_transres2 <- data.frame(scanTabixDataFrame(f, param, col_names = F))
colnames(ff_transres2) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")

nf = "/gpfs/hpc/home/liiskolb/transqtl_final/crediblesets/full-meta-trans/eQTLs/CL_0002057_naive_eQTLres.gz"
my_naive_transres = data.frame(scanTabixDataFrame(nf, param, col_names = F))
colnames(my_naive_transres) <- c("chr", "start", "end", "snp", "gene_id", "statistic", "pvalue", "FDR", "beta")

slc_my_trans2 = my_naive_transres %>% filter(snp == meta_id, gene_id %in% slc_my_trans$gene_id) %>% select(chr, start, end, snp, gene_id, statistic, pvalue, FDR, beta) %>% mutate(group = "Kolberg_2020")
slc_ff_trans2 = ff_transres2 %>% filter(gene_id %in% slc_my_trans$gene_id, snp == meta_id) %>% mutate(group = "Fairfax_2018")

effectres2 = rbind(slc_my_trans2, slc_ff_trans2)
effectres_wide2 = spread(effectres2[,c("gene_id", "group", "beta")], group, beta)

lps24 = ggplot(effectres_wide) + 
  geom_point(aes(x = Fairfax_2018, y = Kolberg_2020), size = 3, alpha = 0.6) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  xlim(-1, 1) + 
  ylim(-1, 1) +
  theme_pubr(base_size = 20) + labs(subtitle = "monocytes LPS 24h")

naive = ggplot(effectres_wide2) + 
  geom_point(aes(x = Fairfax_2018, y = Kolberg_2020), size = 3, alpha = 0.6) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  xlim(-1, 1) + 
  ylim(-1, 1) +
  theme_pubr(base_size = 20) + labs(subtitle = "monocytes naive")

pw1 = (naive | lps24) + plot_annotation(title = "rs75562818 trans effect sizes") & 
  theme(plot.title = element_text(size = 20))
ggsave(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS10.png", pw1, width = 8)

### Table S1 - replication overlaps ####
# LYZ and FF 2012
library("gdata")
ff2012_lyz = read.xls("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3437404/bin/NIHMS40993-supplement-6.xls")
ff2012_lyz = ff2012_lyz %>% filter(P.rs10784774 < 10**(-11), CHR!=12)
ff2012_ensg = gconvert(ff2012_lyz$Gene)
ff2012_ensg = unique(ff2012_ensg$target)
length(intersect(ff2012_ensg, unlist(all_clusters[["chr12_69344099_A_G"]])))

# IFNB1 and FF 2014
library(xlsx)
tmp <- tempfile(fileext = ".xlsx")
download.file(url = "http://science.sciencemag.org/highwire/filestream/595699/field_highwire_adjunct_files/3/1246949stableS4.xlsx", destfile = tmp, mode="wb")

ff2014_ifnb = read.xlsx(tmp, sheetIndex = 2)
ff2014_ifnb = ff2014_ifnb %>% filter(LPS24.FDR <0.05, LPS24 == 1)
ff2014_ifnb_ensg = gconvert(ff2014_ifnb$Gene)
ff2014_ifnb_ensg = unique(ff2014_ifnb_ensg$target)
length(intersect(ff2014_ifnb_ensg, unlist(all_clusters[["chr9_20818520_A_G"]])))

# IFNB1 and Quach
tmp <- tempfile(fileext = ".xlsx")
download.file(url = "http://ars.els-cdn.com/content/image/1-s2.0-S009286741631306X-mmc3.xlsx", destfile = tmp, mode="wb")
quach_ifnb = read.xlsx(tmp, sheetIndex = 1, startRow = 16)
quach_ifnb = quach_ifnb %>% filter(SNPa == "rs12553564", Condition == "LPS") %>% select(trans.regulated.genes..up.e, trans.regulated.genes..down.e)
quach_ifnb_genes = unique(c(as.character(quach_ifnb$trans.regulated.genes..up.e), strsplit(as.character(quach_ifnb$trans.regulated.genes..down.e), " // ")[[1]]))
quach_ifnb_ensg = gconvert(quach_ifnb_genes)
quach_ifnb_ensg = unique(quach_ifnb_ensg$target)
length(intersect(quach_ifnb_ensg, unlist(all_clusters[["chr9_20818520_A_G"]])))

#### ARHGEF3 cluster overlap significance ####
ic68 = all_clusters[["chr3_56815721_T_C"]][["IC68_ICA_integrated"]]
stat5 = all_clusters[["chr3_56815721_T_C"]][["X6.WIERENGA_STAT5A_TARGETS_DN_PLIER_separate_CL_0000233_naive"]]

inter = length(intersect(ic68, stat5))
fisher.test(matrix(c(74, 1000, 844, 16465), nrow = 2), alternative = "greater")

#### Gene level overlap significance plot
metares_filt = read.table(file = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/suppl_data/filt_crediblesets.tsv", sep = "\t", stringsAsFactors = F)
p = ggplot(metares_filt %>% select(cl_id2, cl_method, approach, meta_id, fisher_pval, cl_size, nr_gene_trans, nr_trans_in_cl) %>% distinct()) + 
  geom_histogram(aes(x=fisher_pval), bins = 20) + theme_bw() + facet_grid(cl_method~approach, scales = "free_y") + 
  xlab("P-value (one-sided Fisher’s exact test)") + ylab("Frequency")
ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/figures/FigS_fisher.png", p, height = 6)


### Figure 1 -- figure supplement 6 ####

# generate simulation vectors into a matrix
# samples in columns

# samples in columns
mat = matrix(
  c(rep(1,10), rep(0,10), rep(-1,10),
    rep(1,10), rep(sqrt(2),10), rep(1,10), 
    rep(1,10), rep(-sqrt(2),10), rep(1,10)),
  ncol = 30,
  byrow = TRUE
)

# duplicate data
mat_longer = mat[c(rep(1, 150), rep(2, 100), rep(3, 50)),]

# add noise
for(i in 1:nrow(mat_longer)){
  #noise = rnorm(30, mean = 0, sd = 0.000001)
  noise = rnorm(30, mean = 0, sd = 0.1)
  mat_longer[i,] = mat_longer[i,] + noise
}
# genes are in rows and samples in columns
mat_scaled = t(scale(t(mat_longer), center = F)) # unit variance for genes

# PCA

eigens = prcomp(mat_scaled, center = F, scale = F)
eigengene1 = eigens$rotation[,1]
eigengene2 = eigens$rotation[,2]
eigengene3 = eigens$rotation[,3]

# eigenvectors per module
eigengene11 = prcomp(mat_scaled[1:150,], center = F, scale = F)$rotation[,1]
eigengene22 = prcomp(mat_scaled[151:250,], center = F, scale = F)$rotation[,1]
eigengene33 = prcomp(mat_scaled[251:300,], center = F, scale = F)$rotation[,1]

# heatmap of expression matrix
row.names(mat_scaled) = paste("gene_", seq(1,nrow(mat_scaled)))
row.names(mat_scaled)[c(2:149,151:249,251:299)] = ""
colnames(mat_scaled) = paste("sample_", seq(1, ncol(mat_scaled)))
pheatmap(mat_scaled, gaps_row = c(150,250), cluster_rows = F, cluster_cols = F, angle_col = "90", fontsize = 12)

# eigengene profiles
e1 = ggplot(data.frame("Loadings" = eigengene1, "Sample_id" = seq(1, length(eigengene1))), 
       aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) + 
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene1.png", e1, width = 3, height = 2)


e2 = ggplot(data.frame("Loadings" = eigengene2, "Sample_id" = seq(1, length(eigengene2))), 
            aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) +
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene2.png", e2, width = 3, height = 2)

e3 = ggplot(data.frame("Loadings" = eigengene3, "Sample_id" = seq(1, length(eigengene3))), 
            aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) +
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene3.png", e3, width = 3, height = 2)

e11 = ggplot(data.frame("Loadings" = eigengene11, "Sample_id" = seq(1, length(eigengene11))), 
            aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) + 
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene11.png", e11, width = 3, height = 2)


e22 = ggplot(data.frame("Loadings" = eigengene22, "Sample_id" = seq(1, length(eigengene22))), 
            aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) +
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene22.png", e22, width = 3, height = 2)

e33 = ggplot(data.frame("Loadings" = eigengene33, "Sample_id" = seq(1, length(eigengene33))), 
            aes(x = as.factor(Sample_id), y = Loadings)) + 
  geom_line(aes(group = 1), color = "grey40", size = 1) + 
  xlab("Sample ID") + 
  ylim(c(-0.25,0.25)) +
  theme_pubr(base_size = 10) + 
  theme(axis.text.x = element_blank())

ggsave(filename = "/gpfs/hpc/home/liiskolb/transqtl_final/article_figures/eigenegene33.png", e33, width = 3, height = 2)
                                                                                                           "RdYlBu")))(399))
