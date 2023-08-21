#GENE EXPRESSION BOXPLOTS
install.packages("ggpubr")
library(ggpubr)
library(ggsci)
library(dplyr)
library(DESeq2)
setwd("/Users/rusiq/Desktop/IGR/")

dds = readRDS("all_samples_dds.rds")

dds_norm <- vst(dds)

ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl")
entrez_id <- getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol"),
                   filters = "ensembl_gene_id_version",
                   values = rownames(rowData(dds)),
                   mart = ensembl)

#LOAD DELETIONS ANNOTATION
deletions_dataframe = read.csv("deletions_combined_dataframe.csv")
deletions_dataframe

#ORDER GENE COUNTS
gene_counts = gene_counts[,order(colnames(gene_counts))]
deletions_dataframe = deletions_dataframe[order(deletions_dataframe$sample),]
names(gene_counts) == deletions_dataframe$sample

dds_norm
#calculate wilcoxon significance
compare_means(counts ~ group, data = m, method = "wilcox.test",  p.adjust.method = "BH")
my_comparisons <- list( c("control", "del17p"), c("control", "del3p"), c("control", "del3p_del17p"),
                        c("control", "del3p_del9p"), c("control", "del3p_del9p_del17p"), c("control", "del9p"),
                        c("control", "del9p_del17p"))

gene_counts <- as.matrix(assay(dds_norm[entrez_id[entrez_id$hgnc_symbol == "PLK1","ensembl_gene_id_version"],]))
m <- list(counts = as.numeric(gene_counts), group = as.factor(deletions_dataframe$class))
m <- as_tibble(m)
q <- ggplot(m, aes(group, counts, fill = group)) + geom_boxplot() + geom_jitter(width = 0.1) +
     theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
     scale_fill_npg() +
    labs(x = "Experimental Group", y = "VST counts ", title = "Expression of TPX2") +
     stat_compare_means(comparisons = my_comparisons, label = "p.signif") 
q



del9p_only = del9p[!(del9p$gene_name %in% del9p_del17p) & !(del9p$gene_name %in% del3p_del9p),]
del9p_only[del9p_only$log2FoldChange < 0,"gene_name"]

del9p_only

