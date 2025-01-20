library(apeglm)
library(EnhancedVolcano)
library(DESeq2)
library(tidyverse)
library(tximport)
library(ggsci)
library(ggvenn)

setwd("/Users/rusiq/Desktop/IGR/")
source(paste0(getwd(),"/CustomFunctions.R"))

#PERFORM DIFFERENTIAL EXPRESSION 
setwd("/Users/rusiq/Desktop/IGR/")
#load dataframe with deletion annotations
deletions_dataframe = read.csv("deletions_combined_dataframe.csv")
deletions_dataframe$control <- factor(deletions_dataframe$control, levels = c(1, 0), labels = c("Control", "Tumor")) 

#get file names
files = c()
for (file in unique(deletions_dataframe$name)){
  file_names = GetFileNames("/Users/rusiq/Desktop/IGR/files_quant/", file)
  files = c(files, file_names)
}

#read trancript_to_gene annotation
tx2gene <- read.csv("/Users/rusiq/Desktop/IGR/files_quant/salmon_tx2gene.tsv", 
                    sep = "\t", header = 0)


#Get only unique deletions
deletions_17p_only = deletions_dataframe[deletions_dataframe["binary_3p"] == 0 & deletions_dataframe["binary_9p"] == 0 & deletions_dataframe["binary_17p"] == -1,]
deletions_3p_only = deletions_dataframe[deletions_dataframe["binary_3p"] == -1 & deletions_dataframe["binary_9p"] == 0 & deletions_dataframe["binary_17p"] == 0,]
deletions_9p_only = deletions_dataframe[deletions_dataframe["binary_3p"] == 0 & deletions_dataframe["binary_9p"] == -1 & deletions_dataframe["binary_17p"] == 0,]

#Get only controls 
controls_only = deletions_dataframe[deletions_dataframe["control"] == "Control" & deletions_dataframe["binary_3p"] == 0 & deletions_dataframe["binary_9p"] == 0 & deletions_dataframe["binary_17p"] == 0,]

#Write deletion files
write.csv(deletions_3p_only, "deletions_3p_only_samples.csv")
write.csv(deletions_9p_only, "deletions_9p_only_samples.csv")
write.csv(deletions_17p_only, "deletions_17p_only_samples.csv")
write.csv(controls_only, "controls_only_samples.csv")

#Join deletions with controls for DESeq
deletions_3p_only = rbind(deletions_3p_only, controls_only)
deletions_9p_only = rbind(deletions_9p_only, controls_only)
deletions_17p_only = rbind(deletions_17p_only, controls_only)


#get file paths for tximport
files_17p = files[names(files) %in% deletions_17p_only$sample]
files_3p = files[names(files) %in% deletions_3p_only$sample]
files_9p = files[names(files) %in% deletions_9p_only$sample]

#Import 3p, 9p and 17p files
txi_17p = TximportWrapper(files = files_17p, tx2gene = tx2gene)
txi_9p = TximportWrapper(files = files_9p, tx2gene = tx2gene)
txi_3p = TximportWrapper(files = files_3p, tx2gene = tx2gene)

#Check whether order is normal
deletions_17p_only = deletions_17p_only[order(deletions_17p_only$sample),]
colnames(txi_17p$counts) == deletions_17p_only$sample
deletions_9p_only = deletions_9p_only[order(deletions_9p_only$sample),]
colnames(txi_9p$counts) == deletions_9p_only$sample
deletions_3p_only = deletions_3p_only[order(deletions_3p_only$sample),]
colnames(txi_3p$counts) == deletions_3p_only$sample

#Get DESeq results
txi_17p_dds <- DESeqWrapper(txi_17p, deletions_17p_only, counts_QC_treshold = 10, rowsums_QC_treshold = 50)
txi_9p_dds <- DESeqWrapper(txi_9p, deletions_9p_only, counts_QC_treshold = 10, rowsums_QC_treshold = 50)
txi_3p_dds <- DESeqWrapper(txi_3p, deletions_3p_only, counts_QC_treshold = 10, rowsums_QC_treshold = 50)

#Clean garbage
rm(txi_3p, txi_9p, txi_17p)
gc()

#Save initial objects
saveRDS(txi_17p_dds, "txi_17p_dds.Rds")
saveRDS(txi_9p_dds, "txi_9p_dds.Rds")
saveRDS(txi_3p_dds, "txi_3p_dds.Rds")

#txi_17p_dds = readRDS("txi_17p_dds.Rds")
#txi_9p_dds = readRDS("txi_9p_dds.Rds")
#txi_3p_dds = readRDS("txi_3p_dds.Rds")

#Get the dataframe for gene names conversion
tx2gene <- read.csv("/Users/rusiq/Desktop/IGR/files_quant/salmon_tx2gene.tsv", 
                    sep = "\t", header = 0)
colnames(tx2gene) <- c("transcript", "gene_id", "gene_name")

tx2gene <- tx2gene[!duplicated(tx2gene$gene_id), ]
rownames(tx2gene) <- tx2gene$gene_id

#Not filter DESeq results
padj_threshold = 1
LFC_threshold = 0
txi_17p_initial <- DESeqResultProcessor(txi_17p_dds, padj_threshold, LFC_threshold, tx2gene)
txi_9p_initial <- DESeqResultProcessor(txi_9p_dds, padj_threshold, LFC_threshold, tx2gene)
txi_3p_initial <- DESeqResultProcessor(txi_3p_dds, padj_threshold, LFC_threshold, tx2gene)

#Filter deseq results
padj_threshold = 0.05
LFC_threshold = 2
txi_17p_filtered <- DESeqResultProcessor(txi_17p_dds, padj_threshold, LFC_threshold, tx2gene)
txi_9p_filtered <- DESeqResultProcessor(txi_9p_dds, padj_threshold, LFC_threshold, tx2gene)
txi_3p_filtered <- DESeqResultProcessor(txi_3p_dds, padj_threshold, LFC_threshold, tx2gene)

rm(txi_3p_dds, txi_9p_dds, txi_17p_dds)
gc()

#write initial datasets
write.csv(txi_17p_initial, "/Users/rusiq/Desktop/IGR/txi_17p_res_output_initial.csv")
write.csv(txi_9p_initial, "/Users/rusiq/Desktop/IGR/txi_9p_res_output_initial.csv")
write.csv(txi_3p_initial, "/Users/rusiq/Desktop/IGR/txi_3p_res_output_initial.csv")

write.csv(txi_17p_filtered, "/Users/rusiq/Desktop/IGR/txi_17p_res_output_filtered.csv")
write.csv(txi_9p_filtered, "/Users/rusiq/Desktop/IGR/txi_9p_res_output_filtered.csv")
write.csv(txi_3p_filtered, "/Users/rusiq/Desktop/IGR/txi_3p_res_output_filtered.csv")

#Venn diagram
txi_3p_res <- read.csv("/Users/rusiq/Desktop/IGR/txi_3p_res_output_filtered.csv")
txi_9p_res <- read.csv("/Users/rusiq/Desktop/IGR/txi_9p_res_output_filtered.csv")
txi_17p_res <- read.csv("/Users/rusiq/Desktop/IGR/txi_17p_res_output_filtered.csv")

gene_lists <- list(del3p = txi_3p_res$gene_name, del17p = txi_17p_res$gene_name, del9p = txi_9p_res$gene_name)

ggvenn(
  gene_lists, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)





