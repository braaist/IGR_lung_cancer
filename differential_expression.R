library("DESeq2")
library("tidyverse")
library("tximport")
library("ggsci")

setwd("/Users/rusiq/Desktop/IGR/")

#Read counts
tx2gene <- read.csv("/Users/rusiq/Desktop/IGR/files_quant/salmon_tx2gene.tsv", 
                    sep = "\t", header = 0)
head(tx2gene)

#EMTAB6957
PATH = "/Users/rusiq/Desktop/IGR/files_quant/EMTAB6957/"
files = paste0(PATH, list.files(PATH))
names(files) = as.character(map(strsplit(list.files(PATH), "_"),1))

#SRP162843
PATH = "/Users/rusiq/Desktop/IGR/files_quant/SRP162843/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

#ERP001058
PATH = "/Users/rusiq/Desktop/IGR/files_quant/ERP001058/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

#SRP074349
PATH = "/Users/rusiq/Desktop/IGR/files_quant/SRP074349/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

#SRP090460
PATH = "/Users/rusiq/Desktop/IGR/files_quant/SRP090460/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

#SRP313282
PATH = "/Users/rusiq/Desktop/IGR/files_quant/SRP313282/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

#SRP045225
PATH = "/Users/rusiq/Desktop/IGR/files_quant/SRP045225/"
files2 = paste0(PATH, list.files(PATH))
names(files2) = as.character(map(strsplit(list.files(PATH), "_"),1))
files = c(files, files2)

length(files)

#KEEP ONLY DELETIONS AND CONTROLS
deletions_dataframe = read.csv("deletions_combined_dataframe.csv")

files = files[names(files) %in% deletions_dataframe$sample]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)


coldata = deletions_dataframe[,c("X", "control", "name","binary_3p","binary_9p", "binary_17p","class")]
rownames(coldata) = as.character(map(strsplit(coldata$X, "\\."), 2))

#Reorder columns
txi$counts = txi$counts[,order(colnames(txi$counts))]
coldata = coldata[order(rownames(coldata)),]

#ALL SAMPLES VS ALL CONTROLS
coldata$control = factor(coldata$control)
coldata$name = factor(coldata$name)
dds <- DESeqDataSetFromTximport(txi, coldata, ~name + control)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

coldata$name

#PCA
dds_norm <- vst(dds)

plot <- plotPCA(
  dds_norm,
  intgroup = "control"
)
plot + theme_bw() + ggtitle("Included samples PCA plot") + scale_color_npg()

plot <- plotPCA(
  dds_norm,
  intgroup = "class"
)
plot + theme_bw() + ggtitle("Included samples PCA plot") + scale_color_npg()

saveRDS(dds, "all_samples_dds.rds")


#OTHER COMPONENTS

#DEL 3p VS ALL
coldata$binary_3p = factor(coldata$binary_3p)
coldata$name = factor(coldata$name)
dds <- DESeqDataSetFromTximport(txi, coldata, ~name + binary_3p)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_3p_vs_all.csv")

#DEL 9p VS ALL
coldata$binary_9p = factor(coldata$binary_9p)
coldata$name = factor(coldata$name)
dds <- DESeqDataSetFromTximport(txi, coldata, ~name + binary_9p)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_9p_vs_all.csv")

#DEL 17p VS ALL
coldata$binary_17p = factor(coldata$binary_17p)
coldata$name = factor(coldata$name)
dds <- DESeqDataSetFromTximport(txi, coldata, ~name + binary_17p)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_17p_vs_all.csv")

#SUBSET ONLY PAIRED COMPARISONS
dds <- DESeqDataSetFromTximport(txi, coldata, ~name + control)
## subset the coldata
control_samples = colData(dds)[colData(dds)$control == 1,]

coldata_3p_only = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == -1 & colData(dds)$binary_9p == 0 & colData(dds)$binary_17p == 0),]
coldata_3p_only = rbind(coldata_3p_only, control_samples)
coldata_9p_only = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == 0 & colData(dds)$binary_9p == -1 & colData(dds)$binary_17p == 0),]
coldata_9p_only = rbind(coldata_9p_only, control_samples)
coldata_17p_only = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == 0 & colData(dds)$binary_9p == 0 & colData(dds)$binary_17p == -1),]
coldata_17p_only = rbind(coldata_17p_only, control_samples)
coldata_3p_17p = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == -1 & colData(dds)$binary_9p == 0 & colData(dds)$binary_17p == -1),]
coldata_3p_17p = rbind(coldata_3p_17p, control_samples)
coldata_3p_9p = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == -1 & colData(dds)$binary_9p == -1 & colData(dds)$binary_17p == 0),]
coldata_3p_9p = rbind(coldata_3p_9p, control_samples)
coldata_17p_9p = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == 0 & colData(dds)$binary_9p == -1 & colData(dds)$binary_17p == -1),]
coldata_17p_9p = rbind(coldata_17p_9p, control_samples)
coldata_all = colData(dds)[(colData(dds)$control == 0 & colData(dds)$binary_3p == -1 & colData(dds)$binary_9p == -1 & colData(dds)$binary_17p == -1),]
coldata_all = rbind(coldata_all, control_samples)

#3p

colData(dds)
dds_use <- dds[,coldata_3p_only$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_3p_only.csv")

#9p
dds_use <- dds[,coldata_9p_only$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_9p_only.csv")

#17p
dds_use <- dds[,coldata_17p_only$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_17p_only.csv")

#3p and 9p
dds_use <- dds[,coldata_3p_9p$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_3p_9p_only.csv")

#3p and 17p
dds_use <- dds[,coldata_3p_17p$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_3p_17p_only.csv")

#9p and 17p
dds_use <- dds[,coldata_17p_9p$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_9p_17p_only.csv")

#all
dds_use <- dds[,coldata_all$samples]
dds_use <- DESeq(dds_use)
res <- results(dds_use)
res <- na.omit(res)
res = merge(data.frame(res), gene_name, by = "row.names", all.x = TRUE)
write.csv(res, "/Users/rusiq/Desktop/IGR/DEG_all.csv")


gene_name = read.csv("/Users/rusiq/Desktop/IGR/EMTAB6957_raw_counts.tsv", sep = "\t")
gene_name = gene_name[,c("gene_id", "gene_name")]
rownames(gene_name) = gene_name$gene_id





