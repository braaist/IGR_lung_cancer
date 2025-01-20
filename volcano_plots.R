library(apeglm)
library(EnhancedVolcano)
library(DESeq2)
library(tidyverse)
library(tximport)
library(ggsci)
library(ggvenn)

setwd("/Users/rusiq/Desktop/IGR/")
source(paste0(getwd(),"/CustomFunctions.R"))

txi_3p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_3p_res_output_initial.csv")
txi_9p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_9p_res_output_initial.csv")
txi_17p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_17p_res_output_initial.csv")

EnhancedVolcano(toptable = txi_3p_res,
                x = "log2FoldChange",
                y = "padj",
                lab = txi_3p_res$gene_name,
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 0.05,
                FCcutoff = 2, 
                title = "del_3p vs. control samples"
)

EnhancedVolcano(toptable = txi_9p_res,
                x = "log2FoldChange",
                y = "padj",
                lab = txi_9p_res$gene_name,
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 0.05,
                FCcutoff = 2, 
                title = "del_9p vs. control samples"
)

EnhancedVolcano(toptable = txi_17p_res,
                x = "log2FoldChange",
                y = "padj",
                lab = txi_17p_res$gene_name,
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 0.05,
                FCcutoff = 2, 
                title = "del_17p vs. control samples"
)






