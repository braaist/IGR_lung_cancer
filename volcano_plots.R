library("apeglm")
library("De")
library("EnhancedVolcano")

dds_use <- dds[,coldata_17p_only$samples]
dds_use <- DESeq(dds_use)
dds_result = results(dds_use)

resLFC <- lfcShrink(dds = dds_use, 
                    res = dds_result,
                    type = "normal",
                    coef = "control_1_vs_0") 

??lfcShrink

resLFC = na.omit(resLFC)
resLFC = merge(data.frame(resLFC), gene_name, by = "row.names", all.x = TRUE)
resLFC$log2FoldChange = -resLFC$log2FoldChange

write.csv(resLFC, "/Users/rusiq/Desktop/IGR/DEG_17p_only_shranked.csv")


EnhancedVolcano(toptable = resLFC,              
                x = "log2FoldChange",
                y = "padj",
                lab = resLFC$gene_name
)

results_9p = read.csv("/Users/rusiq/Desktop/IGR/DEG_9p_only_shranked.csv")

EnhancedVolcano(toptable = results_9p,
                x = "log2FoldChange",
                y = "padj",
                lab = results_9p$gene_name,
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 0.001,
                FCcutoff = 2.5, 
                title = "del_3p vs. control samples"
)

results_9p[results_9p$padj < 0.001 & abs(results_9p$log2FoldChange) > 2.5,]



