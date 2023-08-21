if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("enrichplot")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
library("enrichplot")
library("DOSE")
library("clusterProfiler")
library(ggrepel)


#READ FILES AND ADD ANNOTATIONS
results_df = read.csv("/Users/rusiq/Desktop/IGR/DEG_9p_only_shranked.csv")
gene_list = results_df[results_df$padj < 0.001 & abs(results_df$log2FoldChange) > 2.5,c("log2FoldChange", "gene_name")]

ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl")
entrez_id <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                filters = "hgnc_symbol",
                values = gene_list$gene_name,
                mart = ensembl)

gene_list = merge(gene_list, entrez_id, by.x = "gene_name", by.y = "hgnc_symbol")
gene_list <- na.omit(gene_list)

geneList = gene_list$log2FoldChange
names(geneList) = as.character(gene_list$entrezgene_id)
geneList = sort(geneList, decreasing = TRUE)

#ClusterProfiler
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

?gseGO

dotplot(gse, showCategory=15) + ggtitle("Gene ontology for del17p samples")

#FUNCTIONAL CLASSIFICATION HEATMAP
gsex <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
p2 <- heatplot(gsex, foldChange=geneList, showCategory=5) + 
  coord_equal() + coord_flip() +
  ggtitle("Functional classification \n heatmap for del3p")
p2

??heatplot

#NETWORK
gsex <- pairwise_termsim(gsex)
p1 <- emapplot(gsex)
p1



?emapplot

saveRDS(gse, "GSE_object_del9p.rds")

?heatplot



