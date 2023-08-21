install.packages("ggvenn")
library(ggvenn)
setwd("/Users/rusiq/Desktop/IGR/")


del3p = read.csv("DEG_3p_only_shranked.csv")
del3p = del3p[del3p$padj < 0.001 & abs(del3p$log2FoldChange) > 2.5,c("log2FoldChange", "gene_name")]
del9p = read.csv("DEG_9p_only_shranked.csv")
del9p = del9p[del9p$padj < 0.001 & abs(del9p$log2FoldChange) > 2.5,c("log2FoldChange", "gene_name")]
del17p = read.csv("DEG_17p_only_shranked.csv")
del17p = del17p[del17p$padj < 0.001 & abs(del17p$log2FoldChange) > 2.5,c("log2FoldChange", "gene_name")]

"TPX2" %in% del17p$gene_name

gene_lists <- list(del3p = del3p$gene_name, del17p = del17p$gene_name, del9p = del9p$gene_name)

ggvenn(
  gene_lists, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

total_intersection = intersect(del3p$gene_name, intersect(del9p$gene_name, del17p$gene_name))
del9p_del17p = intersect(del9p$gene_name, del17p$gene_name)
del3p_del9p = intersect(del3p$gene_name, del9p$gene_name)
del3p_del17p = intersect(del3p$gene_name, del17p$gene_name)
del3p_only = del3p$gene_name[!(del3p$gene_name %in% del3p_del17p) & !(del3p$gene_name %in% del3p_del9p)]
del9p_only = del9p$gene_name[!(del9p$gene_name %in% del9p_del17p) & !(del9p$gene_name %in% del3p_del9p)]
del17p_only = del17p$gene_name[!(del17p$gene_name %in% del3p_del17p) & !(del17p$gene_name %in% del9p_del17p)]

total_intersection
#ClusterProfiler for unique samples
gene_list = del3p[del3p$gene_name %in% del3p_only,]

dim(deletions_dataframe)
sum(deletions_dataframe$control == 1)

ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl")
entrez_id <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                   filters = "hgnc_symbol",
                   values = total_intersection,
                   mart = ensembl)

gene_list = merge(gene_list, entrez_id, by.x = "gene_name", by.y = "hgnc_symbol")
gene_list <- na.omit(gene_list)

geneList = gene_list$log2FoldChange
names(geneList) = as.character(gene_list$entrezgene_id)
geneList = sort(geneList, decreasing = TRUE)


#ENRICHGO
ego2 <- enrichGO(gene = entrez_id$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               keyType = 'ENTREZID',
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05)

dotplot(ego2, showCategory=15) + ggtitle("Gene enrichment for common genes ")

#GSEGO
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=15) + ggtitle("Gene ontology for del3p samples")

gsex <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
p2 <- heatplot(gsex, foldChange=geneList, showCategory=20) + 
  coord_equal() + 
  ggtitle("Functional classification heatmap for del3p only")
p2



del17p_only

