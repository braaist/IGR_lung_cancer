library(ggplot2)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install(c('biomaRt', 'limma', 'GO.db', 'org.Hs.eg.db', 'GOstats',  'GenomicRanges'))
install.packages("devtools")

require(devtools)
install_github("akdess/CaSpER")

library(CaSpER)
BiocManager::install("GEOquery")
library(GEOquery)

data("yale_meningioma")
yale_meningioma$genoMat

#SRP090460
setwd("/Users/rusiq/Desktop/IGR/SRP090460_quant/")
file_data <- list()
tsv_files <- list.files(pattern = "\\.sf$")
for (file in tsv_files) {
  df <- read.delim(file, sep = "\t", header = TRUE)
  df_subset <- df[c("Name", "TPM")]
  colnames(df_subset) <- c("Name", strsplit(file, "\\.")[[1]][1])
  file_data[[file]] <- df_subset
}

merged_df <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), file_data)
rownames(merged_df) = merged_df$Name
merged_df$Name <- NULL

merged_df = expr_df

View(merged_df)
data(hg38_cytoband)
View(cytoband_hg38)

centromere = read.csv("/Users/rusiq/Desktop/IGR/centromere.tsv", sep = "\t")
centromere

annotation <- generateAnnotation(id_type="ensembl_gene_id_version", 
                                  genes=rownames(merged_df), 
                                  ishg19=F, centromere)

merged_df = merged_df[annotation$Gene, ]

loh <- readBAFExtractOutput(path="/Users/rusiq/Desktop/IGR/SRP090460_snp/",
                            sequencing.type="bulk", suffix = "snp")

names(loh) <- gsub(".snp", "", names(loh))

#Get metadata for dataset
geo_id <- "GSE87340"
gse <- getGEO(geo_id,GSEMatrix=TRUE)
metadata = pData(phenoData(gse[[1]]))

SRX = c()
for (item in metadata$relation.1) {
  split_result <- strsplit(item, "=")
  second_part <- split_result[[1]][2]
  SRX = c(SRX, second_part)
}
rownames(metadata) = SRX
sra_metadata = read.csv("/Users/rusiq/Downloads/SRP090460_run.csv")
rownames(sra_metadata) = sra_metadata$Experiment
merged_data <- merge(metadata, sra_metadata, by = "row.names", all = TRUE)
View(sra_metadata)
View(metadata)
View(merged_data)

control_vec <- merged_data[merged_data$characteristics_ch1.5 == "cell type: NL", "Run"]

loh.name.mapping <- data.frame(
  loh.name = colnames(merged_df),
  sample.name = colnames(merged_df)
)

dim(merged_df)

object <- CreateCasperObject(raw.data=merged_df,
                             loh.name.mapping=loh.name.mapping, 
                             sequencing.type="bulk", 
                             cnv.scale=3, loh.scale=3, 
                             matrix.type="normalized", expr.cutoff=4.5,
                             annotation=annotation, method="iterative", 
                             loh=loh, filter="median",  
                             control.sample.ids=control_vec, 
                             cytoband=cytoband,
                             genomeVersion = "hg38")

final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")
finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75) 

gamma <- 6
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]


obj <- final.objects[[9]]
plotHeatmap(object=obj, fileName="heatmap_soi.png",
            cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

plotLargeScaleEvent(object=obj, fileName="large.scale.events.png") 

plotBAFAllSamples(loh = obj@loh.median.filtered.data,  fileName="LOHAllSamples.png") 

plotBAFOneSample (object, fileName="LOHPlotsAllScales.pdf") 

plotGEAndBAFOneSample(object=obj, cnv.scale=3, loh.scale=3, sample= "SRR4296071")



all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
                IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

test_df = data.frame(t(data.frame(rna.matrix[rownames(rna.matrix) %in% c("NPRL2", "CACNA2D2", "RASSF1", "FUS1", "CDKN2A", "DMRT1", "MTAP", "TP53", "KCTD11", "PTPRH"),])))
test_df$control = (as.numeric(rownames(test_df) %in% control_vec))

deletions_dataframe = read.csv("/Users/rusiq/Desktop/IGR/deletions_combined_dataframe.csv")
test_df = deletions_dataframe
rownames(test_df) = test_df$sample
test_df = test_df[,!colnames(test_df) %in% c("X.2", "X.1", "X", "name", "sum_3p", "sum_9p", "sum_17p", 
          "binary_3p", "binary_9p", "binary_17p", "sample")]
test_df <- test_df[order(test_df$control),]
test_df$class = NULL

test_melt <- melt(as.matrix(test_df))

ggplot(test_melt, aes(x = X1, y = X2, fill = value)) +
  geom_tile() + scale_fill_gradientn(colours = c(low = "blue", mid = "white", high = "red")) +
  coord_fixed() + 
  coord_flip() +
  labs(title = "Binary gene heatmap",
       x = "Samples",
       y = "Genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

View(deletions_dataframe[,c("class", "control")])



?scale_fill_gradient

rna.matrix[rownames(rna.matrix) == "CYB561D2",]

rna.matrix[rownames(rna.matrix) == "TP53",]

sum(data.frame(rna.matrix)$SRR4296064 == "-1")
control_vec

loss.final[loss.final["seqnames"] == "17p",]

?plotHeatmap
