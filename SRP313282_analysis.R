library(ggplot2)
library(CaSpER)
library(GEOquery)
library(GenomicRanges)

setwd("/Users/rusiq/Desktop/IGR/")

#ANALYZE METADATA
geo_id <- "GSE171415"
gse <- getGEO(geo_id,GSEMatrix=TRUE)
metadata = pData(phenoData(gse[[1]]))
View(metadata)

#GET CONTROL SAMPLE IDS
sra_metadata = read.csv("/Users/rusiq/Downloads/SRP313282_run.csv")
sra_metadata = sra_metadata[,c("Run", "SampleName")]
merged_data <- merge(metadata, sra_metadata, 
                     by.x = "row.names", by.y = "SampleName",
                     all = TRUE)

merged_data <- na.omit(merged_data)
control_vec <- merged_data[merged_data$source_name_ch1 == "Normal lung tissue", "Run"]

#READ EXPRESSION
expr_df = read.csv("/Users/rusiq/Desktop/IGR/SRP313282_quant.tsv", sep = "\t")
centromere = read.csv("/Users/rusiq/Desktop/IGR/centromere.tsv", sep = "\t")
annotation <- generateAnnotation(id_type="ensembl_gene_id_version", 
                                 genes=expr_df$gene_id, 
                                 ishg19=F, centromere)
annotation <- annotation[annotation$GeneSymbol != "",]

#Keep genes mapped to chr coords
rownames(expr_df) = expr_df$gene_id
expr_df$gene_id = NULL
expr_df$gene_name = NULL
expr_df <- expr_df[annotation$Gene, ]

#LOAD LOH
loh <- readBAFExtractOutput(path="/Users/rusiq/Desktop/IGR/SRP313282_snp/",
                            sequencing.type="bulk", suffix = "snp")
names(loh) <- gsub(".snp", "", names(loh))

loh.name.mapping <- data.frame(
  loh.name = colnames(expr_df),
  sample.name = colnames(expr_df)
)

dim(loh.name.mapping)
dim(expr_df)
dim(annotation)

names(loh)[!names(loh) %in% colnames(expr_df)]

#CREATE CASPER
object <- CreateCasperObject(raw.data=expr_df,
                             loh.name.mapping=loh.name.mapping, 
                             sequencing.type="bulk", 
                             cnv.scale=3, loh.scale=3, 
                             matrix.type="normalized", expr.cutoff=3,
                             annotation=annotation, method="iterative", 
                             loh=loh, filter="median",  
                             control.sample.ids=control_vec, 
                             cytoband=cytoband,
                             genomeVersion = "hg38")
final.objects <- runCaSpER(object, removeCentromere=T, 
                           cytoband=cytoband, method="iterative")

saveRDS(final.objects, "SRP313282_finals.rds")


obj <- final.objects[[1]]
plotLargeScaleEvent(object=obj, fileName="large.scale.events.png") 

finalChrMat <- extractLargeScaleEvents(final.objects, thr=1) 

gamma <- 8
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]

all.summary<- rbind(loss.final, gain.final) 
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
                IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation, 
                                   keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation[,2])

all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, 
                          all.samples=all.samples, genes.ann=genes.ann)

markers_3p = c("CISH", "HEMK1", "CACNA2D2", "TMEM115", 
               "CYB561D2", "NPRL2", "BLU", "RASSF1", 
               "FUS1", "HYAL2", "HYAL1", "HYAL3", 
               "NAA80", "IFRD2", "SEMA3B", "GNAI1", 
               "GNAT1", "SEMA3F", "RBM6")

markers_9p = c("IFNB1", "IFNW1", "IFNA21", "IFNA4",
               "IFNA7", "IFNA10", "IFNA16", "IFNA17",
               "IFNA14", "IFNA5", "IFNA6", "MTAP",
               "IFNA13", "IFNA2", "IFNA8", "IFNA11",
               "IFNE", "CDKN2A")

markers_17p = c("AIPL1", "SLC13A5", "XAF1", "FBX039",
                "ALOX12P2", "DLG4", "EIF5A", "FGF11",
                "POLR2A", "TNFSF13", "SENP3", "SOX15",
                "SHBG", "ATP1B2", "TP53", "WRAP53",
                "ALOXE3")

markers = c(markers_3p, markers_9p, markers_17p)

test_df = data.frame(t(data.frame(rna.matrix[rownames(rna.matrix) %in% markers,])))
test_df$control = (as.numeric(rownames(test_df) %in% control_vec))
test_df <- test_df[order(test_df$control),]

test_melt <- melt(as.matrix(test_df))

ggplot(test_melt, aes(x = X2, y = X1, fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c(low = "blue", 
                                   mid = "white", 
                                   high = "red")) +
  labs(title = "Binary gene heatmap",
       x = "Variable 1",
       y = "Variable 2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write.csv(test_df, "SRP313282_deletions.csv")

obj <- final.objects[[1]]
plotLargeScaleEvent(object=obj, fileName="large.scale.events.png") 

View(test_df)










