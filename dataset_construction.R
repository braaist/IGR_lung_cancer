library(dplyr)
library(tidyr)
library(ggsci)
library(gplots)
library(ggplot2)
library(tidyverse)

#LOAD MARKERS FOR HEATMAP
markers_3p = c("CISH", "HEMK1", "CACNA2D2", "TMEM115", 
               "CYB561D2", "NPRL2", "BLU", "RASSF1", 
               "FUS1", "HYAL2", "HYAL1", "HYAL3", 
               "NAA80", "IFRD2", "SEMA3B", 
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

#COMBINE DELETIONS DATAFRAMES FROM CASPER
deletions_csv = c("SRP313282_deletions.csv",
                  "ERP001058_deletions.csv",
                  "SRP045225_deletions.csv",
                  "SRP162843_deletions.csv",
                  "SRP090460_deletions.csv",
                  "EMTAB6957_deletions.csv",
                  "SRP074349_deletions.csv")

list_dataframes = list()
names = c()
for (name in deletions_csv){
  file = read.csv(name, row.names = 1) 
  file$name = strsplit(name, "_")[[1]][1]
  file$sum_3p = rowSums(file[,colnames(file) %in% markers_3p])
  file$sum_9p = rowSums(file[,colnames(file) %in% markers_9p])
  file$sum_17p = rowSums(file[,colnames(file) %in% markers_17p])
  
  list_dataframes = append(list_dataframes, list(file))
  names = c(names, strsplit(name, "_")[[1]][1])
}

names(list_dataframes) = names
combined_dataframe <- do.call(rbind, list_dataframes)
combined_dataframe$GNAI1 = NULL
View(combined_dataframe)

#ONLY TUMOR DELETIONS

combined_dataframe$binary_3p = ifelse(combined_dataframe$sum_3p <= -13, -1, 0)
combined_dataframe$binary_9p = ifelse(combined_dataframe$sum_9p <= -11, -1, 0)
combined_dataframe$binary_17p = ifelse(combined_dataframe$sum_17p <= -10, -1, 0) 

tumor_only = combined_dataframe[combined_dataframe$control == 0,]
control_only = combined_dataframe[combined_dataframe$control == 1,]
tumor_only = tumor_only[(tumor_only$sum_3p <= 0) & (tumor_only$sum_9p <= 0) & (tumor_only$sum_17p <= 0),]

#filter 0,0,0 deletions
tumor_only = tumor_only[!(tumor_only$binary_3p == 0 & tumor_only$binary_9p == 0 & tumor_only$binary_17p == 0),]
View(tumor_only)
View(control_only)

#FINAL DATAFRAME WITH CONTROLS AND REQUIRED DELETIONS
final_deletions = rbind(tumor_only, control_only)

#write.csv(final_deletions, "deletions_combined_dataframe.csv")

#LOAD EXPRESSION DATA
list_cols = list()
names_file = c()

for (item in names(list_dataframes)){
  row = rownames(final_deletions)[startsWith(rownames(final_deletions), item)]
  expr_raw_salmon = read.csv(paste0(item, "_raw_counts.tsv"), sep = "\t", row.names = 1)
  current_df = expr_raw_salmon[,colnames(expr_raw_salmon) %in% as.character(sapply(strsplit(row,"\\."), `[`, 2))]
  colnames(current_df) <- paste0(item, "_", colnames(current_df))
  print(item)
  list_cols = append(list_cols, list(current_df))  
}
names(list_cols) = names(list_dataframes)

initial = list_cols$SRP313282
for (item in list_cols[-1]){
#  print(sum(rownames(initial) %in% rownames(item)))
  initial = merge(initial, item, by = "row.names", all = TRUE)
  rownames(initial) = initial$Row.names
  initial$Row.names = NULL
#  print(names(item))
}

initial = merge(initial, expr_raw_salmon[,"gene_name", drop = FALSE], by = "row.names", all = TRUE)
write.csv(initial, "/Users/rusiq/Desktop/IGR/selected_deletions_controls_counts.csv")

