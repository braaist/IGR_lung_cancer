library(dplyr)
library(tidyr)
library(ggsci)
library(gplots)
library(ggplot2)
library(tidyverse)

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

deletions_csv = c("SRP313282_deletions.csv",
                  "ERP001058_deletions.csv",
                  "EMTAB6957_deletions.csv",
                  "SRP045225_deletions.csv",
                  "SRP162843_deletions.csv",
                  "SRP090460_deletions.csv")

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
View(tumor_only)

tumor_only = combined_dataframe[combined_dataframe$control == 0,]
only_deleted = tumor_only[(tumor_only$sum_3p <= 0) & (tumor_only$sum_9p <= 0) & (tumor_only$sum_17p <= 0),]
only_deleted

tumor_only$binary_3p = ifelse(tumor_only$sum_3p <= -13, -1, 0)
tumor_only$binary_9p = ifelse(tumor_only$sum_9p <= -11, -1, 0)
tumor_only$binary_17p = ifelse(tumor_only$sum_17p <= -10, -1, 0)


#PLOT DELETION PROPORTIONS IN DIFFERENT GROUPS
tumor_only = deletions_dataframe[deletions_dataframe$control == 0,]

combinations <- tumor_only %>%
  group_by(binary_3p, binary_9p, binary_17p, name) %>%
  summarize(count = n()) %>%
  filter(!(binary_3p == 0 & binary_9p == 0 & binary_17p == 0))

# Create a new column for combining the binary columns
combinations$combination <- paste(combinations$binary_3p,
                                  combinations$binary_9p,
                                  combinations$binary_17p,
                                  sep = ",")
combinations$count_percent = (combinations$count/dim(tumor_only)[1])*100

data.frame(combinations)

ggplot(combinations, aes(x = combination, y = count_percent, fill = name)) +
  geom_bar(stat = "identity", colour="black") +
  labs(x = "Combinations", y = "Count") +
  scale_fill_discrete(name = "Name") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  
  theme_minimal() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tumor_only_with_dels = tumor_only[!(tumor_only$binary_3p == 0 & tumor_only$binary_9p == 0 & tumor_only$binary_17p == 0),]

View(tumor_only_with_dels)

#CREATE OBJECT WITH FINAL SAMPLES 
list_cols = list()
names_file = c()
for (item in names(list_dataframes)){
  row = rownames(tumor_only_with_dels)[startsWith(rownames(tumor_only_with_dels), item)]
  expr_tpm = read.csv(paste0(item, "_quant.tsv"), sep = "\t", row.names = 1)
  current_df = expr_tpm[,colnames(expr_tpm) %in% as.character(map(strsplit(row, split = "\\."), 2))]
  colnames(current_df) <- paste0(item, "_", colnames(current_df))
  print(item)
  list_cols = append(list_cols, list(current_df))  
}
names(list_cols) = names(list_dataframes)

item = "SRP074349"
row = rownames(tumor_only_with_dels)[startsWith(rownames(tumor_only_with_dels), item)]
expr_tpm = read.csv(paste0(item, "_quant.tsv"), sep = "\t", row.names = 1)
current_df = expr_tpm[,colnames(expr_tpm) %in% as.character(map(strsplit(row, split = "\\."), 2))]

View(expr_tpm)

length(list_cols)

#PLOT DELETION PROPORTIONS IN DIFFERENT CANCER TYPES
deletions_dataframe$cancer_type = NA
deletions_dataframe[deletions_dataframe$name == "SRP045225","cancer_type"] <- "SCLC"
deletions_dataframe[deletions_dataframe$name != "SRP045225","cancer_type"] <- "NSCLC"


tumor_only = deletions_dataframe[deletions_dataframe$control == 0,]
combinations <- tumor_only %>%
  group_by(cancer_type, class, binary_3p, binary_9p, binary_17p) %>%
  summarize(count = n()) 

# Create a new column for combining the binary columns
combinations$combination <- paste(combinations$binary_3p,
                                  combinations$binary_9p,
                                  combinations$binary_17p,
                                  sep = ",")
combinations$count_percent = (combinations$count/dim(tumor_only)[1])*100

ggplot(combinations, aes(x = cancer_type, y = count_percent, fill = class)) +
  geom_bar(stat = "identity", position = "fill", colour="black") +
  labs(x = "Combinations", y = "Count") +
  scale_fill_discrete(name = "Name") +
  theme_minimal() +
  scale_fill_aaas() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

