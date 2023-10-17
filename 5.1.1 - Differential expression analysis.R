### If DEseq2 needs to be insalled ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.16")

## Import libraries ##
library("DESeq2")
library("openxlsx")
library("ggplot2")

## set working directory ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders ##
Folder1 = "Data/Count Tables"
Folder2 = "Results/Differential expression analysis"
dir.create(file.path(Folder2),recursive = TRUE, showWarnings = FALSE)

## Files ##

File1 = "Count table - Readcounts.csv"
File_diff_1_24 = "LRT_diff_input_list.csv"

## Dataset ##
data <- read.table(file.path(Folder1,File1), header=T, sep=";", stringsAsFactors=F, row.names=1)

# create condition / treatment #
sample_info <- data.frame(Age = factor(c(c(rep("1_Months",3)),c(rep("3_Months", 3)),c(rep("6_Months", 3)),c(rep("12_Months", 3)),c(rep("18_Months", 3)),c(rep("24_Months", 3)))),row.names = factor(colnames(data)))

## create DE dataset in DESeq2 ##

DESeq.ds <- DESeqDataSetFromMatrix(countData = round(data), colData = sample_info, design = ~ Age)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]
colData(DESeq.ds)$Age <- relevel(colData(DESeq.ds)$Age, "1_Months")

### Differential expression analysis with DESeq2

#colData(DESeq.ds)$treatment <- relevel(colData(DESeq.ds)$condition, "Zucker_Obese")

DESeq_rat_age <- DESeq(DESeq.ds,test="LRT", reduced=~1)
#design(DESeq_rat_age)
resultsNames(DESeq_rat_age)

res_3_1 <- results(DESeq_rat_age,name="Age_3_Months_vs_1_Months",independentFiltering = TRUE, alpha = 0.05)
res_6_1 <- results(DESeq_rat_age,name="Age_6_Months_vs_1_Months",independentFiltering = TRUE, alpha = 0.05)
res_12_1 <- results(DESeq_rat_age,name="Age_12_Months_vs_1_Months",independentFiltering = TRUE, alpha = 0.05)
res_18_1 <- results(DESeq_rat_age,name="Age_18_Months_vs_1_Months",independentFiltering = TRUE, alpha = 0.05)
res_24_1 <- results(DESeq_rat_age,name="Age_24_Months_vs_1_Months",independentFiltering = TRUE, alpha = 0.05)
summary(res_3_1)
summary(res_6_1)
summary(res_12_1)
summary(res_18_1)
summary(res_24_1)

### store the genes of interest ###
res_3_1.sorted <- res_3_1[order(res_3_1$padj), ]
res_6_1.sorted <- res_6_1[order(res_6_1$padj), ]
res_12_1.sorted <- res_12_1[order(res_12_1$padj), ]
res_18_1.sorted <- res_18_1[order(res_18_1$padj), ]
res_24_1.sorted <- res_24_1[order(res_24_1$padj), ]
# 
DGEgenes_3_1 <- rownames(subset(res_3_1.sorted, padj < 0.05))
DGEgenes_6_1 <- rownames(subset(res_6_1.sorted, padj < 0.05))
DGEgenes_12_1 <- rownames(subset(res_12_1.sorted, padj < 0.05))
DGEgenes_18_1 <- rownames(subset(res_18_1.sorted, padj < 0.05))
DGEgenes_24_1 <- rownames(subset(res_24_1.sorted, padj < 0.05))

# # Write all expressed genes to matrix -> same for 3vs1,6vs1 etc.
write.csv2(res_24_1,file.path(Folder2,"all_genes_24.csv"))

# Write differentially expressed genes to file #
write.csv2(DGEgenes_24_1,file.path(Folder2,File_diff_1_24))
