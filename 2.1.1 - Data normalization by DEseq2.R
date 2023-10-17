### If DEseq2 needs to be insalled ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.16")

## Import libraries ##
library("DESeq2")
library("openxlsx")
library("ggplot2")
library("ggfortify")

## set working directory ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders ##
Folder1 = "Data/Count Tables"
Folder2 = "Data/Normalized"
dir.create(file.path(Folder2),recursive = TRUE, showWarnings = FALSE)

## Files ##

File1 = "Count table - Readcounts.csv"
File2 = "DEseq2 Normalized data.csv"
File3 = "DEseq2 rlog transformed.csv"
File4 = "DEseq2 rlog 2 transformed.csv"

## Dataset ##
data <- read.table(file.path(Folder1,File1), header=T, sep=";", stringsAsFactors=F, row.names=1)

# create condition / treatment #
sample_info <- data.frame(condition = factor(c(c(rep("1_Months",3)),c(rep("3_Months", 3)),c(rep("6_Months", 3)),c(rep("12_Months", 3)),c(rep("18_Months", 3)),c(rep("24_Months", 3)))),row.names = factor(colnames(data)))

## create DE dataset in DESeq2 ##

DESeq.ds <- DESeqDataSetFromMatrix(countData = round(data), colData = sample_info, design = ~ condition)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]

### Data normalization with DESeq2

DESeq.ds <- DESeq(DESeq.ds)
# Get rlogTransformed  #
rlt <- DESeq2::rlogTransformation(DESeq.ds)
rlt2 <- rlog(DESeq.ds, blind=FALSE)
rld <- assay(rlt)
rld2 <- assay(rlt2)

rld <- rlog(DESeq.ds, blind=TRUE)

vsd <- vst(DESeq.ds, blind=FALSE)
DESeq2:::plotPCA.DESeqTransform(vsd)
DESeq2:::plotPCA.DESeqTransform(rlt)

# Estimate size factors
DESeq.ds <- estimateSizeFactors(DESeq.ds)
# Normalize data #
normalized <- counts(DESeq.ds, normalized=TRUE)

## Save normalized data #

write.csv2(normalized,file.path(Folder2, File2))

## Save rlog data ##
write.csv2(rld,file.path(Folder2, File3))
write.csv2(rld2,file.path(Folder2, File4))


# If interested, see PCA of the different normalized data #
## Create PCA plot ##
#p1 <- plotPCA(rlt, intgroup="condition", returnData=TRUE)

## Plot PCA ##
#percentVar <- round(100 * attr(p1, "percentVar"))
#p2 <- ggplot(p1,aes(PC2, PC1))+
#  geom_point(aes(color=condition, fill=condition),size=5,color="black",shape=21,stroke = 1)+
#  geom_text(aes(label = name))+
#  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
#  theme_bw()+
#  theme(panel.grid.minor = element_blank(),legend.title = element_blank())
#p2

#mat <- assay(rld)
#pca<-prcomp(t(mat))

#color_scheme = c("#FF33CC","#00FF00","#009900","#00FFFF","#0066FF","#800080")

#p <-fviz_pca_ind(pca,habillage=DESeq.ds$condition,label="none",geom = "point",col.ind=DESeq.ds$condition,
#                 addEllipses=TRUE, ellipse.level=0.95, palette = color_scheme,repel=TRUE,ellipse.alpha = 0)

#p
