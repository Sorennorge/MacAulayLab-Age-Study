# load Venn diagram package
library("venn")
library("stringr")
library("ggplot2")
library("ggpolypath")

# Set working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders ##
Folder_in = "Data/Venn data"
Folder_out = "Results/Venn diagrams"

dir.create(Folder_out,recursive = TRUE,showWarnings = FALSE)

## Files ##

Set1_file = "Month_1.txt"
Set2_file = "Month_3.txt"
Set3_file = "Month_6.txt"
Set4_file = "Month_12.txt"
Set5_file = "Month_18.txt"
Set6_file = "Month_24.txt"


## Set color list ##

M1 = "#07F2F2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#A422E6"
M24 = "#B3B3B3"

color_list = c(M1,M3,M6,M12,M18,M24)


# Load data #
Set1 <- read.table(file.path(Folder_in,Set1_file), header=T, sep=";")
Set2 <- read.table(file.path(Folder_in,Set2_file), header=T, sep=";")
Set3 <- read.table(file.path(Folder_in,Set3_file), header=T, sep=";")
Set4 <- read.table(file.path(Folder_in,Set4_file), header=T, sep=";")
Set5 <- read.table(file.path(Folder_in,Set5_file), header=T, sep=";")
Set6 <- read.table(file.path(Folder_in,Set6_file), header=T, sep=";")

list1 = c(Set1$Ensembl.ID)
list2 = c(Set2$Ensembl.ID)
list3 = c(Set3$Ensembl.ID)
list4 = c(Set4$Ensembl.ID)
list5 = c(Set5$Ensembl.ID)
list6 = c(Set6$Ensembl.ID)

venn_data <- venn(list(list1,list2,list3,list4,list5,list6),
                  snames=c("1M","3M","6M","12M","18M","24M"),
                  ilabels = FALSE, 
                  ggplot = T,
                  zcolor = color_list,
                  opacity = 0.3,ilcs = 1,sncs = 1, box = FALSE)
venn_data

ggsave(file.path(Folder_out, "VennDiagram_10TMM.png"),plot=venn_data,width=8,height=6,units=c("in"),dpi=800,bg='transparent')
