# MacAulayLab: Age study of the choroid plexus
The work and scripts are done by the MacAulay Lab.\
All programs used are free and open-source.
In the interest of open science and reproducibility, all data and source code used in our research is provided here.\
Feel free to copy and use code, but please cite:\
(coming soon) \
*Remember* rewrite file_names and folder_names suitable for your pipeline.

## The RNAseq and Analysis follows these steps:
## Raw data analysis - Library Build, Mapping and Quantification ##
The analysis uses RNA STAR for mapping and RSEM for TPM quantification.
### RNA-STAR and RSEM Library build and indexing ###

0.1.1 - RNA_STAR_Indexing.sh \
0.2.1 - RSEM_Indexing.sh

### RNA-STAR Mapping and RSEM quantification ###

0.1.2 -RNA_STAR_Analysis.sh \
0.2.2 - RSEM_Analysis.sh

## Gene count tables ##
1.1.1 - Raw data - GeneCounts.py \
1.2.1 - Raw data - RSEM.py
## Data normalization ##
2.1.1 - Data normalization by DEseq2.R \
2.1.2 - Data normalization - Verify normalization.py

## Data visulization ##
### Density plot ###
2.2.1 - Density plots for normalized data.py
### PCA ###
3.1.2 - PCA plots for raw data.py \
3.1.3 - PCA plots for raw data (sorted data).py

### VennDiagram ###
4.1.1 - Create gene lists for venn.py \
4.2.1 - Venndiagram of overlapping genes.R

## Differential expression analysis - LRT method ##
5.1.1 - Differential expression analysis.R

## Panther - GO (Gene Ontology) Analysis ##
### Data cleaning ###
There might be some manual data verification and sorting depending on your own data. \
6.0.1 - MF data cleaning.py \
6.0.2 - PC data cleaning.py \
6.0.3 - BF data cleaning.py
### GO Enrichment analysis ###
6.1.1 - Enrichment analysis - PC.py \
6.1.2 - Enrichment analysis - MF.py \
6.1.3 - Enrichment analysis - BF.py

### Pitchart plots ###
6.2.1 - Enrichment Plot - PC.py \
6.2.2 - Enrichment plot - MF.py \
6.2.3 - Enrichment Plot - BF.py

## Time series plot of MF ##
7.1.1 - MF - Time series - Generate data.py \
7.1.2 - MF - time series - Generate plots.py

## Transport analyses ##
8.1.1 - Transporter - Dataframe - Heatmap.py \
8.1.2 - Transporter DE - Heatmap plot.py \
8.2.1 - Transporter DE - Usual suspects - Time series plot.py \
8.2.2 - Transporter DE - Usual suspects -  Expression barplot.py

## Metabolic processes ##
9.1.1 - Metabolic processes.py \
9.2.1 - Metabolic gene analysis.py

## Toxin removal analysis ##
10.1.1 - Toxin removal analysis.py

## Suppelemtary tables and statistical significance ##
11.1.1 - Supplementary table 1.py \
12.1.1 - Statistical significance - Usual suspects - Time series.py \
12.1.2 - Statistical significance - MF.py \
12.1.3 - Statistical significance - Toxin removal.py \
12.1.4 - Statistical significance - Metabolic processes.py \
13.1.1 - Supplementary table 1 - Sheet 2 and 3.py
