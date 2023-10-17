# -*- coding: utf-8 -*-

### Data normalization ###

## Libraries ##

import os
import pandas as pd
import conorm

## Folders ##

Folder1 = "Data/Count Tables"
Folder2 = "Data/Raw data/RSEM"
Folder3 = "Data/Normalized"

## Files ##

File1 = "Reduced - Count table - Readcounts.csv"
File2 = "RSEM_Sample_M1_R1.txt"

File3 = "Normalized_TMM.csv"
File4 = "Normalized_geTMM.csv"
File5 = "Normalized_MRN.csv"

## Load data ##

# Raw gene count matrix #
df = pd.read_csv(os.path.join(Folder1,File1),sep=";").set_index('Ensembl ID')

# Gene length #

df_length = pd.read_csv(os.path.join(Folder2,File2),sep="\t")

## Modify gene length dataframe ##
df_length = df_length[['gene_id','length']].rename(columns={"gene_id":"Ensembl ID"}).set_index('Ensembl ID')

# add length to df (with length) -> df_wl
df_wl = pd.concat([df,df_length],join='inner',axis=1)

## Normalize to TMM (EdgeR) - Trimmed Means of M-values##
df_tmm = conorm.tmm(df)
nf_tmm  = conorm.tmm_norm_factors(df)

## Normalize to GeTMM - Gene length corrected trimmed mean of M-values ##
df_getmm = conorm.getmm(df_wl, length="length")

## Normalize to mrn (Deseq2) - Median of Ratios Normalization ##

df_mrn = conorm.mrn(df)
nf_mrn  = conorm.mrn_norm_factors(df)

## Save normalized dataframes ##
df_tmm.to_csv(os.path.join(Folder3,File3),sep=";")
df_getmm.to_csv(os.path.join(Folder3,File4),sep=";")
df_mrn.to_csv(os.path.join(Folder3,File5),sep=";")
