# -*- coding: utf-8 -*-

### Supplementary table 1 ###

## Parameter settings ##

mean_decimals = 0
sd_decimals = 0
sem_decimals = 0

## Functions ##

def Pvalue_annotation(df):
    if df['Pvalue'] <= 0.001:
        return '< 0.001'
    elif df['Pvalue'] <= 0.01:
        return '< 0.01'
    elif df['Pvalue'] <= 0.05:
        return '< 0.05'
    else:
        return 'ns'
    
def Padj_annotation(df):
    if df['Padj'] <= 0.001:
        return '< 0.001'
    elif df['Padj'] <= 0.01:
        return '< 0.01'
    elif df['Padj'] <= 0.05:
        return '< 0.05'
    elif df['Padj'] <= 0.1:
        return '< 0.1'
    else:
        return 'ns'
    
# Libaries #

import os
import pandas as pd

## Folders ##

Folder1 = "Data/Normalized"
Folder2 = "Data/Gene info"
Folder3 = "Results/Differential expression analysis"

Folder_out = "Results/Supplementary tables"
os.makedirs(Folder_out,exist_ok=True)

## Files ##

File1 = "Normalized_TMM.xlsx"
File2 = "Rat_biomart_ensembl_gene.txt"
File3 = "diff all genes.csv"

File_out = "Supplementary table 1.xlsx"

## Load data ##


df_data = pd.read_excel(os.path.join(Folder1,File1))
df_gene_info = pd.read_csv(os.path.join(Folder2,File2))
gene_info_mapping = dict(df_gene_info[['Gene stable ID', 'Gene name']].values)

df_data['Gene name'] = df_data['Ensembl ID'].map(gene_info_mapping)

df_diff = pd.read_csv(os.path.join(Folder3,File3),sep=";",decimal=",").rename(columns={"Unnamed: 0":"Ensembl ID"})
df_diff = df_diff[['Ensembl ID','pvalue','padj']].rename(columns={'pvalue':'Pvalue','padj':'Padj'})

## Calculate mean, sd, and SEM
df_data['1M Mean'] = df_data[['Sample M1 R1','Sample M1 R2','Sample M1 R3']].mean(axis=1)
df_data['1M SD'] = df_data[['Sample M1 R1','Sample M1 R2','Sample M1 R3']].std(axis=1)
df_data['1M SEM'] = df_data[['Sample M1 R1','Sample M1 R2','Sample M1 R3']].sem(axis=1)

df_data['3M Mean'] = df_data[['Sample M3 R1','Sample M3 R2','Sample M3 R3']].mean(axis=1)
df_data['3M SD'] = df_data[['Sample M3 R1','Sample M3 R2','Sample M3 R3']].std(axis=1)
df_data['3M SEM'] = df_data[['Sample M3 R1','Sample M3 R2','Sample M3 R3']].sem(axis=1)

df_data['6M Mean'] = df_data[['Sample M6 R1','Sample M6 R2','Sample M6 R3']].mean(axis=1)
df_data['6M SD'] = df_data[['Sample M6 R1','Sample M6 R2','Sample M6 R3']].std(axis=1)
df_data['6M SEM'] = df_data[['Sample M6 R1','Sample M6 R2','Sample M6 R3']].sem(axis=1)

df_data['12M Mean'] = df_data[['Sample M12 R1','Sample M12 R2','Sample M12 R3']].mean(axis=1)
df_data['12M SD'] = df_data[['Sample M12 R1','Sample M12 R2','Sample M12 R3']].std(axis=1)
df_data['12M SEM'] = df_data[['Sample M12 R1','Sample M12 R2','Sample M12 R3']].sem(axis=1)


df_data['18M Mean'] = df_data[['Sample M18 R1','Sample M18 R2','Sample M18 R3']].mean(axis=1)
df_data['18M SD'] = df_data[['Sample M18 R1','Sample M18 R2','Sample M18 R3']].std(axis=1)
df_data['18M SEM'] = df_data[['Sample M18 R1','Sample M18 R2','Sample M18 R3']].sem(axis=1)


df_data['24M Mean'] = df_data[['Sample M24 R1','Sample M24 R2','Sample M24 R3']].mean(axis=1)
df_data['24M SD'] = df_data[['Sample M24 R1','Sample M24 R2','Sample M24 R3']].std(axis=1)
df_data['24M SEM'] = df_data[['Sample M24 R1','Sample M24 R2','Sample M24 R3']].sem(axis=1)

ordered_list = ['Ensembl ID','Gene name',
                '1M Mean','1M SD','1M SEM',
                '3M Mean','3M SD','3M SEM',
                '6M Mean','6M SD','6M SEM',
                '12M Mean','12M SD','12M SEM',
                '18M Mean','18M SD','18M SEM',
                '24M Mean','24M SD','24M SEM']

ordered_list_Mean = ['1M Mean','3M Mean','6M Mean','12M Mean','18M Mean','24M Mean']
ordered_list_SD = ['1M SD','3M SD','6M SD','12M SD','18M SD','24M SD']
ordered_list_SEM = ['1M SEM','3M SEM','6M SEM','12M SEM','18M SEM','24M SEM']


df_overview = df_data[ordered_list].copy()

df_overview = pd.concat([df_overview.set_index('Ensembl ID'),df_diff.set_index('Ensembl ID')],join='inner',axis=1)

df_overview['Pvalue'] = df_overview.apply(Pvalue_annotation,axis=1)
df_overview['Padj'] = df_overview.apply(Padj_annotation,axis=1)

## Round decimals ##

for key in ordered_list_Mean:
    df_overview[key] = df_overview[key].round(mean_decimals).apply('{0:.0f}'.format)
for key in ordered_list_SD:
    df_overview[key] = df_overview[key].round(sd_decimals).apply('{0:.0f}'.format)
for key in ordered_list_SEM:
    df_overview[key] = df_overview[key].round(sem_decimals).apply('{0:.0f}'.format)

# Save files #
df_overview = df_overview.sort_values(['Padj','Pvalue'], ascending=[True,True])
df_overview.to_excel(os.path.join(Folder_out,File_out),index=True, sheet_name="S1. All genes")