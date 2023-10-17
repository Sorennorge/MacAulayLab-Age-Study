# -*- coding: utf-8 -*-

### Convert all RSEM output files to csv, extracting TPM column ###

## Import libraries ##

import os
import pandas as pd

## Folders #
Folder1 = "Data/Raw data/RSEM"
Folder2 = "Data/Raw data/TPM"
os.makedirs(Folder2,exist_ok=True)
Folder3 = "Data/Count Tables"
os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "Count table - TPM.csv"
File2 = "Reduced - Count table - TPM.csv"

## for each sample file from RNAstar GeneCount -> Generate raw count file ## 
for M in [1,3,6,12,18,24]:
    for R in range(1,4,1):
        print("Reading: Sample_M{}_R{}.txt".format(M,R))
        # Initialize file names #
        file_name_in = "RSEM_Sample_M{}_R{}.txt".format(M,R)
        file_name_out = "Sample_M{}_R{}_TPM.csv".format(M,R)
        # Load data #
        df = pd.read_csv(os.path.join(Folder1,file_name_in),sep="\t")
        # Rename columns #
        df = df.rename(columns=({"gene_id":"Ensembl ID"}))
        # Save file #
        df[['Ensembl ID','TPM']].to_csv(os.path.join(Folder2,file_name_out),index=False,sep=";")
        ## Create count table ##
        # If first iteration -> Create Count table #
        if M == 1 and R == 1:
            df_count_table = df[['Ensembl ID','TPM']].rename(columns=({'TPM':'Sample M{} R{} (TPM)'.format(M,R)})).copy()
        # else append samples into count table #
        else:
            df_count_table = pd.merge(df_count_table, df[['Ensembl ID','TPM']].rename(columns=({'TPM':'Sample M{} R{} (TPM)'.format(M,R)})), on="Ensembl ID")

## reduce count table ##

# Remove rows with only zeroes #
df_count_table_reduced = df_count_table.set_index('Ensembl ID').loc[~(df_count_table.set_index('Ensembl ID')==0).all(axis=1)]

## Save count tables to file ##

# Count table #
df_count_table.to_csv(os.path.join(Folder3,File1),index=False,sep=";")

# Redued count table
df_count_table_reduced.to_csv(os.path.join(Folder3,File2),sep=";")