# -*- coding: utf-8 -*-

#### Convert RNA-STAR gene count files to csv ###

## Import libraries ##

import os
import pandas as pd

## Folders #
Folder1 = "Data/Raw data/GeneCounts"
Folder2 = "Data/Raw data/Raw counts"
os.makedirs(Folder2,exist_ok=True)
Folder3 = "Data/Count Tables"
os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "Count table - Readcounts.csv"
File2 = "Reduced - Count table - Readcounts.csv"

## for each sample file from RNAstar GeneCount -> Generate raw count file ## 
for M in [1,3,6,12,18,24]:
    for R in range(1,4,1):
        print("Reading: Sample_M{}_R{}.txt".format(M,R))
        
        # Initialize file names #
        file_name_in = "RawCount_Sample_M{}_R{}.txt".format(M,R)
        file_name_out = "Sample_M{}_R{}_rawcounts.csv".format(M,R)
        # Load data #
        df = pd.read_csv(os.path.join(Folder1,file_name_in),header=None,sep="\t",skiprows=[0,1,2,3])
        # Rename columns #
        df = df.rename(columns=({0:"Ensembl ID",2:"Raw counts"}))
        # Save file #
        df[['Ensembl ID','Raw counts']].to_csv(os.path.join(Folder2,file_name_out),index=False,sep=";")
        ## Create count table ##
        # If first iteration -> Create Count table #
        if M == 1 and R == 1:
            df_count_table = df[['Ensembl ID','Raw counts']].rename(columns=({'Raw counts':'Sample M1 R1'})).copy()
        # else append samples into count table #
        else:
            df_count_table = pd.merge(df_count_table, df[['Ensembl ID','Raw counts']].rename(columns=({'Raw counts':'Sample M{} R{}'.format(M,R)})), on="Ensembl ID")
## reduce count table ##

# Remove rows with only zeroes #
df_count_table_reduced = df_count_table.set_index('Ensembl ID').loc[~(df_count_table.set_index('Ensembl ID')==0).all(axis=1)]

## Save count tables to file ##

# Count table #
df_count_table.to_csv(os.path.join(Folder3,File1),index=False,sep=";")

# Redued count table
df_count_table_reduced.to_csv(os.path.join(Folder3,File2),sep=";")