# -*- coding: utf-8 -*-

## The other sheets ##

import os
import pandas as pd
import numpy as np

## Functions ##
def capitalize(x):
    return_sentence = x
    return_sentence = return_sentence[0].upper()+return_sentence[1:]
    return return_sentence

## Folders ##

Folder_1 = "Results/Supplementary tables"
Folder_2 = "Data/Panther/NEW/Molecular functions"
Folder_3 = "Data/Panther/NEW/Protein classes"
Folder_4 = "Data/Panther/NEW/Biological functions"

## Files ##

File_1 = "Supplementary table 1.xlsx"
File_1_2 = "Supplementary differentially expressed table.xlsx"
File_2 = "Molecular function time series overview.xlsx"
File_3 = "Protein classes overview V2.xlsx"
File_4 = "Biological functions overview.xlsx"

File_out_1 = "Supplementary table 1 - Sheet 2.xlsx"

## Load data ##
"""
# If differential expressed table have not yet been generated #
df_overview = pd.read_excel(os.path.join(Folder_1,File_1))
df_overview_diff = df_overview[(df_overview['Padj'] == "< 0.001") | (df_overview['Padj'] == "< 0.01") | (df_overview['Padj'] == "< 0.05")]
df_overview_diff = df_overview_diff[['Ensembl ID','Gene name']].set_index('Ensembl ID')
df_overview_diff.to_excel(os.path.join(Folder_1,File_1_2))
"""
# If table is generated #
df_overview_diff = pd.read_excel(os.path.join(Folder_1,File_1_2)).set_index("Ensembl ID")

## load MF, PC, and BF ##

# MF #
df_MF = pd.read_excel(os.path.join(Folder_2,File_2))

MF_dict = df_MF.groupby('Ensembl ID')['MF'].agg(list).to_dict()

# PC #

df_PC = pd.read_excel(os.path.join(Folder_3,File_3))

df_PC = df_PC.set_index('Ensembl ID')
PC_dict = df_PC.T.to_dict(orient='list')

cleaned_PC_dict = {}
for key in PC_dict:
    items_cleaned = [x for x in PC_dict[key] if str(x) != 'nan']
    items_cleaned_capped = [capitalize(x) for x in items_cleaned]
    cleaned_PC_dict[key] = items_cleaned_capped
    
# BF #

df_BF = pd.read_excel(os.path.join(Folder_4,File_4)).rename(columns=({'Unnamed: 0':'Ensembl ID'}))
df_BF = df_BF[['Ensembl ID','BF']]
BF_dict = df_BF.set_index('Ensembl ID')['BF'].to_dict()

cleaned_BF_dict = {}

for key in BF_dict:
    temp_list = BF_dict[key].split(";")
    cap_temp_list = [capitalize(x) for x in temp_list]
    if len(cap_temp_list) > 1:
        if 'Unclassified' in cap_temp_list:
            cap_temp_list.remove('Unclassified')
            cleaned_BF_dict[key] = cap_temp_list
        else:
            cleaned_BF_dict[key] = cap_temp_list
    else:
        cleaned_BF_dict[key] = cap_temp_list
        
### Create dicts for dataframe ###

pre_dataframe_dict_MF = {}
pre_dataframe_dict_PC = {}
pre_dataframe_dict_BF = {}

for key in MF_dict:
    if len(MF_dict[key]) > 1:
        sorted_list = sorted(MF_dict[key])
        temp_value = "{}".format("|".join(sorted_list))
        pre_dataframe_dict_MF[key] = temp_value
    else:
        pre_dataframe_dict_MF[key] = MF_dict[key][0]
        
for key in cleaned_PC_dict:
    if len(cleaned_PC_dict[key]) > 1:
        sorted_list = sorted(cleaned_PC_dict[key])
        temp_value = "{}".format("|".join(sorted_list))
        pre_dataframe_dict_PC[key] = temp_value
    else:
        pre_dataframe_dict_PC[key] = cleaned_PC_dict[key][0]
        
for key in cleaned_BF_dict:
    if len(cleaned_BF_dict[key]) > 1:
        sorted_list = sorted(cleaned_BF_dict[key])
        temp_value = "{}".format("|".join(sorted_list))
        pre_dataframe_dict_BF[key] = temp_value
    else:
        pre_dataframe_dict_BF[key] = cleaned_BF_dict[key][0]

## Convert dicts to dataframe ##

pre_dataframe_MF = pd.DataFrame.from_dict(pre_dataframe_dict_MF, orient='index').rename(columns={0:"MF"})
pre_dataframe_PC = pd.DataFrame.from_dict(pre_dataframe_dict_PC, orient='index').rename(columns={0:"PC"})
pre_dataframe_BF = pd.DataFrame.from_dict(pre_dataframe_dict_BF, orient='index').rename(columns={0:"BF"})

## Generate dataframe for sheet 2 ##

df_overview_sheet2 = pd.concat([df_overview_diff,pre_dataframe_MF,pre_dataframe_PC,pre_dataframe_BF],join='inner',axis=1)
df_overview_sheet2 = df_overview_sheet2.reset_index().rename(columns={'index':'Ensembl ID'})
df_overview_sheet2.to_excel(os.path.join(Folder_1,File_out_1),index=False)