# -*- coding: utf-8 -*-

### Generate Molecular function time series data ###

## Libraries ##

import os
import pandas as pd

## Functions ##

def capitalize(x):
    return_sentence = x
    if "_" in return_sentence:
        return_sentence = return_sentence.replace("_","/")
    else:
        pass
    return_sentence = return_sentence[0].upper()+return_sentence[1:]
    return return_sentence

## Folders ##

Folder1 = "Data/Panther/Molecular functions"
Folder2 = "Data/Panther/LRT diff/Enrichment data"

## Files ##

File1 = "Molecular function overview.xlsx"
File2 = "MF_enrichment_data_all.xlsx"
File3 = "Molecular function time series overview.xlsx"

## Load data ##

df_overview = pd.read_excel(os.path.join(Folder1,File1),names=['Ensembl ID',"Gene name", "MF", "Uniprot"])
df_enrichment = pd.read_excel(os.path.join(Folder2,File2))

Molecular_func_list = df_enrichment['Class'].tolist()
Molecular_func_check_list = Molecular_func_list[:-1]

## modify data ##

Sorted_dict = {}
for index,rows in df_overview.iterrows():
    key = rows['Ensembl ID']
    item = rows['MF']
    if key not in Sorted_dict:
        Sorted_dict[key] = []
        if ";" in item:
            items = item.split(";")
            for x in items:
                x = capitalize(x)
                Sorted_dict[key].append(x)
        else:
            Sorted_dict[key].append(capitalize(item))


for key in Sorted_dict:
    if len(Sorted_dict[key]) > 1:
        if Molecular_func_list[-1] in Sorted_dict[key]:
            Sorted_dict[key].remove(Molecular_func_list[-1])

df_overview_modified = pd.DataFrame.from_dict(Sorted_dict,"index").stack().reset_index()
df_overview_modified = df_overview_modified.rename(columns={'level_0':"Ensembl ID",0:"MF"})[["Ensembl ID","MF"]]

df_overview_modified.to_excel(os.path.join(Folder1,File3),index=False)
            
        