# -*- coding: utf-8 -*-

### Enrichment analysis ###

## Clean enrichment analysis data ##

## Libraries ##

import os
import pandas as pd

## Function ##

def capitalize(x):
    return_sentence = x
    if "_" in return_sentence:
        return_sentence = return_sentence.replace("_","/")
    else:
        pass
    return_sentence = return_sentence[0].upper()+return_sentence[1:]
    return return_sentence

def percentage(x,y):
    percentage = round((x/y)*100,1)
    return percentage
    

## Folders ##

Folder1 = "Data/Panther/Molecular functions"
Folder2 = "Data/Panther/LRT diff/Enrichment data"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "Molecular function overview.xlsx"
File2 = "MF_enrichment_data_all.xlsx"

## Load data ##

df_data = pd.read_excel(os.path.join(Folder1,File1),names=['Ensembl ID','Gene','MF','Uniprot'])

MF_list = df_data['MF'].tolist()

MF_dict = {}

# Verify that 'Unclassified' is not assigned entries with MF #
for key in MF_list:
    if ";" in key:
        temp_list = key.split(";")
        if "Unclassified" in temp_list:
            temp_list.remove("Unclassified")
            if len(temp_list) == 0:
                print(key)
            elif len(temp_list) > 1:
                for keys in temp_list:
                    keys = capitalize(keys)
                    if keys not in MF_dict:
                        MF_dict[keys] = {}
                        MF_dict[keys]['Count'] = 1
                    else:
                        MF_dict[keys]['Count'] += 1
            else:
                key = capitalize(temp_list[0])
                if key not in MF_dict:
                    MF_dict[key] = {}
                    MF_dict[key]['Count'] = 1
                else:
                    MF_dict[key]['Count'] += 1
        else:
            for keys in temp_list:
                keys = capitalize(keys)
                if keys not in MF_dict:
                    MF_dict[keys] = {}
                    MF_dict[keys]['Count'] = 1
                else:
                    MF_dict[keys]['Count'] += 1
    else:
        key = capitalize(key)
        if key not in MF_dict:
            MF_dict[key] = {}
            MF_dict[key]['Count'] = 1
        else:
            MF_dict[key]['Count'] += 1
sum_of_all = 0
for key in MF_dict:
    sum_of_all += MF_dict[key]['Count']

# Cal percentage #
for key in MF_dict:
    MF_dict[key]['Percentage'] = percentage(MF_dict[key]['Count'], sum_of_all)

overall_percentage = 0
for key in MF_dict:
    overall_percentage += MF_dict[key]['Percentage']

# Assign percentage 'Unclassified #
MF_dict["Unclassified"]['Percentage'] += round(100-overall_percentage,1)

MF_dict["Unclassified"]['Percentage'] = round(MF_dict["Unclassified"]['Percentage'],1)

overall_percentage = 0
for key in MF_dict:
    overall_percentage += MF_dict[key]['Percentage']

# Create dataframe and resort #
df = pd.DataFrame().from_dict(MF_dict,"index").reset_index().rename(columns={"index":"Class"})

df_Enrichment_all_sorted = df.sort_values(by=['Count'],ascending=False,ignore_index=True)

idx = df_Enrichment_all_sorted.index.tolist()
idx.pop(0)
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reindex(idx+[0])
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reset_index(drop=True)

## Save data to files ##

df_Enrichment_all_sorted.to_excel(os.path.join(Folder2,File2),index=False)
