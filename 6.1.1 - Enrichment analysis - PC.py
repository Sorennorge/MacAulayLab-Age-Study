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

Folder1 = "Data/Panther/Protein classes"
Folder2 = "Data/Panther/LRT diff/Enrichment data"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "Protein classes overview.xlsx"
File2 = "PC_enrichment_data_all.xlsx"

## Load data ##

df_data = pd.read_excel(os.path.join(Folder1,File1),names=['Ensembl ID','Gene','PC','Uniprot'])

PC_list = df_data['PC'].tolist()

PC_dict = {}

# Verify that 'Unclassified' is not assigned entries with PC #
for key in PC_list:
    if ";" in key:
        temp_list = key.split(";")
        if "Unclassified" in temp_list:
            temp_list.remove("Unclassified")
            if len(temp_list) == 0:
                print(key)
            elif len(temp_list) > 1:
                for keys in temp_list:
                    keys = capitalize(keys)
                    if keys not in PC_dict:
                        PC_dict[keys] = {}
                        PC_dict[keys]['Count'] = 1
                    else:
                        PC_dict[keys]['Count'] += 1
            else:
                key = capitalize(temp_list[0])
                if key not in PC_dict:
                    PC_dict[key] = {}
                    PC_dict[key]['Count'] = 1
                else:
                    PC_dict[key]['Count'] += 1
        else:
            for keys in temp_list:
                keys = capitalize(keys)
                if keys not in PC_dict:
                    PC_dict[keys] = {}
                    PC_dict[keys]['Count'] = 1
                else:
                    PC_dict[keys]['Count'] += 1
    else:
        key = capitalize(key)
        if key not in PC_dict:
            PC_dict[key] = {}
            PC_dict[key]['Count'] = 1
        else:
            PC_dict[key]['Count'] += 1
sum_of_all = 0
for key in PC_dict:
    sum_of_all += PC_dict[key]['Count']

for key in PC_dict:
    PC_dict[key]['Percentage'] = percentage(PC_dict[key]['Count'], sum_of_all)

# Cal percentage #
overall_percentage = 0
for key in PC_dict:
    overall_percentage += PC_dict[key]['Percentage']

# Assign percentage 'Unclassified #
PC_dict["Unclassified"]['Percentage'] += round(100-overall_percentage,1)  

# Create dataframe and resort #
df = pd.DataFrame().from_dict(PC_dict,"index").reset_index().rename(columns={"index":"Class"})

df_Enrichment_all_sorted = df.sort_values(by=['Count'],ascending=False,ignore_index=True)

idx = df_Enrichment_all_sorted.index.tolist()
idx.pop(0)
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reindex(idx+[0])
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reset_index(drop=True)

## Save data to files ##

df_Enrichment_all_sorted.to_excel(os.path.join(Folder2,File2),index=False)