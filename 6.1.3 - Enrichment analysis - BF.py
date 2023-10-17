# -*- coding: utf-8 -*-

### Enrichment analysis ###

## Clean enrichment analysis data ##

## Libraries ##

import os
import pandas as pd

## Function ##

def remove_go_term_and_capitalize(x):
    temp = x.split(" (")
    return_sentence = temp[0]
    if return_sentence == "No PANTHER category is assigned":
        return_sentence = 'Unclassified'
    else:
        pass
    return_sentence = return_sentence[0].upper()+return_sentence[1:]
    return return_sentence

## Folders ##

Folder = "Data/Panther/LRT diff"
Folder2 = "Data/Panther/LRT diff/Enrichment data"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "BF.txt"
File2 = "BF_enrichment_data_all.xlsx"
File3 = "BF_enrichment_data_classified.xlsx"

## Load data ##

df = pd.read_csv(os.path.join(Folder,File1),sep="\t",header=None,usecols=[1,2,4],names=['Class','Count','percentage all'])

df['Class'] = df['Class'].apply(remove_go_term_and_capitalize)

df_Enrichment_all_sorted = df.sort_values(by=['Count'],ascending=False,ignore_index=True)

idx = df_Enrichment_all_sorted.index.tolist()
idx.pop(1)
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reindex(idx+[1])
df_Enrichment_all_sorted = df_Enrichment_all_sorted.reset_index(drop=True)

df_Enrichment_all_classified_sorted = df_Enrichment_all_sorted.copy()
df_Enrichment_all_classified_sorted.drop(df_Enrichment_all_classified_sorted.tail(1).index,inplace=True) # drop last n rows

## Save data to files ##

df_Enrichment_all_sorted.to_excel(os.path.join(Folder2,File2),index=False)
df_Enrichment_all_classified_sorted.to_excel(os.path.join(Folder2,File3),index=False)
