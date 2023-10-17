# -*- coding: utf-8 -*-

### Transport dataframe ###

import os
import pandas as pd

### Folders ###

Folder1 = "Data/Panther/NEW/Protein classes"
Folder2 = "Data/Normalized"

### Files ###

File1 = "Protein classes overview.xlsx"
File2 = "DEseq2 Normalized data.csv"

File_out_1 = "Transporter overview.xlsx"
File_out_2 = "Protein classes overview V2.xlsx"

### Load data ###

df = pd.read_excel(os.path.join(Folder1,File1),usecols=[0,1,2],names=['Ensembl ID','Gene','PC'])

Transporter_dict = {}
PC_dict = {}

for index,rows in df.iterrows():
    ensembl = rows['Ensembl ID']
    genes = rows['Gene']
    protein_class = rows['PC']
    if 'transporter' in protein_class:
        Transporter_dict[ensembl] = genes
    if ensembl not in PC_dict:
        PC_dict[ensembl] = []
        if ";" in protein_class:
            protein_class_temp = protein_class.split(";")
            for key in protein_class_temp:
                PC_dict[ensembl].append(key)
        else:
            PC_dict[ensembl].append(protein_class)
    else:
        if ";" in protein_class:
            protein_class_temp = protein_class.split(";")
            for key in protein_class_temp:
                PC_dict[ensembl].append(key)
        else:
            PC_dict[ensembl].append(protein_class)

df_transport = pd.DataFrame.from_dict(Transporter_dict, orient='index',columns=['Gene'])
df_transport.index = df_transport.index.set_names(['Ensembl ID'])
df_transport = df_transport.reset_index()

for key in PC_dict:
    if len(PC_dict[key]) > 1:
        if "Unclassified" in PC_dict[key]:
            PC_dict[key].remove('Unclassified')
        else:
            pass
df_PC = pd.DataFrame.from_dict(PC_dict, orient='index',columns=['PC_1','PC_2'])
df_PC.index = df_PC.index.set_names(['Ensembl ID'])
df_PC = df_PC.reset_index()

## Add Deseq2 normatlized to transport dataframe ##

df_data = pd.read_csv(os.path.join(Folder2,File2),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))

## modulate data ##

df_data['M1'] = df_data[['Sample.M1.R1','Sample.M1.R2','Sample.M1.R3']].mean(axis=1)
df_data['M3'] = df_data[['Sample.M3.R1','Sample.M3.R2','Sample.M3.R3']].mean(axis=1)
df_data['M6'] = df_data[['Sample.M6.R1','Sample.M6.R2','Sample.M6.R3']].mean(axis=1)
df_data['M12'] = df_data[['Sample.M12.R1','Sample.M12.R2','Sample.M12.R3']].mean(axis=1)
df_data['M18'] = df_data[['Sample.M18.R1','Sample.M18.R2','Sample.M18.R3']].mean(axis=1)
df_data['M24'] = df_data[['Sample.M24.R1','Sample.M24.R2','Sample.M24.R3']].mean(axis=1)                 

df_data_mean = df_data[['Ensembl ID',"M1","M3","M6","M12","M18","M24"]]

df_data_transformed = df_data_mean[df_data_mean['Ensembl ID'].isin(df_transport['Ensembl ID'])].set_index('Ensembl ID')
df_data_transformed_norm_row = df_data_transformed.apply(lambda x: (x-x.mean())/x.std(), axis = 1)

df_transport = pd.concat([df_transport.set_index('Ensembl ID'),df_data_transformed_norm_row],join='inner',axis=1)

# Save dataframes #

df_transport.to_excel(os.path.join(Folder1,File_out_1))
df_PC.to_excel(os.path.join(Folder1,File_out_2),index=False)