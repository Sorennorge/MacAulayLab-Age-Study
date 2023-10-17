# -*- coding: utf-8 -*-

### Sort out molecular function data ###

# Libaries #

import os
import pandas as pd
import numpy as np
import mygene
mg = mygene.MyGeneInfo()

## Folders ##

Folder1 = "Data/Panther/gene list"
Folder2 = "Data/Panther/NEW/Molecular functions"
Folder3 = "Data/Gene info"

## Files ##

File1 = "LRT_diff_input_list.xlsx" # Needs to be generated from the "LRT_diff_input_list.csv" file from 5.1.1 #
File2 = "Category_list.txt"
File3 = "Uniprot_mapping_TEST.txt"
File4 = "Handle_file_two.xlsx" # Needs to be generated #
File5 = "Handle_two_manual_check.txt" # Needs to be generated #
File6 = "Molecular function overview.xlsx"

## load data ##

df_diff = pd.read_excel(os.path.join(Folder1,File1),usecols=[1,2])
df_diff = df_diff.replace({np.nan: None})


diff_dict = df_diff.set_index('Ensembl ID')['Gene name'].to_dict()

# Check if all keys are genuine #
"""
for key in diff_dict:
    if diff_dict[key] == None:
        print(key,diff_dict[key])
"""

## Load category list for molecular functions ##

df_Category_list = pd.read_csv(os.path.join(Folder2,File2),header=None)
Category_list = df_Category_list[0].tolist()

def unipro_func(x):
    x = x.split("|")[-1].split("=")[-1]
    return x

df_MF = pd.DataFrame(columns=['Gene name','MF'])

for key in Category_list:
    temp_file = "{}.txt".format(key)
    temp_df = pd.read_csv(os.path.join(Folder2,temp_file),sep="\t",header=None,usecols=[0,1],names=["Gene info","Gene name"])
    temp_df['Uniprot'] = temp_df['Gene info'].apply(unipro_func)
    temp_df['MF'] = key
    df_MF = pd.concat([df_MF,temp_df],ignore_index=True)


df_uniprot = pd.read_csv(os.path.join(Folder3,File3),header=0,names=["Ensembl ID","Gene name",'Swiss','TrEMBL'])
df_uniprot = df_uniprot.replace({np.nan: None})

Uniprot_dict = {}
for index,rows in df_uniprot.iterrows():
    if rows['Ensembl ID'] not in Uniprot_dict:
        Uniprot_dict[rows['Ensembl ID']] = {}
        Uniprot_dict[rows['Ensembl ID']]['Gene'] = rows['Gene name']
        Uniprot_dict[rows['Ensembl ID']]['Uniprot'] = []
        if rows['Swiss'] != None:
            Uniprot_dict[rows['Ensembl ID']]['Uniprot'].append(rows['Swiss'])
        if rows['TrEMBL'] != None:
            Uniprot_dict[rows['Ensembl ID']]['Uniprot'].append(rows['TrEMBL'])
    else:
        if rows['Swiss'] != None:
            Uniprot_dict[rows['Ensembl ID']]['Uniprot'].append(rows['Swiss'])
        if rows['TrEMBL'] != None:
            Uniprot_dict[rows['Ensembl ID']]['Uniprot'].append(rows['TrEMBL'])

# Check if all keys and uniprot align # 
"""
# Check list #
for key in Uniprot_dict:
    if Uniprot_dict[key]['Uniprot'] == []:
        print(key,Uniprot_dict[key]['Gene'],Uniprot_dict[key]['Uniprot'])
"""

# If not, assign missing entries #
#ENSRNOG00000009203 RGD1561795 [] -> "D3ZF18"
Uniprot_dict['ENSRNOG00000009203']['Uniprot'].append("D3ZF18")
#ENSRNOG00000009226 Crygn [] -> "D3ZEG1"
Uniprot_dict['ENSRNOG00000009226']['Uniprot'].append("D3ZEG1")
#ENSRNOG00000051669 AABR07040855.1 [] -> "B0BMV9", "A0A8I5ZXW8", "A0A0G2JZI5"
Uniprot_dict['ENSRNOG00000051669']['Uniprot'].append("B0BMV9")
#ENSRNOG00000008423 Gpr22 [] -> "D4A3U0", "A0A5H1ZRU0","A0A1S3GQF3","A0A8C6RTS4","A0A8C6RQZ0"
Uniprot_dict['ENSRNOG00000008423']['Uniprot'].append("D4A3U0")
#ENSRNOG00000003738 Ush2a [] -> "Q8K3K1", "A0A5H1ZRU3", "A0A1S3F1X5"
Uniprot_dict['ENSRNOG00000003738']['Uniprot'].append("Q8K3K1")
#ENSRNOG00000006921 Rbl1 [] -> "D3ZS28"
Uniprot_dict['ENSRNOG00000006921']['Uniprot'].append("D3ZS28")
#ENSRNOG00000015084 Necab2 [] -> "F1LQY6"
Uniprot_dict['ENSRNOG00000015084']['Uniprot'].append("F1LQY6")

# Check list #
for key in Uniprot_dict:
    if Uniprot_dict[key]['Uniprot'] == []:
        print(key,Uniprot_dict[key]['Gene'],Uniprot_dict[key]['Uniprot'])


# add molecular function to mapped keys and uniprot ids
MF_dict = {}
for index,rows in df_MF.iterrows():
    key = rows['Gene name']
    item = rows['MF']
    uniprot = rows['Uniprot']
    if "," in key:
        temp_key = key.split(",")
        for keys in temp_key:
            if keys not in MF_dict:
                MF_dict[keys] = {}
                MF_dict[keys]['MF'] = []
                MF_dict[keys]['Uniprot'] = []
                MF_dict[keys]['MF'].append(item)
                MF_dict[keys]['Uniprot'].append(uniprot)
            else:
                MF_dict[keys]['MF'].append(item)
                MF_dict[keys]['Uniprot'].append(uniprot)
    else:
        if key not in MF_dict:
            MF_dict[key] = {}
            MF_dict[key]['MF'] = []
            MF_dict[key]['Uniprot'] = []
            MF_dict[key]['MF'].append(item)
            MF_dict[key]['Uniprot'].append(uniprot)
        else:
            MF_dict[key]['MF'].append(item)
            MF_dict[key]['Uniprot'].append(uniprot)

## Remove duplicates from MF_dict ##

for key in MF_dict:
    MF_dict[key]['MF'] = list(set(MF_dict[key]['MF']))
    MF_dict[key]['Uniprot'] = list(set(MF_dict[key]['Uniprot']))
        
overview_dict = {}
Handle_dict_one = {}
Handle_dict_two = {}

counter_1 = 0
counter_2 = 0

### Uniprot -> Handle_dict_one and two ###

# If handle file two have been created #

df_handle_two_uniprot = pd.read_excel(os.path.join(Folder3, File4))
df_handle_two_uniprot_dict = df_handle_two_uniprot.set_index('symbol')['query'].to_dict()

for key in diff_dict:
    if diff_dict[key] in df_handle_two_uniprot_dict:
        if df_handle_two_uniprot_dict[diff_dict[key]] in MF_dict[diff_dict[key]]['Uniprot']:
            Uniprot_dict[key]['Uniprot'].append(df_handle_two_uniprot_dict[diff_dict[key]])

# Add the incorrect uniprot in MF_dict #
#ENSRNOG00000023385;Rpl37a;P61515
MF_dict["Rpl37a"]['Uniprot'].append("P61515")
#ENSRNOG00000031716;LOC100910978;D3ZS41
MF_dict["LOC100910978"]['Uniprot'].append("D3ZS41")


## Handle two -> missing uniprot ids - Manual lookup on uniprot #

df_manual_handling = pd.read_csv(os.path.join(Folder3,File5),sep=";",header=None,names=['Ensembl ID','Gene','Uniprot'])
df_manual_handling_dict = df_manual_handling.set_index("Ensembl ID")['Uniprot'].to_dict()

for key in df_manual_handling_dict:
    Uniprot_dict[key]['Uniprot'].append(df_manual_handling_dict[key])


key_test = "ENSRNOG00000004464"
print(key_test)
print(diff_dict[key_test])

for key in diff_dict:
    if diff_dict[key] not in MF_dict:
        if key not in Handle_dict_one:
            Handle_dict_one[key] = {}
            Handle_dict_one[key]['Gene'] = []
            Handle_dict_one[key]['MF'] = []
            Handle_dict_one[key]['Uniprot'] = []
            Handle_dict_one[key]['Gene'].append(diff_dict[key])
        else:
            Handle_dict_one[key]['Gene'].append(diff_dict[key])
    else:
        """
        ## Check if every gene is in correctly assigned to both dicts ##
        if Uniprot_dict[key]['Gene'] == diff_dict[key]:
            counter_1 += 1
        else:
            counter_2 += 1
        """
        if key not in overview_dict:
            if key == key_test:
                print("Catch 1")
            overview_dict[key] = {}
            overview_dict[key]['Gene'] = diff_dict[key]
            overview_dict[key]['MF'] = []
            overview_dict[key]['Uniprot'] = []
            uni_list_1 = MF_dict[diff_dict[key]]['Uniprot']
            uni_list_2 = Uniprot_dict[key]['Uniprot']
            if any(x in uni_list_1 for x in uni_list_2):
                if key == key_test:
                    print("Catch 2")
                for item in MF_dict[diff_dict[key]]['MF']:
                    overview_dict[key]['MF'].append(item)
                for uni_item in list(set(uni_list_1).intersection(uni_list_2)):
                    overview_dict[key]['Uniprot'].append(uni_item)
            else:
                if key == key_test:
                    print("Catch 3")
                    print(uni_list_1)
                    print(uni_list_2)
                if key not in Handle_dict_two:
                    Handle_dict_two[key] = {}
                    Handle_dict_two[key]['Gene'] = []
                    Handle_dict_two[key]['MF'] = []
                    Handle_dict_two[key]['Uniprot'] = []
                    Handle_dict_two[key]['Gene'].append(diff_dict[key])
                    for item in MF_dict[diff_dict[key]]['MF']:
                        Handle_dict_two[key]['MF'].append(item)
                    for item in MF_dict[diff_dict[key]]['Uniprot']:
                        Handle_dict_two[key]['Uniprot'].append(item)
                else:
                    Handle_dict_two[key]['Gene'].append(diff_dict[key])
                    for item in MF_dict[diff_dict[key]]['MF']:
                        Handle_dict_two[key]['MF'].append(item)
                    for item in MF_dict[diff_dict[key]]['Uniprot']:
                        Handle_dict_two[key]['Uniprot'].append(item)

## For handle one #

for key in Handle_dict_one:
    overview_dict[key] = {}
    overview_dict[key]['Gene'] = diff_dict[key]
    overview_dict[key]['MF'] = ['Unclassified']
    overview_dict[key]['Uniprot'] = []

### transform overview dict to dataframe and save dataframe ###

for key in overview_dict:
    if len(overview_dict[key]['MF']) > 1:
        overview_dict[key]['MF'] = "{}".format(";".join(overview_dict[key]['MF']))
    else:
        overview_dict[key]['MF'] = "{}".format(overview_dict[key]['MF'][0])
    if len(overview_dict[key]['Uniprot']) > 1:
        overview_dict[key]['Uniprot'] = "{}".format(";".join(overview_dict[key]['Uniprot']))
    elif len(overview_dict[key]['Uniprot']) < 1:
        overview_dict[key]['Uniprot'] = ''
    else:
        overview_dict[key]['Uniprot'] = "{}".format(overview_dict[key]['Uniprot'][0])


# Create overview #
df_overview = pd.DataFrame.from_dict(overview_dict,orient='index')

# Save overview #
df_overview.to_excel(os.path.join(Folder2,File6))

"""
# For key_test and manual handling of handle two#
for key in list(Handle_dict_two):
    print("{};{};{}".format(key,Handle_dict_two[key]['Gene'][0],Handle_dict_two[key]['Uniprot'][0]))
"""
"""
### Run this one time and save to data ###
### Handle two - print catcher ###
#mg = mygene.MyGeneInfo()
Query_list_uniprot = []
for key in Handle_dict_two:
    Query_list_uniprot.append(Handle_dict_two[key]['Uniprot'][0])
        
df_handle_two = mg.querymany(Query_list_uniprot, scopes="uniprot", fields=["uniprot", "symbol"], species="Rat", as_dataframe=True)
df_handle_two_subset = df_handle_two[df_handle_two['notfound'] != True]
# Save found datapoints #
df_handle_two_subset.to_excel(os.path.join(Folder3, File4))
"""
