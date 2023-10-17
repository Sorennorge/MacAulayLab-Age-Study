# -*- coding: utf-8 -*-

### Statistical evaluation of time series ###

import os
import pandas as pd
from scipy import stats
import numpy as np
from scipy.stats import f_oneway
from scipy.stats import tukey_hsd

## Files ##

File1 = "Transport of interest.xlsx"
File2 = "DEseq2 Normalized data.csv"

## Folders ##

Folder1 = "Data/Gene info"
Folder2 = "Data/normalized"


# Transport of interest #
df_Transport_of_interest = pd.read_excel(os.path.join(Folder1,File1))

## Normalized data ##

df_data = pd.read_csv(os.path.join(Folder2,File2),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))

df_Transport_of_interest_data2 = pd.concat([df_Transport_of_interest.set_index('Rat Ensembl ID'),df_data.set_index('Ensembl ID')],join="inner",axis=1)
df_Transport_of_interest_data2 = df_Transport_of_interest_data2.reset_index().set_index(['index','Gene_Alias'])
df_Transport_of_interest_zscored = df_Transport_of_interest_data2.apply(stats.zscore,axis=1)
df_Transport_of_interest_zscored = df_Transport_of_interest_zscored.assign(m=df_Transport_of_interest_zscored[['Sample.M1.R1','Sample.M1.R2','Sample.M1.R3']].mean(axis=1)).sort_values('m').drop('m', axis=1)

#('ENSRNOG00000030019', 'ATP1A1 (NKA)')
#('ENSRNOG00000011648', 'AQP1')
#('ENSRNOG00000014347', 'SLC4A2 (AE2)')
#('ENSRNOG00000010378', 'SLC4A5 (NBCE2)')
#('ENSRNOG00000005957', 'SLC4A7 (NBCn1)')
#('ENSRNOG00000005307', 'SLC4A10 (NCBE)')
#('ENSRNOG00000015971', 'SLC12A2 (NKCC1)')

Anova_dict = {}
turkey_res_list = {}
for index, rows in df_Transport_of_interest_data2.iterrows():
    M1_temp = rows[0:3].values.tolist()
    M3_temp = rows[3:6].values.tolist()
    M6_temp = rows[6:9].values.tolist()
    M12_temp = rows[9:12].values.tolist()
    M18_temp = rows[12:15].values.tolist()
    M24_temp = rows[15:18].values.tolist()
    results = f_oneway(M1_temp,M3_temp,M6_temp,M12_temp,M18_temp,M24_temp)
    Anova_dict[index[1]] = results[1]
    # if interested, pvalue can be shown here by turkey
    res = tukey_hsd(M1_temp,M3_temp,M6_temp,M12_temp,M18_temp,M24_temp)
    turkey_res_list[index[1]] = res