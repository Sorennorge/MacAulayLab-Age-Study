# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:59:36 2023

@author: dcs839
"""

### Create barplot for time series ###

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 
from matplotlib.colors import to_rgba
import matplotlib.ticker as ticker

sns.set(font_scale=1.5)
sns.set_style("white")

## Folders ##

Folder1 = "Data/Gene info"
#Folder2 = "Data/Count Tables"
Folder2 = "Data/Normalized"

Folder3 = "Results/Transport/Interest - time series/Barplot"
os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "Transport of interest.xlsx"
#Count_table_file = "Reduced - Count table - TPM.csv"
Count_table_file = "Normalized_TMM.xlsx"
## color codes ##

M1 = "#FF33CC"
M3 = "#00FF00"
M6 = "#FF9900"
M12 = "#00FFFF"
M18 = "#800080"
M24 = "#B2B2B2"

color_list = [M1,M3,M6,M12,M18,M24]

# Convert to RGB colors #

M1_RGB = to_rgba(M1)
M3_RGB = to_rgba(M3)
M6_RGB = to_rgba(M6)
M12_RGB = to_rgba(M12)
M18_RGB = to_rgba(M18)
M24_RGB = to_rgba(M24)

RGB_list = [M1_RGB,M3_RGB,M6_RGB,M12_RGB,M18_RGB,M24_RGB]

## Load data ##

# Transport of interest #
df_Transport_of_interest = pd.read_excel(os.path.join(Folder1,File1))

#ensembl_list = ["ENSRNOG00000030019","ENSRNOG00000015971","ENSRNOG00000014504","ENSRNOG00000013963","ENSRNOG00000031312"]
#gene_alias = ["ATP1A1 (NKA)","SLC12A2 (NKCC1)","IL1R1","IL-6ST","TNFR1"]

#df_Transport_of_interest = pd.DataFrame({'Rat Ensembl ID':ensembl_list,'Gene_Alias':gene_alias})

#ensembl_list = ["ENSRNOG00000016050","ENSRNOG00000016374","ENSRNOG00000017392"]
#gene_alias = ['FGFR1','FGFR2',"FGF2"]

#df_Transport_of_interest = pd.DataFrame({'Rat Ensembl ID':ensembl_list,'Gene_Alias':gene_alias})

#df_data = pd.read_csv(os.path.join(Folder2,Count_table_file),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))
df_data = pd.read_excel(os.path.join(Folder2,Count_table_file))


df_data_interest = pd.concat([df_Transport_of_interest.set_index('Rat Ensembl ID'),df_data.set_index('Ensembl ID')],join="inner",axis=1)


df_data_interest = df_data_interest.set_index('Gene_Alias')
df_data_interest = df_data_interest.astype(float)

Sample_list = ['M1','M3','M6','M12','M18','M24']

for key in Sample_list:
    df_data_interest = df_data_interest.rename(columns=({'Sample {} R1'.format(key):'{}'.format(key),
                                                      'Sample {} R2'.format(key):'{}'.format(key),
                                                      'Sample {} R3'.format(key):'{}'.format(key)}))

gene_list = df_data_interest.index.tolist()

TPM_max_value = df_data_interest.max().max()

i_value = 6
#for i in range(0,7,1):
#i_value = i
bar_data = df_data_interest[df_data_interest.index == gene_list[i_value]]
bar_data_stacked = bar_data.stack().reset_index().rename(columns=({'level_1':"Months",0:'TPM'}))
bar_data_T = bar_data.T

bar_data_stacked['x_value'] = 1

#for key in Sample_list:
plt.figure(figsize=(3,10))
g = sns.pointplot(data=bar_data_stacked,x='x_value',hue='Months',y='TPM',palette=RGB_list,dodge=0.20,scale=0.70,errwidth=2,orient="v")
sns.despine(top=True, right=True, left=False, bottom=True)
plt.xlabel("")
plt.xticks([])
plt.yticks(rotation=90,va="center")
plt.legend([],[], frameon=False)
plt.ylim([0,45000])
plt.xlim([-0.15,0.75])
#plt.title(gene_list[i_value])
#plt.subplots_adjust(bottom=-1)
g.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1000) + 'K'))

File_out = "{}_{}_time_series.png".format(i_value,gene_list[i_value])
#plt.savefig(os.path.join(Folder3,File_out),dpi=600,bbox_inches='tight',transparent=True)
print(gene_list[i_value])

"""
plt.figure(figsize=(2,10))
ax = sns.barplot(data=bar_data_stacked)
plt.title(gene_list[i_value])
plt.ylabel('TPM')
ax.yaxis.set_label_coords(-0.5, 0.5)
plt.ylim([0,600])
File_out = "{}_{}_time_series.png".format(i_value,gene_list[i_value])
#plt.savefig(os.path.join(Folder3,File_out),dpi=600,bbox_inches='tight')
plt.show()
"""