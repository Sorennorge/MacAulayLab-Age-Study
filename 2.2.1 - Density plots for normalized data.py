# -*- coding: utf-8 -*-

### Density plot of normalized data ##

## Libraries ##

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from natsort import natsort_keygen

## Folders ##

Folder1 = "Data/Normalized"
Folder2 = "Results/Density plots"
os.makedirs(Folder2,exist_ok=True)
Folder3 = "Results/Density plots/Genes of interest"
os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "DEseq2 Normalized data.csv"
File2 = "DEseq2 rlog transformed.csv"
File3 = "Normalized_TMM.csv"
File4 = "Normalized_geTMM.csv"
File5 = "Normalized_MRN.csv"

File_out_1 = "Density_plot_DEseq2.png"
File_out_2 = "Density_plot_DEseq2_zoomed.png"
File_out_genes = "Gene_list_of_interest.xlsx"

## Load data ##

df_norm = pd.read_csv(os.path.join(Folder1,File1),sep=";",decimal=",").rename(columns={"Unnamed: 0":"Ensembl ID"}).set_index("Ensembl ID")

# For interest, density can be created using different normalized data
#df_rlog = pd.read_csv(os.path.join(Folder1,File2),sep=";",decimal=",").rename(columns={"Unnamed: 0":"Ensembl ID"}).set_index("Ensembl ID")
#df_TMM = pd.read_csv(os.path.join(Folder1,File3),sep=";").set_index("Ensembl ID")
#df_geTMM = pd.read_csv(os.path.join(Folder1,File4),sep=";").set_index("Ensembl ID")
#df_MRN = pd.read_csv(os.path.join(Folder1,File5),sep=";").set_index("Ensembl ID")


## Color scheme ##

M1 = "#07F2F2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#A422E6"
M24 = "#B3B3B3"

colorpallette = [M1,M3,M6,M12,M18,M24]

### DEseq2 Norm ###

## Average month within dataframes ##

for i in [1,3,6,12,18,24]:
    df_norm['M{}'.format(i)] =  df_norm[['Sample.M{}.R1'.format(i),
                                    'Sample.M{}.R2'.format(i),
                                    'Sample.M{}.R3'.format(i)]].mean(axis=1)

## Create mean dataframe ##
df_norm_mean = df_norm[['M1','M3','M6','M12','M18','M24']]

## Stack dataframe with mean values ##
df_norm_stacked = df_norm_mean.stack().reset_index().sort_values(by='level_1',key=natsort_keygen()).reset_index(drop=True)
df_norm_stacked2 = df_norm_mean.stack().reset_index().sort_values(by='level_1',key=natsort_keygen()).reset_index(drop=False)

df_norm_stacked = df_norm_stacked.rename(columns={"level_1":'Month',0:"MRN"}).drop(columns=['Ensembl ID'])
df_norm_stacked2 = df_norm_stacked2.rename(columns={"level_1":'Month',0:"MRN"}).drop(columns=['index'])
## remove all zero values ##
df_norm_stacked = df_norm_stacked[df_norm_stacked['MRN'] > 0]
df_norm_stacked2 = df_norm_stacked2[df_norm_stacked2['MRN'] > 0]
# add log10 to dataframe #
df_norm_stacked['log10(MRN)'] = np.log10(df_norm_stacked['MRN'])
df_norm_stacked2['log10(MRN)'] = np.log10(df_norm_stacked2['MRN'])

hue_order_list = ['M24','M3','M6','M12','M18','M1']

df_range_of_interest = df_norm_stacked2[(df_norm_stacked2['log10(MRN)'] >= -0.5) & (df_norm_stacked2['log10(MRN)'] <= 2)]


# Create figure #

fig, ax = plt.subplots(figsize=(12,10),sharex=True)

sns.kdeplot(data=df_norm_stacked, ax=ax, x="log10(MRN)",hue="Month",fill=False,legend=False,lw=2,palette=colorpallette)
sns.despine(right=True,top=True,left=False,bottom=False,ax=ax)
ax.set_ylim([0,0.07])
plt.hlines(y=[0.024,0.042], xmin=-0.5, xmax=2, colors='black', ls='dashed', lw=2,alpha=0.5)
plt.vlines(x=[-0.5,2], ymin=0.024, ymax=0.042, colors='black', ls='dashed', lw=2,alpha=0.5)
# Add dotted lines between plots #
plt.plot([2.05,3.95], [0.0241,0.03],linestyle='--',lw=2,color='black',alpha=0.5)
plt.plot([2.05,3.95], [0.042,0.067],linestyle='--',lw=2,color='black',alpha=0.5)
# set labels for main plot #
plt.xlabel("Gene expression ($log_{10}(MRN)$)", size=30)
ax.xaxis.labelpad = 20
plt.ylabel("Gene density", size=30)
ax.yaxis.labelpad = 20
# set label ticks font size #

ax.tick_params(axis='both', which='major', labelsize=25)

# plot zoomed #
ax2 = plt.axes([0.60, 0.45, .4, .4])
sns.kdeplot(data=df_norm_stacked, ax=ax2, x="log10(MRN)",hue="Month",fill=False,lw=2,legend=False,palette=colorpallette)
sns.despine(top=True,bottom=True,right=True,left=True,ax=ax2)
ax2.set_ylim([0.024,0.042])
ax2.set_xlim([-0.5,2])
ax2.tick_params(bottom=False,left=False)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_ylabel('')
ax2.set_xlabel('')
plt.hlines(y=[0.024,0.042], xmin=-0.5, xmax=2, colors='black', ls='dashed', lw=3.2,alpha=0.5)
plt.vlines(x=[-0.5,2], ymin=0.024, ymax=0.042, colors='black', ls='dashed', lw=3.2,alpha=0.5)

plt.savefig(os.path.join(Folder2,File_out_1),dpi=600,bbox_inches='tight')
#plt.show()


## Plotting for different normalization ##
"""
## Density plot ##
plt.figure(figsize=(12,4))
sns.displot(df_norm_stacked, x="log10(MRN)",hue="Month", kind="kde",fill=False,legend=False,palette=colorpallette)
plt.vlines(x=[-0.5,2], ymin=0.024, ymax=0.042, colors='black', ls='dashed', lw=2,alpha=0.5)
plt.hlines(y=[0.024,0.042], xmin=-0.5, xmax=2, colors='black', ls='dashed', lw=2,alpha=0.5)
plt.xlabel("Gene expression ($log_{10}(MRN)$)", size=20)
plt.ylabel("Gene density", size=20)
plt.xticks([-2,0,2,4,6])
plt.ylim([0,0.08])
plt.tick_params(bottom=True)
plt.tick_params(left=True)
#plt.title("Age study - Gene expression density", size=20)
plt.savefig(os.path.join(Folder2,File_out_1),dpi=600,bbox_inches='tight')
plt.show()


### Zoomed ###
figure, ax = plt.subplots(figsize=(10,8))
ax = sns.displot(df_norm_stacked, x="log10(MRN)",hue="Month", kind="kde",fill=False,legend=False,palette=colorpallette)
plt.ylim([0.024,0.042])
plt.xlim([-0.5,2])
sns.despine(top=True, right=True, left=True, bottom=True)
plt.vlines(x=[-0.5,2], ymin=0.024, ymax=0.042, colors='black', ls='dashed', lw=6,alpha=0.5)
plt.hlines(y=[0.024,0.042], xmin=-0.5, xmax=2, colors='black', ls='dashed', lw=6,alpha=0.5)
plt.xlabel("")
plt.ylabel("")
plt.xticks([])
plt.yticks([])
#plt.xlabel("Gene expression $log_{10}(MRN)$", size=20)
#plt.ylabel("Density", size=20)
#plt.title("Age study - Gene expression density", size=20)
plt.savefig(os.path.join(Folder2,File_out_2),dpi=600,bbox_inches='tight')
#plt.show()

### MRN ###

## Average month within dataframes ##
for i in [1,3,6,12,18,24]:
    df_MRN['M{}'.format(i)] =  df_MRN[['Sample M{} R1'.format(i),
                                    'Sample M{} R2'.format(i),
                                    'Sample M{} R3'.format(i)]].mean(axis=1)

## Create mean dataframe ##
df_MRN_mean = df_MRN[['M1','M3','M6','M12','M18','M24']]

## Stack dataframe with mean values ##
df_MRN_stacked = df_MRN_mean.stack().reset_index().sort_values(by='level_1',key=natsort_keygen()).reset_index(drop=True)
df_MRN_stacked = df_MRN_stacked.rename(columns={"level_1":'Month',0:"MRN"}).drop(columns=['Ensembl ID'])
## remove all zero values ##
df_MRN_stacked = df_MRN_stacked[df_MRN_stacked['MRN'] > 0]
# add log10 to dataframe #
df_MRN_stacked['log10(MRN)'] = np.log10(df_MRN_stacked['MRN'])

## Density plot ##
plt.figure(figsize=(10,8))
sns.displot(df_MRN_stacked, x="log10(MRN)",hue="Month", kind="kde",fill=False)
plt.xlabel("Expression gene log10(MRN)", size=20)
plt.ylabel("Density", size=14)
plt.title("Age study distribution (MRN)", size=20)
plt.savefig(os.path.join(Folder2,File_out_2),dpi=600,bbox_inches='tight')
plt.show()


### TMM ###

## Average month within dataframes ##
for i in [1,3,6,12,18,24]:
    df_TMM['M{}'.format(i)] =  df_TMM[['Sample M{} R1'.format(i),
                                    'Sample M{} R2'.format(i),
                                    'Sample M{} R3'.format(i)]].mean(axis=1)

## Create mean dataframe ##
df_TMM_mean = df_TMM[['M1','M3','M6','M12','M18','M24']]

## Stack dataframe with mean values ##
df_TMM_stacked = df_TMM_mean.stack().reset_index().sort_values(by='level_1',key=natsort_keygen()).reset_index(drop=True)
df_TMM_stacked = df_TMM_stacked.rename(columns={"level_1":'Month',0:"TMM"}).drop(columns=['Ensembl ID'])
## remove all zero values ##
df_TMM_stacked = df_TMM_stacked[df_TMM_stacked['TMM'] > 0]
# add log10 to dataframe #
df_TMM_stacked['log10(TMM)'] = np.log10(df_TMM_stacked['TMM'])

## Density plot ##
plt.figure(figsize=(10,8))
sns.displot(df_TMM_stacked, x="log10(TMM)",hue="Month", kind="kde",fill=False)
plt.xlabel("Expression gene log10(TMM)", size=20)
plt.ylabel("Density", size=14)
plt.title("Age study distribution (TMM)", size=20)
plt.savefig(os.path.join(Folder2,File_out_3),dpi=600,bbox_inches='tight')
plt.show()


### geTMM ###

## Average month within dataframes ##
for i in [1,3,6,12,18,24]:
    df_geTMM['M{}'.format(i)] =  df_geTMM[['Sample M{} R1'.format(i),
                                    'Sample M{} R2'.format(i),
                                    'Sample M{} R3'.format(i)]].mean(axis=1)

## Create mean dataframe ##
df_geTMM_mean = df_geTMM[['M1','M3','M6','M12','M18','M24']]

## Stack dataframe with mean values ##
df_geTMM_stacked = df_geTMM_mean.stack().reset_index().sort_values(by='level_1',key=natsort_keygen()).reset_index(drop=True)
df_geTMM_stacked = df_geTMM_stacked.rename(columns={"level_1":'Month',0:"geTMM"}).drop(columns=['Ensembl ID'])
## remove all zero values ##
df_geTMM_stacked = df_geTMM_stacked[df_geTMM_stacked['geTMM'] > 0]
# add log10 to dataframe #
df_geTMM_stacked['log10(geTMM)'] = np.log10(df_geTMM_stacked['geTMM'])

## Density plot ##
plt.figure(figsize=(10,8))
sns.displot(df_geTMM_stacked, x="log10(geTMM)",hue="Month", kind="kde",fill=False)
plt.xlabel("Expression gene log10(geTMM)", size=20)
plt.ylabel("Density", size=14)
plt.title("Age study distribution (geTMM)", size=20)
plt.savefig(os.path.join(Folder2,File_out_4),dpi=600,bbox_inches='tight')
plt.show()
"""