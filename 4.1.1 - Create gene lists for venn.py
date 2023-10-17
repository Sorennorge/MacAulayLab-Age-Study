# -*- coding: utf-8 -*-

### Create gene lists for venn ###

import os
import pandas as pd


## Input parameters ##

TMM = 10

## Folders ##

Folder1 = "Data/Normalized"
Folder2 = "Data/Venn data"

os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "Normalized_TMM.xlsx"

## Load data ##

df = pd.read_excel(os.path.join(Folder1,File1))

for i in [1,3,6,12,18,24]:
    df['M{}_mean'.format(i)] =  df[['Sample M{} R1'.format(i),
                                    'Sample M{} R2'.format(i),
                                   'Sample M{} R3'.format(i)]].mean(axis=1)
df_mean = df2[['Ensembl ID','M1_mean','M3_mean','M6_mean','M12_mean','M18_mean','M24_mean']].set_index('Ensembl ID')


month_1 = df_mean.loc[df_mean['M1_mean'] > TMM].index.to_frame()
month_3 = df_mean.loc[df_mean['M3_mean'] > TMM].index.to_frame()
month_6 = df_mean.loc[df_mean['M6_mean'] > TMM].index.to_frame()
month_12 = df_mean.loc[df_mean['M12_mean'] > TMM].index.to_frame()
month_18 = df_mean.loc[df_mean['M18_mean'] > TMM].index.to_frame()
month_24 = df_mean.loc[df_mean['M24_mean'] > TMM].index.to_frame()

# Cal overlap #
Shared_all = len(list(set(month_1.index) 
               & set(month_3.index) 
               & set(month_6.index) 
               & set(month_12.index) 
               & set(month_18.index) 
               & set(month_24.index)))
# 1 month #
Exclusive_1 = len(list(set(month_1.index) 
               - set(month_3.index) 
               - set(month_6.index) 
               - set(month_12.index) 
               - set(month_18.index) 
               - set(month_24.index)))

# 3 month #
Exclusive_3 = len(list(set(month_3.index) 
               - set(month_1.index) 
               - set(month_6.index) 
               - set(month_12.index) 
               - set(month_18.index) 
               - set(month_24.index)))

# 6 month #
Exclusive_6 = len(list(set(month_6.index) 
               - set(month_1.index) 
               - set(month_3.index) 
               - set(month_12.index) 
               - set(month_18.index) 
               - set(month_24.index)))

# 12 month #
Exclusive_12 = len(list(set(month_12.index) 
               - set(month_1.index) 
               - set(month_3.index) 
               - set(month_6.index) 
               - set(month_18.index) 
               - set(month_24.index)))

# 18 month #
Exclusive_18 = len(list(set(month_18.index) 
               - set(month_1.index) 
               - set(month_3.index) 
               - set(month_6.index) 
               - set(month_12.index) 
               - set(month_24.index)))

# 24 month #
Exclusive_24 = (len(list(set(month_24.index) 
               - set(month_1.index) 
               - set(month_3.index) 
               - set(month_6.index) 
               - set(month_12.index) 
               - set(month_18.index))))


Exclusive_1_P = Exclusive_1/len(month_1.index)*100
Exclusive_3_P = Exclusive_3/len(month_3.index)*100
Exclusive_6_P = Exclusive_6/len(month_6.index)*100
Exclusive_12_P = Exclusive_12/len(month_12.index)*100
Exclusive_18_P = Exclusive_18/len(month_18.index)*100
Exclusive_24_P = Exclusive_24/len(month_24.index)*100

print("TMM: {}".format(TMM))
print("Number of genes included:")
print("1M: {} - 3M: {} - 6M: {} - 12M: {} - 18M: {} - 24M: {}".format(len(month_1.index),
                                                                      len(month_3.index),
                                                                      len(month_6.index),
                                                                      len(month_12.index),
                                                                      len(month_18.index),
                                                                      len(month_24.index)))
print("shared all: {} \nExclusive 1M: {} ({:.2f} %)\nExclusive 3M: {} ({:.2f} %)\nExclusive 6M: {} ({:.2f} %)\nExclusive 12M: {} ({:.2f} %)\nExclusive 18M: {} ({:.2f} %)\nExclusive 24M: {} ({:.2f} %)".format(
    Shared_all,Exclusive_1,Exclusive_1_P,
    Exclusive_3,Exclusive_3_P,
    Exclusive_6,Exclusive_6_P,
    Exclusive_12,Exclusive_12_P,
    Exclusive_18,Exclusive_18_P,
    Exclusive_24,Exclusive_24_P))

# Cal percentage #
all_1_P = Shared_all/len(month_1.index)*100
all_3_P = Shared_all/len(month_3.index)*100
all_6_P = Shared_all/len(month_6.index)*100
all_12_P = Shared_all/len(month_12.index)*100
all_18_P = Shared_all/len(month_18.index)*100
all_24_P = Shared_all/len(month_24.index)*100

print('\nShared all perentage:')
print("1M: {:.2f}\n3M: {:.2f}\n6M: {:.2f}\n12M: {:.2f}\n18M: {:.2f}\n24M: {:.2f}".format(
    all_1_P,all_3_P,all_6_P,all_12_P,all_18_P,all_24_P))

# Save data to files #
month_1.to_csv(os.path.join(Folder2,"Month_1.txt"),index=False)
month_3.to_csv(os.path.join(Folder2,"Month_3.txt"),index=False)
month_6.to_csv(os.path.join(Folder2,"Month_6.txt"),index=False)
month_12.to_csv(os.path.join(Folder2,"Month_12.txt"),index=False)
month_18.to_csv(os.path.join(Folder2,"Month_18.txt"),index=False)
month_24.to_csv(os.path.join(Folder2,"Month_24.txt"),index=False)