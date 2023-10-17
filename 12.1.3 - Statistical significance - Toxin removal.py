# -*- coding: utf-8 -*-

### P_values for Toxin removal ###

import pandas as pd
import os
from scipy import stats
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


## Folders ##

Folder_1 = "Data/Pvalue data/Toxin removal"

# Files ##

File_1 = "File_list.xlsx"

# Load file data #

File_df = pd.read_excel(os.path.join(Folder_1,File_1))
File_list = File_df['File names'].tolist()

# Categories #

## Global variables ##

Anova_dict = {}

for File in File_list:
#File = File_list[0]
    Category_name = File.replace(".xlsx",'')
    File_out_1 = "Tukey data - {}.xlsx".format(Category_name)
    File_out_2 = "Tukey data Sinificant - {}.xlsx".format(Category_name)
    
    temp_df_stacked = pd.read_excel(os.path.join(Folder_1,File)).set_index(["Ensembl ID",'Months'])
    
    data_unstacked = temp_df_stacked.unstack(level=-1)
    data_unstacked = data_unstacked.droplevel(None, axis=1)
    
    
    Anova = f_oneway(data_unstacked['M1'],
                     data_unstacked['M3'],
                     data_unstacked['M6'],
                     data_unstacked['M12'],
                     data_unstacked['M18'],
                     data_unstacked['M24'])
    
    Anova_dict[Category_name] = Anova[1]
    if Anova[1] < 0.1:
        print(Category_name)

        temp_df_stacked = temp_df_stacked.reset_index()
        res_tukey_hsd = pairwise_tukeyhsd(temp_df_stacked['Values'],temp_df_stacked['Months'],alpha=0.05)
        
        Tukey_df_results = pd.DataFrame(data=res_tukey_hsd._results_table.data[1:], columns=res_tukey_hsd._results_table.data[0])
        
        Tukey_df_results_significant = Tukey_df_results[Tukey_df_results['p-adj'] < 0.1]
    
        Tukey_df_results.to_excel(os.path.join(Folder_1,File_out_1),index=False)
        Tukey_df_results_significant.to_excel(os.path.join(Folder_1,File_out_2),index=False)
