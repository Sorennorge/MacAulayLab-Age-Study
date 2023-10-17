# -*- coding: utf-8 -*-

### P_values for MF ###

import pandas as pd
import os
from scipy import stats
import numpy as np
from scipy.stats import f_oneway
from scipy.stats import tukey_hsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd

## Folders ##

Folder_1 = "Data/Panther/LRT diff/Enrichment data"
Folder_2 = "Data/Pvalue data/MF"


## Files ##
File_1 = "MF_enrichment_data_all.xlsx"

File_out_1 = "Tukey_results_MF.xlsx"
File_out_2 = "Tukey_results_MF_significant.xlsx"

# Categories #

df_enrichment = pd.read_excel(os.path.join(Folder_1,File_1))

df_enrichment_subset = df_enrichment[df_enrichment['Percentage'] >= 1.0]

MF_list = df_enrichment_subset['Class'].tolist()
MF_list = MF_list[:-1]

p_value_dict = {}
p_value_dict_second = {}
Turkey_dict = {}

for i_value in MF_list:
    data_file = "{} pvalue data.xlsx".format(i_value)
    data = pd.read_excel(os.path.join(Folder_2,data_file))
    data_stacked = data.set_index('MF').stack().to_frame().reset_index()
    data_stacked = data_stacked.rename(columns={'level_1':'Months',0:"Values"})
    
    Anova = f_oneway(data['M1'],data['M3'],data['M6'],data['M12'],data['M18'],data['M24'])
    res_tukey_hsd = pairwise_tukeyhsd(data_stacked['Values'],data_stacked['Months'],alpha=0.05)
    
    p_value_dict[i_value] = Anova[1]
    Turkey_dict[i_value] = res_tukey_hsd
    # Create and store results #
    df_temp = pd.DataFrame(data=res_tukey_hsd._results_table.data[1:], columns=res_tukey_hsd._results_table.data[0])
    df_temp['MF'] = i_value
    if i_value == MF_list[0]:
        Tukey_df_results = df_temp
    else:
        Tukey_df_results = pd.concat([Tukey_df_results,df_temp])

Tukey_df_results_significant = Tukey_df_results[Tukey_df_results['p-adj'] < 0.05]

Tukey_df_results.to_excel(os.path.join(Folder_2,File_out_1),index=False)
Tukey_df_results_significant.to_excel(os.path.join(Folder_2,File_out_2),index=False)
