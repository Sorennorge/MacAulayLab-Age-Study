# -*- coding: utf-8 -*-

### P_values for Metabolic processes ###

import pandas as pd
import os
from scipy import stats
import numpy as np
from scipy.stats import f_oneway
from scipy.stats import tukey_hsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd


## Folders ##

Folder_1 = "Data/Pvalue data/Metabolic processes"

# Files ##

File_1 = "TCA data stacked.xlsx"

File_out_1 = "Tukey_results_TCA.xlsx"
File_out_2 = "Tukey_results_TCA_significant.xlsx"

# Categories #


#i_value = MF_list[0]

data_stacked = pd.read_excel(os.path.join(Folder_1,File_1)).set_index(["Ensembl ID",'Months'])

data_unstacked = data_stacked.unstack(level=-1)
data_unstacked = data_unstacked.droplevel(None, axis=1)

print(data_unstacked.columns)

Anova = f_oneway(data_unstacked['M1'],
                 data_unstacked['M3'],
                 data_unstacked['M6'],
                 data_unstacked['M12'],
                 data_unstacked['M18'],
                 data_unstacked['M24'])

data_stacked = data_stacked.reset_index()
res_tukey_hsd = pairwise_tukeyhsd(data_stacked['Values'],data_stacked['Months'],alpha=0.05)

Tukey_df_results = pd.DataFrame(data=res_tukey_hsd._results_table.data[1:], columns=res_tukey_hsd._results_table.data[0])

Tukey_df_results_significant = Tukey_df_results[Tukey_df_results['p-adj'] < 0.05]

Tukey_df_results.to_excel(os.path.join(Folder_1,File_out_1),index=False)
Tukey_df_results_significant.to_excel(os.path.join(Folder_1,File_out_2),index=False)


