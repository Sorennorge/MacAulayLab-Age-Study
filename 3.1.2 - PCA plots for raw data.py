# -*- coding: utf-8 -*-

import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
from matplotlib.patches import Ellipse
sns.set(style="white")
sns.set(style="ticks",font_scale=1.70)

### Function ###

def get_cov_ellipse(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.

    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by 
    vx, vy = eigvecs[:,0][0], eigvecs[:,0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals)
    return Ellipse(xy=centre, width=width, height=height,
                   angle=np.degrees(theta), **kwargs)


def fat_get_cov_ellipse(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.

    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by 
    vx, vy = eigvecs[:,0][0], eigvecs[:,0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals)
    return Ellipse(xy=centre, width=width, height=height*5,
                   angle=np.degrees(theta), **kwargs)
## Program parameter settings ##

x_pc = 1
y_pc = 2

## Folders ##

Folder1 = "Data/Normalized"
Folder2 = "Results/PCA"
Folder3 = "Results/PCA/data"

os.makedirs(Folder2, exist_ok=True)
os.makedirs(Folder3, exist_ok=True)

## Files ##

file_data = "DEseq2 rlog transformed.csv"

PCA_data_file = "PC_data_sorted.xlsx"

PCA_plot_file = "PCA_plot_age_PC{}_vs_PC{}.png".format(x_pc,y_pc)


## Load data ##
df_data = pd.read_csv(os.path.join(Folder1,file_data), sep=";",decimal=",")
df_data = df_data.rename(columns=({"Unnamed: 0":"Ensembl ID"})).set_index("Ensembl ID")
df_data_T = df_data.T

# data scaling
x_scaled = StandardScaler().fit_transform(df_data_T)

# set principal compnents #
n = 4
pca = PCA(n_components=n)

# transport data
pca_features = pca.fit_transform(x_scaled)

columns_n = []
for i in range(1,n+1,1):
    columns_n.append("PC{}".format(i))
# create dataframe with the n PC
pca_df = pd.DataFrame(
    data=pca_features, 
    columns=columns_n)

## Variance for plot 

Variance_PC_array = pca.explained_variance_ratio_

PC_1_V = round(Variance_PC_array[0]*100,2)
PC_2_V = round(Variance_PC_array[1]*100,2)
PC_3_V = round(Variance_PC_array[2]*100,2)
PC_4_V = round(Variance_PC_array[3]*100,2)

pc_list = [PC_1_V,PC_2_V,PC_3_V,PC_4_V]
# map target names to PCA features   


### Color scheme ###

M1 = "#07f2f2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#a422e6"
M24 = "#b3b3b3"

colorpallette = [M1,M3,M6,M12,M18,M24]

age_list = ['1 month old','3 months old','6 months old','12 months old','18 months old','24 months old']

pca_df['Age'] = ['1 month old','1 month old','1 month old',
                    '3 months old','3 months old','3 months old',
                    '6 months old','6 months old','6 months old',
                    '12 months old','12 months old','12 months old',
                    '18 months old','18 months old','18 months old',
                    '24 months old','24 months old','24 months old']

fig, ax = plt.subplots(figsize=(6, 6))

ax.set_facecolor('none')
plt.xticks([-100,-50,0,50,100,150])
plt.yticks([-150,-100,-50,0,50,100,150])
for i in range(0,6,1):
    sdata = pca_df[pca_df['Age']==age_list[i]]
    height_mean = np.mean(sdata['PC{}'.format(x_pc)])
    mass_mean = np.mean(sdata['PC{}'.format(y_pc)])
    cov = np.cov(sdata['PC{}'.format(x_pc)], sdata['PC{}'.format(y_pc)])
    if i == 1 or i == 3 or i == 6:
        e = fat_get_cov_ellipse(cov, (height_mean, mass_mean), 1.3,
                            fc=colorpallette[i], alpha=0.2)
    else:
        e = get_cov_ellipse(cov, (height_mean, mass_mean), 1.3,
                            fc=colorpallette[i], alpha=0.2)   
    ax.add_artist(e)

sns.scatterplot(ax = ax,
    x='PC{}'.format(x_pc), 
    y='PC{}'.format(y_pc), 
    data=pca_df, 
    hue='Age',
    palette=colorpallette,
    legend=False,linewidth=0.8, alpha = 0.7,edgecolor="black",s=75
    )
sns.despine(top=False, right=False, left=False, bottom=False)
#sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5))
plt.xlabel( "PC {} ({} %)".format(x_pc,pc_list[x_pc-1]))
plt.ylabel( "PC {} ({} %)".format(y_pc,pc_list[y_pc-1]))
#plt.xticks([-15,-5,5,15,25,35])
file_out = "PCA_PC{}_vs_PC{}_ellipse.png".format(x_pc,y_pc)

#plt.savefig(os.path.join(Folder2,PCA_plot_file),dpi=600,bbox_inches='tight')

plt.xticks([-100,-50,0,50,100,150])
plt.yticks([-150,-100,-50,0,50,100,150])
#plt.savefig(os.path.join(Folder2,file_out),dpi=600,bbox_inches='tight')
plt.show()


# If you want sorted the pca data #
#pca_df.to_excel(os.path.join(Folder3,PCA_data_file),index=False)

#df_Variance = pd.DataFrame({ 'Variance': Variance_PC_array })
#df_Variance.to_excel(os.path.join(Folder3,"Variance.xlsx"),index=False)