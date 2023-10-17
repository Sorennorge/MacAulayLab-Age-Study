# -*- coding: utf-8 -*-

### PCA from sorted data ###

import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.text as mtext
sns.set(style="white")
sns.set(style="ticks",font_scale=2)

### Function ###

class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, orig_handle, usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title

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

x_pc = 2
y_pc = 3

## Folders ##

Folder1 = "Data/Normalized"
Folder2 = "Results/PCA"
Folder3 = "Results/PCA/data"

os.makedirs(Folder2, exist_ok=True)
os.makedirs(Folder3, exist_ok=True)

## Files ##

#file_data = "DEseq2 Normalized data.csv"
#file_data = "DEseq2 rlog transformed.csv"
#file_data = "DEseq2 rlog 2 transformed.csv"
#file_targets = "PCA_targets.xlsx"

PCA_data_file = "PC_data_sorted_V1.xlsx"

PCA_plot_file = "Stored_PCA_plot_age_PC{}_vs_PC{}.png".format(x_pc,y_pc)

file_out = "Stored_PCA_PC{}_vs_PC{}_ellipse.png".format(x_pc,y_pc)

## Load data ##
pca_df = pd.read_excel(os.path.join(Folder3,PCA_data_file))
df_variance = pd.read_excel(os.path.join(Folder3,"Variance_V1.xlsx"))

## Global variables ##


age_list = ['1 month old','3 months old','6 months old','12 months old','18 months old','24 months old']


## Variance for plot 

Variance_PC_array = df_variance['Variance']

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

fig, ax = plt.subplots(figsize=(12, 10))

ax.set_facecolor('none')
plt.xticks([-100,-50,0,50,100,150])
plt.yticks([-150,-100,-50,0,50,100,150])
for i in range(0,6,1):
    sdata = pca_df[pca_df['Age']==age_list[i]]
    height_mean = np.mean(sdata['PC{}'.format(x_pc)])
    mass_mean = np.mean(sdata['PC{}'.format(y_pc)])
    cov = np.cov(sdata['PC{}'.format(x_pc)], sdata['PC{}'.format(y_pc)])
    #print(cov)
    if i == 1 or i == 3 or i == 6:
        #print(age_list[i])
        e = fat_get_cov_ellipse(cov, (height_mean, mass_mean), 1.3,
                            fc=colorpallette[i], alpha=0.3)
    else:
        e = get_cov_ellipse(cov, (height_mean, mass_mean), 1.3,
                            fc=colorpallette[i], alpha=0.3)   
    ax.add_artist(e)

g = sns.scatterplot(ax = ax,
    x='PC{}'.format(x_pc), 
    y='PC{}'.format(y_pc), 
    data=pca_df, 
    hue='Age',
    palette=colorpallette,
    legend=None,linewidth=0.9, alpha = 0.90,edgecolor="black",s=80)
sns.despine(top=True, right=True, left=False, bottom=False)

# coordinates for text #
M1_coord = [100,-5]
M3_coord = [0,-60]
M6_coord = [22,-60]
M12_coord = [43,-60]
M18_coord = [-30,-115]
M24_coord = [-3,-115]

# add text and dots #

# M1 #
plt.text(M1_coord[0],M1_coord[1],"1M",fontsize=25)
plt.scatter(M1_coord[0]-2,M1_coord[1]+4.5, s=80, facecolors=colorpallette[0],alpha = 0.90,linewidth=0.9, edgecolors='Black')
# M3 #
plt.text(M3_coord[0],M3_coord[1],"3M",fontsize=25)
plt.scatter(M3_coord[0]-3,M3_coord[1]+4.5, s=80, facecolors=colorpallette[1],alpha = 0.90,linewidth=0.9, edgecolors='Black')
# M6 #
plt.text(M6_coord[0],M6_coord[1],"6M",fontsize=25)
plt.scatter(M6_coord[0]-3,M6_coord[1]+4.5, s=80, facecolors=colorpallette[2],alpha = 0.90,linewidth=0.9, edgecolors='Black')
# M12 #
plt.text(M12_coord[0],M12_coord[1],"12M",fontsize=25)
plt.scatter(M12_coord[0]-2,M12_coord[1]+4.5, s=80, facecolors=colorpallette[3],alpha = 0.90,linewidth=0.9, edgecolors='Black')
# M18 #
plt.text(M18_coord[0],M18_coord[1],"18M",fontsize=25)
plt.scatter(M18_coord[0]-2,M18_coord[1]+4.5, s=80, facecolors=colorpallette[4],alpha = 0.90,linewidth=0.9, edgecolors='Black')
# M24 #
plt.text(M24_coord[0],M24_coord[1],"24M",fontsize=25)
plt.scatter(M24_coord[0]-3,M24_coord[1]+4.5, s=80, facecolors=colorpallette[5],alpha = 0.90,linewidth=0.9, edgecolors='Black')


## Set axis label and ticks ##

plt.xlabel( "PC {}".format(x_pc),labelpad=10)
plt.ylabel( "PC {}".format(y_pc),labelpad=-10)


plt.xticks([-100,-50,0,50,100,150])
plt.yticks([-150,-100,-50,0,50,100,150])

## Add three group elipses ##

ellipse_1 = Ellipse(xy=(52.5, 30.5), width=168, height=50, 
                        edgecolor='Grey', fill=False,ls='--', lw=1,angle=5)
ax.add_patch(ellipse_1)

ellipse_2 = Ellipse(xy=(21.5, -18), width=94, height=52, 
                        edgecolor='Grey', fill=False,ls='--', lw=1,angle=-2)
ax.add_patch(ellipse_2)

ellipse_3 = Ellipse(xy=(-60, -3), width=260, height=63, 
                        edgecolor='Grey', fill=False,ls='--', lw=1,angle=90)
ax.add_patch(ellipse_3)

plt.savefig(os.path.join(Folder2,file_out),dpi=600,bbox_inches='tight')
#plt.show()
