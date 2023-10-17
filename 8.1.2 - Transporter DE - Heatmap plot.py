# -*- coding: utf-8 -*-
"""
Created on Sun May 21 12:36:30 2023

@author: dcs839
"""

### Transport panther long list ###

import os
import pandas as pd
import numpy as np
from natsort import index_natsorted
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

## input ##

n_bins = 100

#fmt_value = '.3g'
#fmt_value = '.2g'
fmt_value = '.0%'
fmt_value_TT = ".0f"

## Folders #

Folder1 = "Data/Panther/NEW/Protein classes"
Folder2 = "Results/Transport/Heatmap"
Folder3 = "Data/Normalized"

## Files ##

#File_1 = "Transport_panther_dataframe.xlsx"
File_1 = "Transporter overview.xlsx"
#File_2 = "Transport_mystudy_dataframe.xlsx"
File_3 = "Normalized_TMM.xlsx"

File_out_1 = "Heatmap_Panther_transport_percentage.png"

## Creating colors ##

M1 = "#07F2F2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#A422E6"
M24 = "#B3B3B3"
white = "#FFFFFF"

M1_RGB = to_rgba(M1)
M3_RGB = to_rgba(M3)
M6_RGB = to_rgba(M6)
M12_RGB = to_rgba(M12)
M18_RGB = to_rgba(M18)
M24_RGB = to_rgba(M24)
white_RGB = to_rgba(white)

colors_M1 = [white_RGB,M1_RGB]
colors_M3 = [white_RGB,M3_RGB]
colors_M6 = [white_RGB,M6_RGB]
colors_M12 = [white_RGB,M12_RGB]
colors_M18 = [white_RGB,M18_RGB]
colors_M24 = [white_RGB,M24_RGB]


cmap_M1 = LinearSegmentedColormap.from_list("cmap_name", colors_M1, N=n_bins)
cmap_M3 = LinearSegmentedColormap.from_list("cmap_name", colors_M3, N=n_bins)
cmap_M6 = LinearSegmentedColormap.from_list("cmap_name", colors_M6, N=n_bins)
cmap_M12 = LinearSegmentedColormap.from_list("cmap_name", colors_M12, N=n_bins)
cmap_M18 = LinearSegmentedColormap.from_list("cmap_name", colors_M18, N=n_bins)
cmap_M24 = LinearSegmentedColormap.from_list("cmap_name", colors_M24, N=n_bins)

## Load data ##

df_panther = pd.read_excel(os.path.join(Folder1,File_1))

# rename AABR07029596.1 to Slc22a4
#df_panther.loc[df_panther["Gene_name"] == "AABR07029596.1", "Gene_name"] = "Slc22a4"

df_expression = pd.read_excel(os.path.join(Folder3,File_3))

df_express_subset = df_expression[df_expression['Ensembl ID'].isin(df_panther['Ensembl ID'])].copy().set_index('Ensembl ID')

df_express_subset['TMM'] = df_express_subset.mean(axis=1)

df_express_subset = df_express_subset[['TMM']]

df_panther_mean = pd.concat([df_panther.set_index('Ensembl ID'),df_express_subset],join='inner',axis=1)

df_panther_sorted = df_panther_mean.sort_values(by=['TMM'],ascending=False,ignore_index=True)

df_panther_sorted.TMM = df_panther_sorted.TMM.round(0)

df_panther_sorted_subset = df_panther_sorted[df_panther_sorted['TMM'] >= 100]

df_panther_M1 = df_panther_sorted_subset.set_index('Gene')[['M1']]
df_panther_M3 = df_panther_sorted_subset.set_index('Gene')[['M3']]
df_panther_M6 = df_panther_sorted_subset.set_index('Gene')[['M6']]
df_panther_M12 = df_panther_sorted_subset.set_index('Gene')[['M12']]
df_panther_M18 = df_panther_sorted_subset.set_index('Gene')[['M18']]
df_panther_M24 = df_panther_sorted_subset.set_index('Gene')[['M24']]

df_panther_TMM = df_panther_sorted_subset.set_index('Gene')[['TMM']]

sns.set(font_scale=3)
fig = plt.figure(figsize=(20,50))

gs = fig.add_gridspec(8,8,width_ratios=[1,1,1,1,1,1,0.1,1],height_ratios=[1,1,1,1,1,1,0,0.1],hspace=0.1,wspace=0)
gs1 = gs[-8:-1].subgridspec(1, 6, hspace=0, wspace=0.3)
ax1 = fig.add_subplot(gs[:-1, 0])
ax2 = fig.add_subplot(gs[:-1, 1])
ax3 = fig.add_subplot(gs[:-1, 2])
ax4 = fig.add_subplot(gs[:-1, 3])
ax5 = fig.add_subplot(gs[:-1, 4])
ax6 = fig.add_subplot(gs[:-1, 5])
ax7 = fig.add_subplot(gs[:-1, 6])
ax7.set_visible(False)
ax8 = fig.add_subplot(gs[:-1, 7])

# color bars 
ax_C1 = fig.add_subplot(gs1[0, 0])
ax_C2 = fig.add_subplot(gs1[0, 1])
ax_C3 = fig.add_subplot(gs1[0, 2])
ax_C4 = fig.add_subplot(gs1[0, 3])
ax_C5 = fig.add_subplot(gs1[0, 4])
ax_C6 = fig.add_subplot(gs1[0, 5])

# Month 1 #
sns.heatmap(df_panther_M1,cbar=False,
                 annot_kws={'color':'black'},
                 cbar_kws={"shrink": .25},
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M1,
                 vmin=-2.5, vmax=5,ax=ax1)

# Month 3 #
sns.heatmap(df_panther_M3,cbar=False,
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M3,
                 vmin=-2.5, vmax=5,ax=ax2)

# Month 6 #
sns.heatmap(df_panther_M6,cbar=False,
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M6,
                 vmin=-2.5, vmax=5,ax=ax3)

# Month 12 #
sns.heatmap(df_panther_M12,cbar=False,
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M12,
                 vmin=-2.5, vmax=5,ax=ax4)

# Month 18 #
sns.heatmap(df_panther_M18,cbar=False,
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M18,
                 vmin=-2.5, vmax=5,ax=ax5)

# Month 24 #
sns.heatmap(df_panther_M24,cbar=False,
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value,
                 annot=True,cmap=cmap_M24,
                 vmin=-2.5, vmax=5,ax=ax6)

# TMM #

sns.heatmap(df_panther_TMM,cbar=False,
                 annot_kws={'color':'black'},
                 cbar_kws={"shrink": .25},
                 linecolor='black',linewidths=1.5,
                 clip_on=False,
                 fmt=fmt_value_TT,
                 annot=True,cmap=ListedColormap(['white']),ax=ax8)

# ylabels #
ax1.set_ylabel('')
ax2.set_ylabel('')
ax3.set_ylabel('')
ax4.set_ylabel('')
ax5.set_ylabel('')
ax6.set_ylabel('')
ax8.set_ylabel('')
# xlabels #
ax1.set_xticks([])
ax2.set_xticks([])
ax3.set_xticks([])
ax4.set_xticks([])
ax5.set_xticks([])
ax6.set_xticks([])
ax8.set_xticks([])
# remove genes -> yticks #
ax2.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
ax5.set_yticks([])
ax6.set_yticks([])
ax8.set_yticks([])
# set title #
title_pad = 10
ax1.set_title('1M',pad=title_pad)
ax2.set_title('3M',pad=title_pad)
ax3.set_title('6M',pad=title_pad)
ax4.set_title('12M',pad=title_pad)
ax5.set_title('18M',pad=title_pad)
ax6.set_title('24M',pad=title_pad)
ax8.set_title('TMM',pad=title_pad)



# Color bars #
cbar1 = fig.colorbar(ax1.get_children()[0], ax=ax_C1, cax=ax_C1,orientation="horizontal",ticks=[0, 5])
cbar1.outline.set_color('Black')
cbar1.outline.set_linewidth(2)
cbar1.ax.set_xticklabels(['0%', "500%"])

cbar2 = fig.colorbar(ax2.get_children()[0], ax=ax_C2, cax=ax_C2,orientation="horizontal",ticks=[0, 5])
cbar2.outline.set_color('Black')
cbar2.outline.set_linewidth(2)
cbar2.ax.set_xticklabels(['0%', "500%"])

cbar3 = fig.colorbar(ax3.get_children()[0], ax=ax_C3, cax=ax_C3,orientation="horizontal",ticks=[0, 5])
cbar3.outline.set_color('Black')
cbar3.outline.set_linewidth(2)
cbar3.ax.set_xticklabels(['0%', "500%"])

cbar4 = fig.colorbar(ax4.get_children()[0], ax=ax_C4, cax=ax_C4,orientation="horizontal",ticks=[0, 5])
cbar4.outline.set_color('Black')
cbar4.outline.set_linewidth(2)
cbar4.ax.set_xticklabels(['0%', "500%"])

cbar5 = fig.colorbar(ax5.get_children()[0], ax=ax_C5, cax=ax_C5,orientation="horizontal",ticks=[0, 5])
cbar5.outline.set_color('Black')
cbar5.outline.set_linewidth(2)
cbar5.ax.set_xticklabels(['0%', "500%"])

cbar6 = fig.colorbar(ax6.get_children()[0], ax=ax_C6, cax=ax_C6,orientation="horizontal",ticks=[0, 5])
cbar6.outline.set_color('Black')
cbar6.outline.set_linewidth(2)
cbar6.ax.set_xticklabels(['0%', "500%"])

#plt.show()
plt.savefig(os.path.join(Folder2,File_out_1),dpi=600,bbox_inches='tight')