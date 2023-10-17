# -*- coding: utf-8 -*-

### Metabolic processes ###

## Libraries ##

import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba
from matplotlib.cm import register_cmap
from matplotlib.collections import LineCollection
from matplotlib.path import Path
from matplotlib.patches import PathPatch

## input ##

#Save_fig = "yes"
Save_fig = "no"

# Resolution #
n_bins = 50

# figure dimensions #

fig_x = 12
fig_y = 5

## limits ##

y_lim = [-0.5,1]

# outer line #
outer_line_width = 1

# outer points #
outer_point_size = 2

# Circles #
circle_edge_width = 2

# Line between points #
line_thickness = 6
line_edge_thickness = line_thickness+2

# Alpha settings #
Alpha_line = 1.0
alpha_fill = 0.1 # Needs changing depending on the resolution 

## Font sizes ##
x_axis_font_size = 20
y_axis_font_size = 20
y_label_fontsize = 26
x_label_fontsize = 26
## Functions ##

title_y_coord = 1
title_x_coord = 0.5
title_fontsize = 40

def split(start, end, segments):
    x_delta = (end[0] - start[0]) / float(segments)
    y_delta = (end[1] - start[1]) / float(segments)
    points = []
    for i in range(1, segments):
        points.append([start[0] + i * x_delta, start[1] + i * y_delta])
    return [start] + points + [end]

def linspace(start, stop, step=1.):
  return np.linspace(start, stop, int((stop - start) / step + 1))

## Folders ##

Folder1 = "Data/Panther/Biological functions"
Folder2 = "Data/Normalized"

Folder3 = "Results/Metabolic processes"
os.makedirs(Folder3,exist_ok=True)

Folder4 = "Data/Pvalue data/Metabolic processes"
os.makedirs(Folder4,exist_ok=True)

Folder_cmap = "Data/cmaps"
os.makedirs(Folder_cmap,exist_ok=True)

## Files ##

File1 = "Biological functions overview.xlsx"
File2 = "BF_Metabolic_process_all_genes.txt"
File3 = "DEseq2 Normalized data.csv"
File4 = "Metabolic_data_stack.xlsx"

File_out_1 = "Metabolic processes zscore.png"

## Create color scheme ###

# relative positions of colors in cmap/palette 

# colors used #
M1 = "#07F2F2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#A422E6"
M24 = "#B3B3B3"

pos = [0.0,1.0]

## Different segments ##

# Convert to RGB colors #

M1_RGB = to_rgba(M1)
M3_RGB = to_rgba(M3)
M6_RGB = to_rgba(M6)
M12_RGB = to_rgba(M12)
M18_RGB = to_rgba(M18)
M24_RGB = to_rgba(M24)

RGB_list = [M1_RGB,M3_RGB,M6_RGB,M12_RGB,M18_RGB,M24_RGB]

# color lists #

colors_M1_M3 = [M1_RGB,M3_RGB]
colors_M3_M6 = [M3_RGB,M6_RGB]
colors_M6_M12 = [M6_RGB,M12_RGB]
colors_M12_M18 = [M12_RGB,M18_RGB]
colors_M18_M24 = [M18_RGB,M24_RGB]

# Create cmap for the different segments #

cmap_M1_M3 = LinearSegmentedColormap.from_list("cmap_name", colors_M1_M3, N=n_bins)
cmap_M3_M6 = LinearSegmentedColormap.from_list("cmap_name", colors_M3_M6, N=n_bins)
cmap_M6_M12 = LinearSegmentedColormap.from_list("cmap_name", colors_M6_M12, N=n_bins)
cmap_M12_M18 = LinearSegmentedColormap.from_list("cmap_name", colors_M12_M18, N=n_bins)
cmap_M18_M24 = LinearSegmentedColormap.from_list("cmap_name", colors_M18_M24, N=n_bins)

# Create cmap color list #
cmap_stack = [cmap_M1_M3,cmap_M3_M6,cmap_M6_M12,cmap_M12_M18,cmap_M18_M24]

## Global variables ##

Sample_list = ['M1','M3','M6','M12','M18','M24']

## Load data ##

## Metabolic processes - LRT ##

df_metabo_LRT = pd.read_excel(os.path.join(Folder1,File1),usecols=[0,1,2],names=['Ensembl ID','Gene','BF'])

BF_dict = {}
df_metabo_LRT_list = []
for index,rows in df_metabo_LRT.iterrows():
    key = rows['Ensembl ID']
    item = rows['BF']
    BF_temp = item.split(";")
    if "metabolic process" in BF_temp:
        BF_dict[key] =  BF_temp
        df_metabo_LRT_list.append(key)

## Normalized data ##

df_data = pd.read_csv(os.path.join(Folder2,File3),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))

### subset data ###

df_data_LRT = df_data[df_data['Ensembl ID'].isin(df_metabo_LRT_list)].copy()


for key in Sample_list:
    df_data_LRT[key] = df_data_LRT[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)

## Reduce dataframe subsets ##
df_data_LRT_mean = df_data_LRT[['Ensembl ID']+Sample_list].copy()
df_data_LRT_mean = df_data_LRT_mean.set_index('Ensembl ID')
df_data_LRT_mean = df_data_LRT_mean.T
df_data_LRT_zscored_T = df_data_LRT_mean.apply(stats.zscore)
df_data_LRT_zscored = df_data_LRT_zscored_T.T

df_zscored_stacked = df_data_LRT_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})

## Save data for statistical analysis ##

df_zscored_stacked.to_excel(os.path.join(Folder4,File4),index=False)


df_zscored_std = df_zscored_stacked.groupby('Months')['Values'].std()

df_zscored_stacked_1_3 = df_zscored_stacked[(df_zscored_stacked['Months'] == 'M1') | (df_zscored_stacked['Months'] == 'M3')]
df_zscored_stacked_3_6 = df_zscored_stacked[(df_zscored_stacked['Months'] == 'M3') | (df_zscored_stacked['Months'] == 'M6')]
df_zscored_stacked_6_12 = df_zscored_stacked[(df_zscored_stacked['Months'] == 'M6') | (df_zscored_stacked['Months'] == 'M12')]
df_zscored_stacked_12_18 = df_zscored_stacked[(df_zscored_stacked['Months'] == 'M12') | (df_zscored_stacked['Months'] == 'M18')]
df_zscored_stacked_18_24 = df_zscored_stacked[(df_zscored_stacked['Months'] == 'M18') | (df_zscored_stacked['Months'] == 'M24')]

data_stack = [df_zscored_stacked_1_3,df_zscored_stacked_3_6,df_zscored_stacked_6_12,df_zscored_stacked_12_18,df_zscored_stacked_18_24]

## confidence interval - for fill ##

confidence_interval_M1 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M1']['Values']))
confidence_interval_M3 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M3']['Values']))
confidence_interval_M6 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M6']['Values']))
confidence_interval_M12 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M12']['Values']))
confidence_interval_M18 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M18']['Values']))
confidence_interval_M24 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked[df_zscored_stacked['Months'] == 'M24']['Values']))

# Create lists of confidence intervals #
confidence_list = [confidence_interval_M1,
                   confidence_interval_M3,
                   confidence_interval_M6,
                   confidence_interval_M12,
                   confidence_interval_M18,
                   confidence_interval_M24]


## Initiation of figure
scatter_array = []
line_array = []
fig, ax = plt.subplots(figsize=(fig_x, fig_y))
for i in range(0,5,1):
    ax = sns.lineplot(data=data_stack[i], x="Months",color='black',
                       y="Values",alpha=0.75,
                       marker="o",
                       dashes=False,legend=None, markersize=10,errorbar=None)
    x, y = ax.get_lines()[0].get_data()
    line_array.append([[x[0],x[1]],[y[0],y[1]]])
    ax.get_lines()[0].remove()
    if i == 0:
        scatter_array.append([x[0], y[0]])
        scatter_array.append([x[1], y[1]])
    else:
        scatter_array.append([x[1], y[1]])
plt.ylim(y_lim)

### create black lines ###
for i in range(0,5,1):
    x1 = line_array[i][0][0]
    x2 = line_array[i][0][1]
    y1 = line_array[i][1][0]
    y2 = line_array[i][1][1]
    ax.plot([x1,x2],[y1,y2],lw=line_edge_thickness,color='Black')
    ax.plot([x1,x2],[y1,y2],lw=line_thickness,color='White')

for i in range(0,6,1):
    if i == 0:
        # Plot upper line #
        ax.plot([i,i+1],[confidence_list[i][1],confidence_list[i+1][1]],lw=outer_line_width,color='Black')
        # Plot lower line #
        ax.plot([i,i+1],[confidence_list[i][0],confidence_list[i+1][0]],lw=outer_line_width,color='Black')
        #Plot begining line #
        ax.plot([i,i],[confidence_list[i][0],confidence_list[i][1]],lw=outer_line_width,color='Black')    
    elif i == 5:
        #Plot end line #
        ax.plot([i,i],[confidence_list[i][0],confidence_list[i][1]],lw=outer_line_width,color='Black')
    else:
        # Plot upper line #
        ax.plot([i,i+1],[confidence_list[i][1],confidence_list[i+1][1]],lw=outer_line_width,color='Black')
        # Plot lower line #
        ax.plot([i,i+1],[confidence_list[i][0],confidence_list[i+1][0]],lw=outer_line_width,color='Black')

for i in range(0,5,1):
    # Segment lines into 50 bits #
    x1 = line_array[i][0][0]
    x2 = line_array[i][0][1]
    y1 = line_array[i][1][0]
    y2 = line_array[i][1][1]
    x_seg = np.linspace(x1,x2, n_bins)
    y_seg = np.linspace(y1,y2, n_bins)
    points = np.array([x_seg, y_seg]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    cols = np.linspace(0,1,len(x_seg))
    lc = LineCollection(segments, cmap=cmap_stack[i],linestyle="solid")
    lc.set_array(cols)
    lc.set_linewidth(line_thickness)
    line = ax.add_collection(lc)
  
for i in range(0,len(scatter_array),1):
    point_x = scatter_array[i][0]
    point_y = scatter_array[i][1]
    ax.plot(point_x, point_y, marker='o',color=RGB_list[i],markeredgecolor='black',markeredgewidth=circle_edge_width, ms=25)

for i in range(0,5,1):
    if i < 5:
        x1 = line_array[i][0][0]
        x2 = line_array[i][0][1]
        y1 = line_array[i][1][0]
        y2 = line_array[i][1][1]
        x_seg = np.linspace(x1,x2, n_bins)
        y_seg = np.linspace(y1,y2, n_bins)
        ## Set backgroup color ##
        # set colormap #
        cmap = cmap_stack[i]
        cols = np.linspace(0,1,n_bins)
        y_upper = np.linspace(confidence_list[i][1],confidence_list[i+1][1], n_bins)
        y_lower = np.linspace(confidence_list[i][0],confidence_list[i+1][0], n_bins)
        for ii in range(0,len(x_seg)-1,1):
            ax.fill_between([x_seg[ii], 
                             x_seg[ii+1]],
                            [y_lower[ii],
                             y_lower[ii+1]],
                            [y_upper[ii],
                             y_upper[ii+1]],
                            color=cmap(cols[ii]),
                            alpha=alpha_fill)


## Format plot axes ##
plt.title("Metabolic processes",ha='center',va='center',y=title_y_coord,x=title_x_coord,fontsize=title_fontsize)
ax.tick_params(bottom=False)
plt.ylim(y_lim)
plt.yticks([y_lim[0],0,0.5,y_lim[1]],fontsize=y_axis_font_size)
sns.despine(fig=fig, ax=ax, top=True, right=True, left=False, bottom=True)
plt.xticks(Sample_list,['1M','3M','6M','12M','18M','24M'],fontsize=x_axis_font_size)
plt.ylabel("Z-score", fontsize=y_label_fontsize)
plt.xlabel("", fontsize=x_label_fontsize)
#adjust x-axis label position 
ax.xaxis.set_label_coords(0.5, -.1)

if Save_fig == "yes":
    plt.savefig(os.path.join(Folder3,File_out_1),dpi=600,bbox_inches='tight')
else:
    plt.show()