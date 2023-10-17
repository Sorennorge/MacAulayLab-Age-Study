# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:53:40 2023

@author: dcs839
"""

### Transport usual suspects ###

import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection



## input ##


i_value = 1

#z_aksis_print = "yes"
z_aksis_print = "no"

#Save_fig = "yes"
Save_fig = "no"

#last_plot = "yes"
last_plot = "no"

# Resolution #
n_bins = 50

# figure dimensions #

fig_x = 12
fig_y = 2

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

y_limit_upper = 3
y_limit_lower = -3




## Files ##

File1 = "Transport of interest.xlsx"
File2 = "DEseq2 Normalized data.csv"

## Folders ##

Folder1 = "Data/Gene info"
Folder2 = "Data/normalized"

Folder_out = "Results/Transport/Interest - time series"
os.makedirs(Folder_out,exist_ok=True)

## Load data ##

# Transport of interest #
df_Transport_of_interest = pd.read_excel(os.path.join(Folder1,File1))

#ensembl_list = ["ENSRNOG00000030019","ENSRNOG00000015971","ENSRNOG00000014504","ENSRNOG00000013963","ENSRNOG00000031312"]
#gene_alias = ["ATP1A1 (NKA)","SLC12A2 (NKCC1)","IL1R1","IL-6ST","TNFR1"]
#ensembl_list = ["ENSRNOG00000016050","ENSRNOG00000016374","ENSRNOG00000017392"]
#gene_alias = ['FGFR1','FGFR2',"FGF2"]

#df_Transport_of_interest = pd.DataFrame({'Rat Ensembl ID':ensembl_list,'Gene_Alias':gene_alias})

## Normalized data ##

df_data = pd.read_csv(os.path.join(Folder2,File2),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))

### Color mapping ##

M1 = "#07F2F2"
M3 = "#00F000"
M6 = "#FFA500"
M12 = "#FF0000"
M18 = "#A422E6"
M24 = "#B3B3B3"

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


## Modulate data to zscored values ##

df_Transport_of_interest_data2 = pd.concat([df_Transport_of_interest.set_index('Rat Ensembl ID'),df_data.set_index('Ensembl ID')],join="inner",axis=1)
df_Transport_of_interest_data2 = df_Transport_of_interest_data2.reset_index().set_index(['index','Gene_Alias'])
df_Transport_of_interest_zscored = df_Transport_of_interest_data2.apply(stats.zscore,axis=1)
df_Transport_of_interest_zscored = df_Transport_of_interest_zscored.assign(m=df_Transport_of_interest_zscored[['Sample.M1.R1','Sample.M1.R2','Sample.M1.R3']].mean(axis=1)).sort_values('m').drop('m', axis=1)


Sample_list = ['M1','M3','M6','M12','M18','M24']

for key in Sample_list:
    df_Transport_of_interest_zscored = df_Transport_of_interest_zscored.rename(columns=({'Sample.{}.R1'.format(key):'{}'.format(key),
                                                      'Sample.{}.R2'.format(key):'{}'.format(key),
                                                      'Sample.{}.R3'.format(key):'{}'.format(key)}))
## Drop ensembl ID ##
df_Transport_of_interest_zscored = df_Transport_of_interest_zscored.droplevel('index')
#df_Transport_of_interest_zscored = df_Transport_of_interest_zscored.droplevel('Gene_Alias')

transport_of_interest_list_run = df_Transport_of_interest_zscored.index.tolist()
#transport_of_interest_list_run = ["ENSRNOG00000030019","ENSRNOG00000015971"]


## create data for plot ##

#df_Transport_of_interest_zscored = df_Transport_of_interest_zscored[df_Transport_of_interest_zscored.index[0] == transport_of_interest_list_run[i_value]]
df_Transport_of_interest_zscored = df_Transport_of_interest_zscored[df_Transport_of_interest_zscored.index == transport_of_interest_list_run[i_value]]

df_zscored_stacked = df_Transport_of_interest_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})

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
fig1, ax = plt.subplots(figsize=(fig_x, fig_y))
for i in range(0,5,1):
    ax = sns.lineplot(data=data_stack[i], x="Months",color='black',
                       y="Values",alpha=0.75,
                       marker="o",
                       dashes=False,legend=None, markersize=10,ci=None)
    x, y = ax.get_lines()[0].get_data()
    line_array.append([[x[0],x[1]],[y[0],y[1]]])
    ax.get_lines()[0].remove()
    if i == 0:
        scatter_array.append([x[0], y[0]])
        scatter_array.append([x[1], y[1]])
    else:
        scatter_array.append([x[1], y[1]])
plt.ylim([y_limit_lower,y_limit_upper])


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
    #print(i)
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
    #ax.plot([x1,x2],[y1,y2],lw=6,color='Black')
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

plt.title(transport_of_interest_list_run[i_value], fontsize=24)
if last_plot == "yes":
    ## For last plot ##
    ## Format plot axes ##
    ax.tick_params(bottom=False)
    plt.ylim([y_limit_lower,y_limit_upper])
    plt.yticks([-2,0,2],fontsize=y_axis_font_size)
    sns.despine(fig=fig1, ax=ax, top=True, right=True, left=False, bottom=True)
    plt.xticks(Sample_list,['1','3','6','12','18','24'],fontsize=x_axis_font_size)
    plt.ylabel("Z-score", fontsize=y_label_fontsize)
    plt.xlabel("Months", fontsize=x_label_fontsize)
    #adjust x-axis label position 
    ax.xaxis.set_label_coords(0.5, -.2)
else:
    if z_aksis_print == "yes":
        ax.tick_params(bottom=False)
        ax.tick_params(left=True)
        plt.ylim([y_limit_lower,y_limit_upper])
        ax.spines['left'].set_bounds(0, 1)
        plt.yticks([0,1],fontsize=y_axis_font_size)
        #plt.yticks([-2,0,2],fontsize=y_axis_font_size)
        sns.despine(fig=fig1, ax=ax, top=True, right=True, left=False, bottom=True)
        plt.xticks("")
        plt.ylabel("Δ z-score", fontsize=y_label_fontsize,y=0.76)
        plt.xlabel("")
        #adjust x-axis label position 
        ax.xaxis.set_label_coords(0.5, -.2)
        
    else:
        ax.tick_params(bottom=False)
        ax.tick_params(left=True)
        plt.ylim([y_limit_lower,y_limit_upper])
        plt.yticks([],fontsize=y_axis_font_size)
        sns.despine(fig=fig1, ax=ax, top=True, right=True, left=True, bottom=True)
        plt.xticks("")
        plt.ylabel("", fontsize=y_label_fontsize)
        plt.xlabel("")
        #adjust x-axis label position 
        ax.xaxis.set_label_coords(0.5, -.2)

if Save_fig == "yes":
    File_out = "{}_{}_time_series.png".format(i_value,transport_of_interest_list_run[i_value])
    plt.savefig(os.path.join(Folder_out,File_out),dpi=600,bbox_inches='tight')
else:
    plt.show()
    
    