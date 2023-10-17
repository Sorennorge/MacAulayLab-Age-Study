# -*- coding: utf-8 -*-

### Molecular function time series enrichment ###

## libraries ##

import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt 
import seaborn as sns
from natsort import index_natsorted
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection


# input #

#sns.set_style("white")
plt.set_loglevel('WARNING') 

## Save figure ##

#Save_fig = "yes"
Save_fig = "no"

## Last plot ##

last_plot = "no"

## Add z-aksis ##

#z_aksis_print = "yes"
z_aksis_print = "no"

# Resolution #
n_bins = 10

# figure dimensions #

fig_x = 12
fig_y = 2.5

# outer line #
outer_line_width = 1

# outer points #
outer_point_size = 2

# Circles #
circle_edge_width = 2
font_size_circle = 14

# Line between points #
line_thickness = 6
line_edge_thickness = line_thickness+2

# Alpha settings #
Alpha_line = 1.0
alpha_fill = 0.1 # Needs changing depending on the resolution 
pie_alpha_all = 0.2
pie_alpha_focus = 1.0

## Font sizes ##
x_axis_font_size = 25
y_axis_font_size = 25
y_label_fontsize = 30
x_label_fontsize = 26

y_limit_upper = 2.1
y_limit_lower = -2

## Title coords ##

title_y_coord = 0.9
title_x_coord = 0.4


## Folders ##

Folder_1 = "Data/Panther/LRT diff/Enrichment data"
Folder_2 = "Data/Panther/Molecular functions"

Folder_3 = "Data/normalized"

Folder_4 = "Results/Time series/Molecular function/main titled"
os.makedirs(Folder_4,exist_ok=True)

Folder_5 = "Data/Pvalue data/MF"
os.makedirs(Folder_5,exist_ok=True)

## Files ##

File_1 = "MF_enrichment_data_all.xlsx"
File_2 = "Molecular function time series overview.xlsx"
File_3 = "DEseq2 Normalized data.csv"

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

# Global variables #

Sample_list = ['M1','M3','M6','M12','M18','M24']

## Load data ##

# Categories #

df_enrichment = pd.read_excel(os.path.join(Folder_1,File_1))

df_enrichment_subset = df_enrichment[df_enrichment['Percentage'] >= 1.0]

MF_list = df_enrichment_subset['Class'].tolist()
MF_list = MF_list[:-1]

df_overview = pd.read_excel(os.path.join(Folder_2,File_2))

    
df_overview = df_overview[['Ensembl ID','MF']]

df_overview = df_overview.set_index('Ensembl ID')

for i_value in MF_list:

    df_overview_subset = df_overview[df_overview['MF'] == i_value]

    ## load normalized data for each genes ##
    
    df_data = pd.read_csv(os.path.join(Folder_3,File_3),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))
    
    df_data['M1'] = df_data[['Sample.M1.R1','Sample.M1.R2','Sample.M1.R3']].mean(axis=1)
    df_data['M3'] = df_data[['Sample.M3.R1','Sample.M3.R2','Sample.M3.R3']].mean(axis=1)
    df_data['M6'] = df_data[['Sample.M6.R1','Sample.M6.R2','Sample.M6.R3']].mean(axis=1)
    df_data['M12'] = df_data[['Sample.M12.R1','Sample.M12.R2','Sample.M12.R3']].mean(axis=1)
    df_data['M18'] = df_data[['Sample.M18.R1','Sample.M18.R2','Sample.M18.R3']].mean(axis=1)
    df_data['M24'] = df_data[['Sample.M24.R1','Sample.M24.R2','Sample.M24.R3']].mean(axis=1)
    
    df_data_mean = df_data[['Ensembl ID',"M1","M3","M6","M12","M18","M24"]]
    df_data_mean = df_data_mean.set_index('Ensembl ID')
    
    df_overview_with_data = df_overview_subset.merge(df_data_mean, left_index=True, right_index=True)
    
    df_MF = df_overview_with_data.set_index("MF").T
    df_MF_zscored_T = df_MF.apply(stats.zscore)
    df_MF_zscored = df_MF_zscored_T.T
    
    df_MF_mean = df_overview_with_data.groupby('MF').mean()
    df_MF_mean_T = df_MF_mean.T
    
    
    Overview_file = "{} pvalue data.xlsx".format(i_value)
    df_MF_zscored.to_excel(os.path.join(Folder_5,Overview_file))
    
    df_MF_zscored_stacked = df_MF_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    df_MF_zscored_stacked = df_MF_zscored_stacked.sort_values(
        by="MF",
        key=lambda x: np.argsort(index_natsorted(df_MF_zscored_stacked["MF"])))
    
    print(i_value)
    ## Stack zscored values ##
    df_zscored_stacked = df_MF_zscored_stacked[df_MF_zscored_stacked.MF == i_value].copy()
    
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
    
    fig = plt.figure(figsize=(fig_x, fig_y),constrained_layout=True)
    
    gs = fig.add_gridspec(1,2,width_ratios=[1,0.3])
    
    ax1 = fig.add_subplot(gs[0, 0])
    
    for i in range(0,5,1):
        ax1 = sns.lineplot(data=data_stack[i], x="Months",color='black',
                           y="Values",alpha=0.75,
                           marker="o",
                           dashes=False,legend=None, markersize=10,errorbar=None)
        x, y = ax1.get_lines()[0].get_data()
        line_array.append([[x[0],x[1]],[y[0],y[1]]])
        ax1.get_lines()[0].remove()
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
        ax1.plot([x1,x2],[y1,y2],lw=line_edge_thickness,color='Black')
        ax1.plot([x1,x2],[y1,y2],lw=line_thickness,color='White')
    
    for i in range(0,6,1):
        if i == 0:
            # Plot upper line #
            ax1.plot([i,i+1],[confidence_list[i][1],confidence_list[i+1][1]],lw=outer_line_width,color='Black')
            # Plot lower line #
            ax1.plot([i,i+1],[confidence_list[i][0],confidence_list[i+1][0]],lw=outer_line_width,color='Black')
            #Plot begining line #
            ax1.plot([i,i],[confidence_list[i][0],confidence_list[i][1]],lw=outer_line_width,color='Black')    
        elif i == 5:
            #Plot end line #
            ax1.plot([i,i],[confidence_list[i][0],confidence_list[i][1]],lw=outer_line_width,color='Black')
        else:
            # Plot upper line #
            ax1.plot([i,i+1],[confidence_list[i][1],confidence_list[i+1][1]],lw=outer_line_width,color='Black')
            # Plot lower line #
            ax1.plot([i,i+1],[confidence_list[i][0],confidence_list[i+1][0]],lw=outer_line_width,color='Black')
    
    for i in range(0,5,1):
        # Segment lines into "n" bits #
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
        line = ax1.add_collection(lc)
      
    for i in range(0,len(scatter_array),1):
        point_x = scatter_array[i][0]
        point_y = scatter_array[i][1]
        ax1.plot(point_x, point_y, marker='o',color=RGB_list[i],markeredgecolor='black',markeredgewidth=circle_edge_width, ms=25)
    
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
                ax1.fill_between([x_seg[ii], 
                                 x_seg[ii+1]],
                                [y_lower[ii],
                                 y_lower[ii+1]],
                                [y_upper[ii],
                                 y_upper[ii+1]],
                                color=cmap(cols[ii]),
                                alpha=alpha_fill)
    
    
    if last_plot == "yes":
        ## For last plot ##
        ## Format plot axes ##
        ax1.tick_params(bottom=False)
        plt.ylim([y_limit_lower,y_limit_upper])
        plt.yticks([],fontsize=y_axis_font_size)
        sns.despine(fig=fig, ax=ax1, top=True, right=True, left=True, bottom=True)
        plt.xticks(Sample_list,['M1','M3','M6','M12','M18','M24'],fontsize=x_axis_font_size)
        plt.ylabel("", fontsize=y_label_fontsize)
        plt.xlabel("", fontsize=x_label_fontsize)
        #adjust x-axis label position 
        ax1.xaxis.set_label_coords(0.5, -.2)
    else:
		# If z-aksis is wanted this if-statement handles #
        if z_aksis_print == "yes":
            ax1.tick_params(bottom=False)
            ax1.tick_params(left=True)
            plt.ylim([y_limit_lower,y_limit_upper])
            plt.yticks([0,1],fontsize=y_axis_font_size)
            sns.despine(fig=fig, ax=ax1, top=True, right=True, left=False, bottom=True)
            ax1.spines['left'].set_bounds(0, 1)
            ax1.set_yticks([0,1])
            plt.xticks("")
            plt.ylabel("Δ z-score", fontsize=y_label_fontsize,y=0.76)
            plt.xlabel("")
            #adjust x-axis label position 
            ax1.xaxis.set_label_coords(0.5, -.2)
        else:
            ax1.tick_params(bottom=False)
            ax1.tick_params(left=False)
            plt.ylim([y_limit_lower,y_limit_upper])
            plt.yticks([],fontsize=y_axis_font_size)
            sns.despine(fig=fig, ax=ax1, top=True, right=True, left=True, bottom=True)
            plt.xticks("")
            plt.ylabel("", fontsize=y_label_fontsize)
            plt.xlabel("")
            #adjust x-axis label position 
            ax1.xaxis.set_label_coords(0.5, -.2)
    
    ### Subplot ###
    
    ## Variables ##
    
    data_all = df_enrichment['Count']
    labels_all = df_enrichment['Class']
    pal = sns.color_palette("crest_r",n_colors=12)
    
    label_index = list(labels_all).index(i_value)
    
    explode_set_all = [0.0]*len(labels_all)
    explode_set_all[label_index] = 0.15
    
    for n in range(0,len(pal),1):
        if n == label_index:
            pass
        else:
            pal[n] = "lightgrey"
    
    ax2 = fig.add_subplot(gs[0, 1])
    n = plt.pie(data_all,
            radius=1.0,
            explode=explode_set_all,
            colors  = pal[0:],
            wedgeprops = {"edgecolor":"black",'linewidth': 1,"alpha": 0.85},
            counterclock=False,
            startangle=45)
    ## set alpha ##
    
    hole = plt.Circle((0, 0), 0.6, facecolor='white',edgecolor=((0,0,0,pie_alpha_all)), linewidth=1,)
    plt.gcf().gca().add_artist(hole)
    
    for i in range(len(n[0])):
        if i == label_index:
            pass
        else:
            n[0][i].set_alpha(pie_alpha_all)
    ## Set title #
    fig.suptitle(i_value, fontsize=30,ha='center',va='center',y=title_y_coord,x=title_x_coord)

    if Save_fig == "yes":
        File_out = "{}_{}_time_series.png".format(label_index,i_value)
        plt.savefig(os.path.join(Folder_4,File_out),dpi=600,bbox_inches='tight',transparent=True)
    else:
        plt.show()