# -*- coding: utf-8 -*-

### Toxin removal analysis ###

import os
import pandas as pd
from scipy import stats
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection
from natsort import natsorted
from textwrap import wrap

## Program input ##

# Save figure #
#Save_fig = "yes"
Save_fig = "no"

# figure dimensions #
fig_x = 12
fig_y = 4

group_start = 1
groups = 11 # For all set to 11

# Resolution #
n_bins = 10

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
pie_alpha_all = 0.2
pie_alpha_focus = 1.0

## Font sizes ##
x_axis_font_size = 24
y_axis_font_size = 24
y_label_fontsize = 26
x_label_fontsize = 26

title_label_fontsize = 28
legend_label_fontsize = 26

y_limit_upper = 2
y_limit_lower = -2

label_sep_value = 0.075

# Title adjustment #
title_length = 35
X_adjust = 0


## Colormapping ##

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


## Folders ##

Folder1 = "Data/Toxin removal"
Folder2 = "Data/Normalized"

Folder3 = "Results/Toxin removal"
os.makedirs(Folder3,exist_ok=True)

Folder4 = "Data/Pvalue data/Toxin removal"
os.makedirs(Folder4,exist_ok=True)

## Files ##

File1 = "Toxin_removal_Gene_overview.xlsx"
File2 = "DEseq2 Normalized data.csv"

## Global variables ##

Sample_list = ['M1','M3','M6','M12','M18','M24']


## Load data ##

df_toxin_removal_overview = pd.read_excel(os.path.join(Folder1,File1))
df_data = pd.read_csv(os.path.join(Folder2,File2),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))


## plot titles ##

df_toxin_removal_overview_titles = df_toxin_removal_overview[['Group','Description']]
df_toxin_removal_overview_titles = df_toxin_removal_overview_titles.drop_duplicates()
title_dict = df_toxin_removal_overview_titles.set_index('Group')['Description'].to_dict()

for i in range(group_start,groups,1):
    i_value = i
    df_Group = df_toxin_removal_overview[df_toxin_removal_overview['Group'] == 'Group {}'.format(i_value)]
    df_data_removal = df_data[df_data['Ensembl ID'].isin(df_Group['Ensembl ID'])].copy()
    for key in Sample_list:
        df_data_removal[key] = df_data_removal[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)
        
    df_data_removal_mean = df_data_removal[['Ensembl ID']+Sample_list].copy()
    df_data_removal_mean = df_data_removal_mean.set_index('Ensembl ID').T
    df_data_removal_zscored = df_data_removal_mean.apply(stats.zscore).T
    
    df_zscored_stacked = df_data_removal_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    
    # Save data for statistical analysis #
    name = title_dict['Group {}'.format(i_value)]
    temp_file = "{} data stacked.xlsx".format(name.replace('/',''))
    if i == 1:
        file_data_names = pd.DataFrame(data=[{'File names':temp_file}],index=[i_value])
    else:
        df_new_row = pd.DataFrame(data=[{'File names':temp_file}],index=[i_value])
        file_data_names = pd.concat([file_data_names, df_new_row])
    df_zscored_stacked.to_excel(os.path.join(Folder4,temp_file),index=False)
    file_data_names.to_excel(os.path.join(Folder4,"File_list.xlsx"))
    
    # stack zscore #
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
        #print(i)
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
    
    #  #
    ax1.tick_params(bottom=False)
    ax1.tick_params(left=True)
    plt.ylim([y_limit_lower,y_limit_upper])
    plt.yticks([-2,-1,0,1,2],fontsize=y_axis_font_size)
    sns.despine(fig=fig, ax=ax1, top=True, right=True, left=False, bottom=True)
    plt.xticks("")
    plt.ylabel("Z-score", fontsize=y_label_fontsize)
    plt.xlabel("")
    #adjust x-axis label position 
    ax1.xaxis.set_label_coords(0.5, -.2)

    # Title # 
    ### in case of group 9, the title will be split in three lines, and to avoid this a brute force if statement is implemented. ###
    if i_value == 9:    
        if len(title_dict['Group {}'.format(i_value)]) > 37:
            plt.title("{}".format('\n'.join(wrap(title_dict['Group {}'.format(i_value)], 37))),fontsize=title_label_fontsize,y=1.05,x=0.5-X_adjust,wrap=True,ha='center')
        else:
            plt.title("\n{}".format('\n'.join(wrap(title_dict['Group {}'.format(i_value)], 37))),fontsize=title_label_fontsize,y=1.05,x=0.5-X_adjust,wrap=True,ha='center')
    else:
        if len(title_dict['Group {}'.format(i_value)]) > 37:
            plt.title("{}".format('\n'.join(wrap(title_dict['Group {}'.format(i_value)], 37))),fontsize=title_label_fontsize,y=1.05,x=0.5-X_adjust,wrap=True,ha='center')
        else:
            plt.title("\n{}".format('\n'.join(wrap(title_dict['Group {}'.format(i_value)], 37))),fontsize=title_label_fontsize,y=1.05,x=0.5-X_adjust,wrap=True,ha='center')
    # save figure #
    if Save_fig == "yes":
        File_out = "Group_{}_Figure.png".format(i_value)
        plt.savefig(os.path.join(Folder3,File_out),dpi=600,bbox_inches='tight')
    else:
        plt.show()
    
## If wanted - right label for all the genes in the different analyses ##
for i in range(group_start,groups,1):
    i_value = i
    
    df_Group = df_toxin_removal_overview[df_toxin_removal_overview['Group'] == 'Group {}'.format(i_value)]
        
    fig = plt.figure(figsize=(int(fig_x*0.4), fig_y),constrained_layout=True)
    
    gs = fig.add_gridspec(1,1,width_ratios=[1])
    
    ax1 = fig.add_subplot(gs[0, 0])
    
    ## Text plots ##
    ax1.spines[['right', 'top','left','bottom']].set_visible(False)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_ylim(-0.2,1.2)
    ax1.set_xlim(-0.1,1)
    Gene_list = df_Group['Gene'].tolist()
    Gene_list_sorted = natsorted(Gene_list)
    dot_dist = 0.1

    # evaluate length of list #

    # If length of list is uneven #
    if len(Gene_list_sorted) % 2 == 1:
        start_index = int(len(Gene_list_sorted)/2)
        start_point = 0.5
        
        for n in range(start_index,-1,-1):
            ax1.scatter([0.0],[start_point],color='black')
            ax1.text(0.0+dot_dist,start_point,"{}".format(Gene_list_sorted[n]),fontsize=title_label_fontsize,horizontalalignment='left',
                verticalalignment='center_baseline')
            start_point = start_point+label_sep_value
        start_point = 0.5
        for n in range(start_index+1,len(Gene_list_sorted),1):
            start_point = start_point-label_sep_value
            ax1.scatter([0.0],[start_point],color='black')
            ax1.text(0.0+dot_dist,start_point,"{}".format(Gene_list_sorted[n]),fontsize=title_label_fontsize,horizontalalignment='left',
                verticalalignment='center_baseline')
    else:
        # If length of list is even #
        first_startpoint = int(len(Gene_list_sorted) / 2) - 1
        second_startpoint = int(len(Gene_list_sorted) / 2)
        start_point_first = 0.5+(label_sep_value/15)
        for n in range(first_startpoint,-1,-1):
            ax1.scatter([0.0],[start_point_first],color='black')
            ax1.text(0.0+dot_dist,start_point_first,"{}".format(Gene_list_sorted[n]),fontsize=title_label_fontsize,horizontalalignment='left',
                verticalalignment='center_baseline')
            start_point_first = start_point_first+label_sep_value
        start_point_second = 0.5-(label_sep_value/15)
        for n in range(second_startpoint,len(Gene_list_sorted),1):
            start_point_second = start_point_second-label_sep_value
            ax1.scatter([0.0],[start_point_second],color='black')
            ax1.text(0.0+dot_dist,start_point_second,"{}".format(Gene_list_sorted[n]),fontsize=title_label_fontsize,horizontalalignment='left',
                verticalalignment='center_baseline')
    if Save_fig == "yes":
        File_out_legend = "Group_{}_Legend.png".format(i_value)
        plt.savefig(os.path.join(Folder3,File_out_legend),dpi=600,bbox_inches='tight')
    else:
        plt.show()