# -*- coding: utf-8 -*-

### Metabolic genes analysis ###

## Libraries ##

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

## Program input ##

#Save_fig = "yes"
Save_fig = "no"

last_plot = "no"

#analysis = "Glycolysis"
analysis = "TCA"


# figure dimensions #

fig_x = 17
fig_y = 6

# Resolution #
n_bins = 10

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

y_limit_upper = 1.5
y_limit_lower = -1.5

## Title coords ##

title_fontsize = 40

title_y_coord = 1
title_x_coord = 0.4


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


## Folders ##

Folder1 = "Data/Metabolic genes"

Folder2 = "Data/Normalized"
Folder3 = "Results/Differential expression analysis/gene lists LRT"

Folder4 = "Data/Pvalue data/Metabolic processes"

Folder_out = "Results/Metabolic processes/Gene of interest"
os.makedirs(Folder_out,exist_ok=True)

## Files ##

File1 = "Metabolic genes overview.xlsx"
File2 = "DEseq2 Normalized data.csv"
File3 = "LRT differential expressed gene list.txt"

File_out_1 = "Glycolysis.png"
File_out_2 = "TCA cycle.png"
File_out_3 = "Other.png"
File_out_4 = "Transporters.png"
File_out_5 = "{} data stacked.xlsx".format(analysis)


## Global variables ##

Sample_list = ['M1','M3','M6','M12','M18','M24']

## Load data ##

df_metabolic_overview = pd.read_excel(os.path.join(Folder1,File1))
df_data = pd.read_csv(os.path.join(Folder2,File2),sep=";",decimal=",").rename(columns=({"Unnamed: 0":"Ensembl ID"}))

## Subset data ##

df_Glycolysis = df_metabolic_overview[df_metabolic_overview['Category'] == 'Glycolysis']
df_TCA = df_metabolic_overview[df_metabolic_overview['Category'] == 'TCA cycle']
df_Other = df_metabolic_overview[df_metabolic_overview['Category'] == 'Other']
df_Transporters = df_metabolic_overview[df_metabolic_overview['Category'] == 'Transporters']

df_data_metabolic_glycolysis = df_data[df_data['Ensembl ID'].isin(df_Glycolysis['Rat Ensembl'])].copy()
df_data_metabolic_TCA = df_data[df_data['Ensembl ID'].isin(df_TCA['Rat Ensembl'])].copy()
df_data_metabolic_other = df_data[df_data['Ensembl ID'].isin(df_Other['Rat Ensembl'])].copy()
df_data_metabolic_Transporters = df_data[df_data['Ensembl ID'].isin(df_Transporters['Rat Ensembl'])].copy()

## create mean of the 3 biological replicate ##
for key in Sample_list:
    df_data_metabolic_glycolysis[key] = df_data_metabolic_glycolysis[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)
    df_data_metabolic_TCA[key] = df_data_metabolic_TCA[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)
    df_data_metabolic_other[key] = df_data_metabolic_other[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)
    df_data_metabolic_Transporters[key] = df_data_metabolic_Transporters[['Sample.{}.R1'.format(key),
                                    'Sample.{}.R2'.format(key),
                                    'Sample.{}.R3'.format(key)]].mean(axis=1)

if analysis == "Glycolysis":
    ## Reduce dataframe subsets ##
    
    # Glycolysis #
    
    df_data_metabolic_glycolysis_mean = df_data_metabolic_glycolysis[['Ensembl ID']+Sample_list].copy()
    df_data_metabolic_glycolysis_mean = df_data_metabolic_glycolysis_mean.set_index('Ensembl ID').T
    df_data_metabolic_glycolysis_zscored = df_data_metabolic_glycolysis_mean.apply(stats.zscore).T
    
    df_zscored_stacked_glycolysis = df_data_metabolic_glycolysis_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    
    # save stacked for statistical analysis #
    df_zscored_stacked_glycolysis.to_excel(os.path.join(Folder4,File_out_5),index=False)
    
    df_zscored_stacked_glycolysis_1_3 = df_zscored_stacked_glycolysis[(df_zscored_stacked_glycolysis['Months'] == 'M1') | (df_zscored_stacked_glycolysis['Months'] == 'M3')]
    df_zscored_stacked_glycolysis_3_6 = df_zscored_stacked_glycolysis[(df_zscored_stacked_glycolysis['Months'] == 'M3') | (df_zscored_stacked_glycolysis['Months'] == 'M6')]
    df_zscored_stacked_glycolysis_6_12 = df_zscored_stacked_glycolysis[(df_zscored_stacked_glycolysis['Months'] == 'M6') | (df_zscored_stacked_glycolysis['Months'] == 'M12')]
    df_zscored_stacked_glycolysis_12_18 = df_zscored_stacked_glycolysis[(df_zscored_stacked_glycolysis['Months'] == 'M12') | (df_zscored_stacked_glycolysis['Months'] == 'M18')]
    df_zscored_stacked_glycolysis_18_24 = df_zscored_stacked_glycolysis[(df_zscored_stacked_glycolysis['Months'] == 'M18') | (df_zscored_stacked_glycolysis['Months'] == 'M24')]
    
    data_stack_glycolysis = [df_zscored_stacked_glycolysis_1_3,df_zscored_stacked_glycolysis_3_6,df_zscored_stacked_glycolysis_6_12,df_zscored_stacked_glycolysis_12_18,df_zscored_stacked_glycolysis_18_24]
    
     ## confidence interval - for fill ##
     
    confidence_interval_glycolysis_M1 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M1']['Values']))
    confidence_interval_glycolysis_M3 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M3']['Values']))
    confidence_interval_glycolysis_M6 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M6']['Values']))
    confidence_interval_glycolysis_M12 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M12']['Values']))
    confidence_interval_glycolysis_M18 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M18']['Values']))
    confidence_interval_glycolysis_M24 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_glycolysis[df_zscored_stacked_glycolysis['Months'] == 'M24']['Values']))
    
    # Create lists of confidence intervals #
    confidence_list = [confidence_interval_glycolysis_M1,
                       confidence_interval_glycolysis_M3,
                       confidence_interval_glycolysis_M6,
                       confidence_interval_glycolysis_M12,
                       confidence_interval_glycolysis_M18,
                       confidence_interval_glycolysis_M24]
     
    ## Initiation of figure
    scatter_array = []
    line_array = []
    
    fig = plt.figure(figsize=(fig_x, fig_y),constrained_layout=True)
    
    gs = fig.add_gridspec(1,2,width_ratios=[1,0.3])
    
    ax1 = fig.add_subplot(gs[0, 0])
    
    
    for i in range(0,5,1):
        ax1 = sns.lineplot(data=data_stack_glycolysis[i], x="Months",color='black',
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
    if last_plot == "yes":
        ## For last plot ##
        ## Format plot axes ##
        ax1.tick_params(bottom=False)
        plt.ylim([y_limit_lower,y_limit_upper])
        plt.yticks([y_limit_lower,-1,-0.5,0,0.5,1,y_limit_upper],fontsize=y_axis_font_size)
        sns.despine(fig=fig, ax=ax1, top=True, right=True, left=False, bottom=True)
        plt.xticks(Sample_list,['M1','M3','M6','M12','M18','M24'],fontsize=x_axis_font_size)
        plt.ylabel("Z-score", fontsize=y_label_fontsize)
        plt.xlabel("", fontsize=x_label_fontsize)
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
    
    fig.suptitle("Glycolysis", fontsize=title_fontsize,ha='center',va='center',y=title_y_coord,x=title_x_coord)
    
    if Save_fig == "yes":
        plt.savefig(os.path.join(Folder_out,File_out_1),dpi=600,bbox_inches='tight')
    else:
        plt.show()

elif analysis == "TCA":
    ## Reduce dataframe subsets ##
    
    # TCA #
    
    df_data_metabolic_TCA_mean = df_data_metabolic_TCA[['Ensembl ID']+Sample_list].copy()
    df_data_metabolic_TCA_mean = df_data_metabolic_TCA_mean.set_index('Ensembl ID').T
    df_data_metabolic_TCA_zscored = df_data_metabolic_TCA_mean.apply(stats.zscore).T
    
    df_zscored_stacked_TCA = df_data_metabolic_TCA_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    
    # save stacked for statistical analysis #
    df_zscored_stacked_TCA.to_excel(os.path.join(Folder4,File_out_5),index=False)
    
    df_zscored_stacked_TCA_1_3 = df_zscored_stacked_TCA[(df_zscored_stacked_TCA['Months'] == 'M1') | (df_zscored_stacked_TCA['Months'] == 'M3')]
    df_zscored_stacked_TCA_3_6 = df_zscored_stacked_TCA[(df_zscored_stacked_TCA['Months'] == 'M3') | (df_zscored_stacked_TCA['Months'] == 'M6')]
    df_zscored_stacked_TCA_6_12 = df_zscored_stacked_TCA[(df_zscored_stacked_TCA['Months'] == 'M6') | (df_zscored_stacked_TCA['Months'] == 'M12')]
    df_zscored_stacked_TCA_12_18 = df_zscored_stacked_TCA[(df_zscored_stacked_TCA['Months'] == 'M12') | (df_zscored_stacked_TCA['Months'] == 'M18')]
    df_zscored_stacked_TCA_18_24 = df_zscored_stacked_TCA[(df_zscored_stacked_TCA['Months'] == 'M18') | (df_zscored_stacked_TCA['Months'] == 'M24')]
    
    data_stack_TCA = [df_zscored_stacked_TCA_1_3,df_zscored_stacked_TCA_3_6,df_zscored_stacked_TCA_6_12,df_zscored_stacked_TCA_12_18,df_zscored_stacked_TCA_18_24]
    
     ## confidence interval - for fill ##
     
    confidence_interval_TCA_M1 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M1']['Values']))
    confidence_interval_TCA_M3 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M3']['Values']))
    confidence_interval_TCA_M6 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M6']['Values']))
    confidence_interval_TCA_M12 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M12']['Values']))
    confidence_interval_TCA_M18 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M18']['Values']))
    confidence_interval_TCA_M24 = sns.utils.ci(sns.algorithms.bootstrap(df_zscored_stacked_TCA[df_zscored_stacked_TCA['Months'] == 'M24']['Values']))
    
    # Create lists of confidence intervals #
    confidence_list = [confidence_interval_TCA_M1,
                       confidence_interval_TCA_M3,
                       confidence_interval_TCA_M6,
                       confidence_interval_TCA_M12,
                       confidence_interval_TCA_M18,
                       confidence_interval_TCA_M24]
    # TCA #
    
    df_data_metabolic_TCA_mean = df_data_metabolic_TCA[['Ensembl ID']+Sample_list].copy()
    df_data_metabolic_TCA_mean = df_data_metabolic_TCA_mean.set_index('Ensembl ID')
    df_data_metabolic_TCA_mean = df_data_metabolic_TCA_mean.T
    df_data_metabolic_TCA_zscored_T = df_data_metabolic_TCA_mean.apply(stats.zscore)
    df_data_metabolic_TCA_zscored = df_data_metabolic_TCA_zscored_T.T
    
    
    # Others #
    
    df_data_metabolic_other_mean = df_data_metabolic_other[['Ensembl ID']+Sample_list].copy()
    df_data_metabolic_other_mean = df_data_metabolic_other_mean.set_index('Ensembl ID')
    df_data_metabolic_other_mean = df_data_metabolic_other_mean.T
    df_data_metabolic_other_zscored_T = df_data_metabolic_other_mean.apply(stats.zscore)
    df_data_metabolic_other_zscored = df_data_metabolic_other_zscored_T.T
    
    
    # Transporters  #
    df_data_metabolic_Transporters_mean = df_data_metabolic_Transporters[['Ensembl ID']+Sample_list].copy()
    df_data_metabolic_Transporters_mean = df_data_metabolic_Transporters_mean.set_index('Ensembl ID')
    df_data_metabolic_Transporters_mean = df_data_metabolic_Transporters_mean.T
    df_data_metabolic_Transporters_zscored_T = df_data_metabolic_Transporters_mean.apply(stats.zscore)
    df_data_metabolic_Transporters_zscored = df_data_metabolic_Transporters_zscored_T.T
    
    ## Stack data ##
    
    df_TCA_zscored_stacked = df_data_metabolic_TCA_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    df_TCA_zscored_stacked = df_data_metabolic_TCA_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    df_Other_zscored_stacked = df_data_metabolic_other_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    df_Transporters_zscored_stacked = df_data_metabolic_Transporters_zscored.stack().reset_index().rename(columns={"level_1":"Months",0:"Values"})
    
    
    ## Initiation of figure
    scatter_array = []
    line_array = []
    
    fig = plt.figure(figsize=(fig_x, fig_y),constrained_layout=True)
    
    gs = fig.add_gridspec(1,2,width_ratios=[1,0.3])
    
    ax1 = fig.add_subplot(gs[0, 0])
    
    
    for i in range(0,5,1):
        ax1 = sns.lineplot(data=data_stack_TCA[i], x="Months",color='black',
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
        plt.yticks([y_limit_lower,-1,-0.5,0,0.5,1,y_limit_upper],fontsize=y_axis_font_size)
        sns.despine(fig=fig, ax=ax1, top=True, right=True, left=False, bottom=True)
        plt.xticks(Sample_list,['M1','M3','M6','M12','M18','M24'],fontsize=x_axis_font_size)
        plt.ylabel("Z-score", fontsize=y_label_fontsize)
        plt.xlabel("", fontsize=x_label_fontsize)
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
    
    fig.suptitle("Tricarboxylic acid (TCA) cycle", fontsize=title_fontsize,ha='center',va='center',y=title_y_coord,x=title_x_coord)
    
    if Save_fig == "yes":
        plt.savefig(os.path.join(Folder_out,File_out_2),dpi=600,bbox_inches='tight')
    else:
        plt.show()