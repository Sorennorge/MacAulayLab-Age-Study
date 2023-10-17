# -*- coding: utf-8 -*-

### Enrichment plot - Biological function ###

## Libraries ##

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

## functions ##

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{:.1f}%\n({v:d})'.format(pct, v=val)
    return my_format

# Point connection function #

def draw_connection_line(placement_orientation,pie_chart_number, number_of_breaks, placement_number, length_start, length_orientation,explode_length):
    
    index_pie_chart_number = pie_chart_number-1
    
    length = length_start
    exploded_length = explode_length
    ang = numpy_array_list[index_pie_chart_number][0]

    x_start = numpy_array_list[index_pie_chart_number][1] + exploded_length * np.cos(np.deg2rad(ang))
    y_start = numpy_array_list[index_pie_chart_number][2] + exploded_length * np.sin(np.deg2rad(ang))

    x_end = x_start + length * np.cos(np.deg2rad(ang))
    y_end = y_start + length * np.sin(np.deg2rad(ang))
    
    # Placement_orientation needs to be right, or left
    if placement_orientation == 'right':
        length_right = length_orientation
        # number_of_breaks needs to be 1-3 #
        if number_of_breaks == 1:
            # Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
            while_placement = placement_number
            while_placement_break_counter = 0  # ensures no infinity loops
            if while_placement < 0:
                while (y_end > while_placement+placement_correction_unit) or (y_end < while_placement-placement_correction_unit):
                    if y_end > while_placement:
                        length = length+placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    elif y_end < while_placement:
                        length = length-placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    else:
                        print(while_placement,y_end,length)
                    # ensures no infinity loops #
                    while_placement_break_counter += 1 
                    if while_placement_break_counter > 100000:
                        break
            else:
                while (y_end > while_placement+placement_correction_unit) or (y_end < while_placement-placement_correction_unit):
                    if y_end < while_placement:
                        length = length+placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    elif y_end > while_placement:
                        length = length-placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    else:
                        print(while_placement,y_end,length)
                    # ensures no infinity loops
                    while_placement_break_counter += 1 
                    if while_placement_break_counter > 100000:
                        break
            # Plot 2 lines (one break)
            plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([x_end,length_right-0.05], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
            
        elif number_of_breaks == 2:
            plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([x_end,length_right-0.125], [y_end,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([length_right-0.125,length_right-0.05], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
        else:
            pass
    elif placement_orientation == 'left':
        length_left = length_orientation
        # Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
        if number_of_breaks == 1:
            while_placement = placement_number
            if while_placement < 0:
                while (y_end > while_placement+placement_correction_unit) or (y_end < while_placement-placement_correction_unit):
                    if y_end > while_placement:
                        length = length+placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    elif y_end < while_placement:
                        length = length-placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    else:
                        print(while_placement,y_end,length)
            else:
                while (y_end > while_placement+placement_correction_unit) or (y_end < while_placement-placement_correction_unit):
                    if y_end < while_placement:
                        length = length+placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    elif y_end > while_placement:
                        length = length-placement_correction_unit
                        x_end = x_start + length * np.cos(np.deg2rad(ang))
                        y_end = y_start + length * np.sin(np.deg2rad(ang))
                    else:
                        print(while_placement,y_end,length)
            # Plot 2 lines (one break)
            plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([x_end,length_left+0.04], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
        elif number_of_breaks == 2:
            plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([x_end,length_left+0.2], [y_end,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([length_left+0.2,length_left+0.04], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)

        elif number_of_breaks == 3:
            # Four points
            plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([x_end,length_left+0.3], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([length_left+0.3,length_left+0.2], [y_end,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
            plt.plot([length_left+0.2,length_left+0.04], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)

        else:
            pass
    else:
        pass
    
def draw_connection_line_small_group(placement_orientation,pie_chart_number_start,pie_chart_number_end, placement_number,length_start,collection_point):
    # set global #
    small_group_collection_point = collection_point
    length = length_start
    # for start end set collection point #
    for i in range(pie_chart_number_start,pie_chart_number_end+1,1):
        index_pie_chart_number = i-1
        
        ang = numpy_array_list[index_pie_chart_number][0]
    
        x_start = numpy_array_list[index_pie_chart_number][1]
        y_start = numpy_array_list[index_pie_chart_number][2]
    
        x_end = x_start + length * np.cos(np.deg2rad(ang))
        y_end = y_start + length * np.sin(np.deg2rad(ang))
        
        plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
        plt.plot([x_end,length_left+small_group_collection_point], [y_end,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)
    # Draw end line
    plt.plot([length_left+small_group_collection_point,length_left+0.04], [placement_number,placement_number], ls='dashed',color="grey", linewidth=2,clip_on = False)


## Global variables ##

save_file = "yes"
#save_file = "no"

font_size_text = 34
fig_x = 12
fig_y = 10

## Folders ##

Folder1 = "Data/Panther/LRT diff/Enrichment data"
Folder2 = "Results/Enrichment/Biological functions"
os.makedirs(Folder2,exist_ok=True)


## Files ##

File1 = "BF_enrichment_data_all.xlsx"
File3 = "BF_Enrichment_all.png"

## Load data ##

df_all = pd.read_excel(os.path.join(Folder1,File1))

## Variables ##

data_all = df_all['Count']
labels_all = df_all['Class']
percent_list = df_all['Percentage'].tolist()

others_percentage = 0

for key in percent_list:
    key_value = float(key)
    if key_value < 1.0:
        others_percentage += key_value
others_percentage = round(others_percentage,1)

## Set colors ##

pal = sns.color_palette("crest_r",n_colors=21)

explode_set_all = [0]*len(labels_all)
explode_set_all[2] = 0.1


## Enrichment plot - classified ##
print("Creating enrichment plot for all metabolites")
plt.figure(figsize=(fig_x,fig_y))
wedges, texts = plt.pie(data_all,
        radius=1.0,
        explode=explode_set_all,
        colors  = pal[0:],
        wedgeprops = {"edgecolor":"black",'linewidth': 1,"alpha": 0.85},
        counterclock=False,
        startangle=90)


### Add a hole in the pie  
hole = plt.Circle((0, 0), 0.6, facecolor='white',edgecolor="black", linewidth=1)
plt.gcf().gca().add_artist(hole)

# Protein classes
plt.text(0,0,"Biological\nfunction",fontsize=44,ha='center',va='center',fontweight="bold")

# Length of left #
length_right = 1.3
length_left = -1.5

# Placement positions #

Placement_right_1 = float(0.80)
Placement_right_2 = float(-0.3)
Placement_right_3 = float(-1.4)

Placement_left_1 = float(-1.4)
Placement_left_2 = float(-1.15)
Placement_left_3 = float(-0.9)
Placement_left_4 = float(-0.65)
Placement_left_5 = float(-0.4)
Placement_left_6 = float(-0.15)
Placement_left_7 = float(0.1)
Placement_left_8 = float(0.35)
Placement_left_9 = float(0.60)
Placement_left_10 = float(0.9)

placement_correction_unit = 0.0001

# Cellular process #
Cellular_label = labels_all[0].replace("lar proc","lar\nproc")
plt.text(length_right,Placement_right_1,"{} ({}%)".format(Cellular_label,percent_list[0]),fontsize=font_size_text,ha='left',va='center')
# Biological regulation #
regulation_label = labels_all[1].replace("gical reg","gical\nreg")
plt.text(length_right,Placement_right_2,"{} ({}%)".format(regulation_label,percent_list[1]),fontsize=font_size_text,ha='left',va='center')
# Metabolic process #
Metabolic_label = labels_all[2].replace("lic proc","lic\nproc")
plt.text(length_right,Placement_right_3,"{} ({}%)".format(Metabolic_label,percent_list[2]),fontsize=font_size_text,weight='bold',ha='left',va='center')

### Break left ##

# Response to stimulus #
plt.text(length_left,Placement_left_1,"{} ({}%)".format(labels_all[3],percent_list[3]),fontsize=font_size_text,ha='right',va='center')

# Localization #
plt.text(length_left,Placement_left_2,"{} ({}%)".format(labels_all[4],percent_list[4]),fontsize=font_size_text,ha='right',va='center')


# Signaling #
plt.text(length_left,Placement_left_3,"{} ({}%)".format(labels_all[5],percent_list[5]),fontsize=font_size_text,ha='right',va='center')

# Developmental process #
plt.text(length_left,Placement_left_4,"{} ({}%)".format(labels_all[6],percent_list[6]),fontsize=font_size_text,ha='right',va='center')

# Multicellular organismal process #
plt.text(length_left,Placement_left_5,"{} ({}%)".format(labels_all[7],percent_list[7]),fontsize=font_size_text,ha='right',va='center')

# Biological adhesion #
plt.text(length_left,Placement_left_6,"{} ({}%)".format(labels_all[8],percent_list[8]),fontsize=font_size_text,ha='right',va='center')

# Locomotion #
plt.text(length_left,Placement_left_7,"{} ({}%)".format(labels_all[9],percent_list[9]),fontsize=font_size_text,ha='right',va='center')

# Immune system process #
plt.text(length_left,Placement_left_8,"{} ({}%)".format(labels_all[10],percent_list[10]),fontsize=font_size_text,ha='right',va='center')

# Small group collection #

plt.text(length_left,Placement_left_9,"Small group collection ({}%)".format(others_percentage),fontsize=font_size_text,ha='right',va='center')


# Unclassified #

plt.text(length_left,Placement_left_10,"Unclassified ({}%)".format(percent_list[-1]),fontsize=font_size_text,ha='right',va='center')


### Draw connecting lines ###


# Cellular process #
# Biological regulation #
# Metabolic process #
# Response to stimulus #
# Localization #
# Signaling #
# Developmental process #
# Multicellular organismal process #
# Biological adhesion #
# Locomotion #
# Immune system process #
# Reproductive process #
# Reproduction #
# Biological process involved in interspecies interaction between organisms #
# Growth #
# Biological phase #
# Rhythmic process #
# Biomineralization #
# Pigmentation #
# Unclassified #

# Calculate degress #
numpy_array_list = []
for i, p in enumerate(wedges):
    if i == 23:
        ang = (p.theta2 - p.theta1)/1.25 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        numpy_array_list.append([ang,x,y]) 
    else:
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        numpy_array_list.append([ang,x,y])

# Cellular process #

# Length of left #
#length_right = 1.3
#length_left = -1.5

draw_connection_line("right",1, 1, Placement_right_1, 3.1, length_right,0)
# Biological regulation #
draw_connection_line("right",2, 1, Placement_right_2, 0.1, length_right,0)
# Metabolic process #
draw_connection_line("right",3, 1, Placement_right_3, 0.1, length_right,0.1)
# Response to stimulus #
draw_connection_line("left",4, 1, Placement_left_1, 0.2, length_left,0)
# Localization #
draw_connection_line("left",5, 1, Placement_left_2, 0.2, length_left,0)
# Signaling #
draw_connection_line("left",6, 1, Placement_left_3, 0.2, length_left,0)
# Developmental process #
draw_connection_line("left",7, 3, Placement_left_4, 0.1, length_left,0)
# Multicellular organismal process #
draw_connection_line("left",8, 3, Placement_left_5, 0.1, length_left,0)
# Biological adhesion #
draw_connection_line("left",9, 2, Placement_left_6, 0.2, length_left,0)
# Locomotion #
draw_connection_line("left",10, 2, Placement_left_7, 0.2, length_left,0)
# Immune system process #
draw_connection_line("left",11, 2, Placement_left_8, 0.2, length_left,0)

# Small group collection #
small_group_length = 0.10
small_group_collection_point = 0.20

draw_connection_line_small_group("left",12,19, Placement_left_9,small_group_length,small_group_collection_point)

# Unclassified #
draw_connection_line("left",20, 1, Placement_left_10, 0.2, length_left,0)


if save_file == 'yes':
    print("Saving plot...")
    plt.savefig(os.path.join(Folder2,File3),dpi=600,bbox_inches='tight')
    print("Done.")
else:
    plt.show()
