# -*- coding: utf-8 -*-

### Enrichment plot - Protein classes ###

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

## Global variables - Program settings ##

save_file = "yes"
#save_file = "no"

font_size_text = 34
fig_x = 12
fig_y = 10

## Folders ##

Folder1 = "Data/Panther/LRT diff/Enrichment data"
Folder2 = "Results/Enrichment/Protein classes"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "PC_enrichment_data_all.xlsx"

File3 = "Enrichment_all.png"

## Load data ##

df_all = pd.read_excel(os.path.join(Folder1,File1))

## Variables ##

data_all = df_all['Count']
labels_all = df_all['Class']
percent_list = df_all['Percentage'].tolist()

others_percentage = 0

for key in percent_list:
    key_value = float(key)
    if key_value < 2.0:
        others_percentage += key_value
others_percentage = round(others_percentage,1)

## Set colors ##

pal = sns.color_palette("crest_r",n_colors=30)

## Enrichment plot - classified ##
print("Creating enrichment plot for all metabolites")
plt.figure(figsize=(fig_x,fig_y))
wedges, texts = plt.pie(data_all,
        radius=1.0,
        colors  = pal[0:],
        wedgeprops = {"edgecolor":"black",'linewidth': 1,"alpha": 0.85},
        counterclock=False,
        startangle=60)


### Add a hole in the pie  
hole = plt.Circle((0, 0), 0.6, facecolor='white',edgecolor="black", linewidth=1)
plt.gcf().gca().add_artist(hole)

# Protein classes
plt.text(0,0,"Protein\nclasses",fontsize=44,ha='center',va='center',fontweight="bold")

# Length of left #
length_right = 1.3
length_left = -1.5

# Placement positions #

Placement_right_1 = float(0.9)
Placement_right_2 = float(0.5)
Placement_right_3 = float(0.15)
Placement_right_4 = float(-0.3)
Placement_right_5 = float(-0.75)
Placement_right_6 = float(-1.1)
Placement_right_7 = float(-1.4)

Placement_left_1 = float(-1.4)
Placement_left_2 = float(-1.15)
Placement_left_3 = float(-0.9)
Placement_left_4 = float(-0.65)
Placement_left_5 = float(-0.4)
Placement_left_6 = float(-0.15)
Placement_left_7 = float(0.1)
Placement_left_8 = float(0.4)
Placement_left_9 = float(0.9)

placement_correction_unit = 0.0001

# Metabolite #
Metabolite_label = labels_all[0].replace("sion enzy","sion\nenzy")
plt.text(length_right,Placement_right_1,"{} ({}%)".format(Metabolite_label,percent_list[0]),fontsize=font_size_text,ha='left',va='center')
# protein modifying #
plt.text(length_right,Placement_right_2,"{} ({}%)".format(labels_all[1],percent_list[1]),fontsize=font_size_text,ha='left',va='center')
# Transporter #
plt.text(length_right,Placement_right_3,"{} ({}%)".format(labels_all[2],percent_list[2]),fontsize=font_size_text,ha='left',va='center')
# Gene specific #
protein_binding = labels_all[3].replace("vity modu","vity\nmodu")
plt.text(length_right,Placement_right_4,"{} ({}%)".format(protein_binding,percent_list[3]),fontsize=font_size_text,ha='left',va='center')
# protein binding #
gene_specific = labels_all[4].replace("nal reg","al\nreg")
plt.text(length_right,Placement_right_5,"{} ({}%)".format(gene_specific,percent_list[4]),fontsize=font_size_text,ha='left',va='center')
# Cytoskeletal protein #
plt.text(length_right,Placement_right_6,"{} ({}%)".format(labels_all[5],percent_list[5]),fontsize=font_size_text,ha='left',va='center')
# scaffold #
plt.text(length_right,Placement_right_7,"{} ({}%)".format(labels_all[6],percent_list[6]),fontsize=font_size_text,ha='left',va='center')

## Break left ##

# Transmembrane #
plt.text(length_left,Placement_left_1,"{} ({}%)".format(labels_all[7],percent_list[7]),fontsize=font_size_text,ha='right',va='center')
# cell adhesion #
plt.text(length_left,Placement_left_2,"{} ({}%)".format(labels_all[8],percent_list[8]),fontsize=font_size_text,ha='right',va='center')
# Intercellular signal #
plt.text(length_left,Placement_left_3,"{} ({}%)".format(labels_all[9],percent_list[9]),fontsize=font_size_text,ha='right',va='center')
# Translation #
plt.text(length_left,Placement_left_4,"{} ({}%)".format(labels_all[10],percent_list[10]),fontsize=font_size_text,ha='right',va='center')

# membrane traffic #
plt.text(length_left,Placement_left_5,"{} ({}%)".format(labels_all[11],percent_list[11]),fontsize=font_size_text,ha='right',va='center')

# Extracelllular #
plt.text(length_left,Placement_left_6,"{} ({}%)".format(labels_all[12],percent_list[12]),fontsize=font_size_text,ha='right',va='center')

# RNA metabolism #

chrom_label = labels_all[13].replace(",",",\n")
plt.text(length_left,Placement_left_7,"{} ({}%)".format(chrom_label,percent_list[13]),fontsize=font_size_text,ha='right',va='center')

# Small group collection #

plt.text(length_left,Placement_left_8,"Small group collection ({}%)".format(others_percentage),fontsize=font_size_text,ha='right',va='center')

# Unclassified #
plt.text(length_left+0.02,Placement_left_9,"Unclassified\n({})".format(percent_list[-1]),fontsize=font_size_text,ha='right',va='center')


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


# Metabolite #
length = 0.2
x_start = numpy_array_list[0][1]
y_start = numpy_array_list[0][2]
ang = numpy_array_list[0][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))

# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_right_1
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
    
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.05], [Placement_right_1,Placement_right_1], ls='dashed',color="grey", linewidth=2,clip_on = False)

## Label and line placement ##

# protein modifying #

length = 0.075
x_start = numpy_array_list[1][1]
y_start = numpy_array_list[1][2]
ang = numpy_array_list[1][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.125], [y_end,Placement_right_2], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.125,length_right-0.05], [Placement_right_2,Placement_right_2], ls='dashed',color="grey", linewidth=2,clip_on = False)


# Transporter #

length = 0.075
x_start = numpy_array_list[2][1]
y_start = numpy_array_list[2][2]
ang = numpy_array_list[2][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.25], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.25,length_right-0.125], [y_end,Placement_right_3], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.125,length_right-0.05], [Placement_right_3,Placement_right_3], ls='dashed',color="grey", linewidth=2,clip_on = False)

# protein binding #

length = 0.15
x_start = numpy_array_list[3][1]
y_start = numpy_array_list[3][2]
ang = numpy_array_list[3][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.275], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.275,length_right-0.125], [y_end,Placement_right_4], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.125,length_right-0.05], [Placement_right_4,Placement_right_4], ls='dashed',color="grey", linewidth=2,clip_on = False)

# Gene specific #

length = 0.15
x_start = numpy_array_list[4][1]
y_start = numpy_array_list[4][2]
ang = numpy_array_list[4][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Four points #
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.2], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.2,length_right-0.125], [y_end,Placement_right_5], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_right-0.125,length_right-0.05], [Placement_right_5,Placement_right_5], ls='dashed',color="grey", linewidth=2,clip_on = False)

# Cytoskeletal protein #

length = 0.01
x_start = numpy_array_list[5][1]
y_start = numpy_array_list[5][2]
ang = numpy_array_list[5][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))

# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_right_6
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
            
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.04], [y_end,Placement_right_6], ls='dashed',color="grey", linewidth=2,clip_on = False)

# scaffold #

length = 0.55
x_start = numpy_array_list[6][1]
y_start = numpy_array_list[6][2]
ang = numpy_array_list[6][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))

# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_right_7
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
            
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_right-0.04], [y_end,Placement_right_7], ls='dashed',color="grey", linewidth=2,clip_on = False)



## Break left ##

# Transmembrane #

length = 0.1
x_start = numpy_array_list[7][1]
y_start = numpy_array_list[7][2]
ang = numpy_array_list[7][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))

# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_left_1
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
            
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.04], [y_end,Placement_left_1], ls='dashed',color="grey", linewidth=2,clip_on = False)

# Cell adhesion #
length = 0.1
x_start = numpy_array_list[8][1]
y_start = numpy_array_list[8][2]
ang = numpy_array_list[8][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))

# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_left_2
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
            
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.04], [y_end,Placement_left_2], ls='dashed',color="grey", linewidth=2,clip_on = False)



# Intercelluar signal #
length = 0.20
x_start = numpy_array_list[9][1]
y_start = numpy_array_list[9][2]
ang = numpy_array_list[9][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Four points
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.3], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.3,length_left+0.2], [y_end,Placement_left_3], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.2,length_left+0.04], [Placement_left_3,Placement_left_3], ls='dashed',color="grey", linewidth=2,clip_on = False)



# Translational #
length = 0.2
x_start = numpy_array_list[10][1]
y_start = numpy_array_list[10][2]
ang = numpy_array_list[10][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Four points
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.35], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.35,length_left+0.2], [y_end,Placement_left_4], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.2,length_left+0.04], [Placement_left_4,Placement_left_4], ls='dashed',color="grey", linewidth=2,clip_on = False)


# Membrane traffic #

length = 0.2
x_start = numpy_array_list[11][1]
y_start = numpy_array_list[11][2]
ang = numpy_array_list[11][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Four points
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.4], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.4,length_left+0.2], [y_end,Placement_left_5], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.2,length_left+0.04], [Placement_left_5,Placement_left_5], ls='dashed',color="grey", linewidth=2,clip_on = False)

# Extracellular matrix #

length = 0.20
x_start = numpy_array_list[12][1]
y_start = numpy_array_list[12][2]
ang = numpy_array_list[12][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Four points
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.45], [y_end,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.45,length_left+0.2], [y_end,Placement_left_6], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([length_left+0.2,length_left+0.04], [Placement_left_6,Placement_left_6], ls='dashed',color="grey", linewidth=2,clip_on = False)

# RNA metabolism #

length = 0.15
x_start = numpy_array_list[13][1]
y_start = numpy_array_list[13][2]
ang = numpy_array_list[13][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)

plt.plot([x_end,length_left+0.2], [y_end,Placement_left_7], ls='dashed',color="grey", linewidth=2,clip_on = False)

plt.plot([length_left+0.2,length_left+0.04], [Placement_left_7,Placement_left_7], ls='dashed',color="grey", linewidth=2,clip_on = False)

# Small group collection #

small_group_length = 0.10
small_group_collection_point = 0.20


length = small_group_length
x_start = numpy_array_list[14][1]
y_start = numpy_array_list[14][2]
ang = numpy_array_list[14][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)

length = small_group_length
x_start = numpy_array_list[15][1]
y_start = numpy_array_list[15][2]
ang = numpy_array_list[15][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[16][1]
y_start = numpy_array_list[16][2]
ang = numpy_array_list[16][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[17][1]
y_start = numpy_array_list[17][2]
ang = numpy_array_list[17][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[18][1]
y_start = numpy_array_list[18][2]
ang = numpy_array_list[18][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[19][1]
y_start = numpy_array_list[19][2]
ang = numpy_array_list[19][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[20][1]
y_start = numpy_array_list[20][2]
ang = numpy_array_list[20][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)
# Next #
length = small_group_length
x_start = numpy_array_list[21][1]
y_start = numpy_array_list[21][2]
ang = numpy_array_list[21][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+small_group_collection_point], [y_end,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)

plt.plot([length_left+small_group_collection_point,length_left+0.025], [Placement_left_8,Placement_left_8], ls='dashed',color="grey", linewidth=2,clip_on = False)

## unclassified
length = 0.01
x_start = numpy_array_list[22][1]
y_start = numpy_array_list[22][2]
ang = numpy_array_list[22][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=2,clip_on = False)
plt.plot([x_end,length_left+0.04], [Placement_left_9,Placement_left_9], ls='dashed',color="grey", linewidth=2,clip_on = False)

if save_file == 'yes':
    print("Saving plot...")
    plt.savefig(os.path.join(Folder2,File3),dpi=600,bbox_inches='tight')
    print("Done.")
else:
    plt.show()