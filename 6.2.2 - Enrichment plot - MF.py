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

## Global variables ##

save_file = "yes"
#save_file = "no"

font_size_text = 42
title_font = 42
fig_x = 10
fig_y = 10

placement_correction_unit = 0.001

#font_size = 42
length_middle = 0.15

# Length of left #
length_left = -1.7

line_width = 3.0


## Folders ##

Folder1 = "Data/Panther/LRT diff/Enrichment data"
Folder2 = "Results/Enrichment/Molecular function"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "MF_enrichment_data_all.xlsx"

File2 = "Enrichment_all.png"


## Load data ##

df_all = pd.read_excel(os.path.join(Folder1,File1))
#df_classified = pd.read_excel(os.path.join(Folder1,File2))

## Variables ##

data_all = df_all['Count']
labels_all = df_all['Class']
percent_list = df_all['Percentage'].tolist()

others_percentage = 0

for key in percent_list:
    #key = key.replace("%","")
    key_value = float(key)
    if key_value < 1.0:
        others_percentage += key_value
others_percentage = round(others_percentage,1)

## Set colors ##

pal = sns.color_palette("crest_r",n_colors=12)

## Enrichment plot - classified ##
print("Creating enrichment plot for all metabolites")
plt.figure(figsize=(fig_x,fig_y))
wedges, texts = plt.pie(data_all,
        radius=1.0,
        colors  = pal[0:],
        wedgeprops = {"edgecolor":"black",'linewidth': 1,"alpha": 1},
        counterclock=False,
        startangle=45)

#print("Saving plot...")

### Add a hole in the pie  
hole = plt.Circle((0, 0), 0.6, facecolor='white',edgecolor="black", linewidth=1)
plt.gcf().gca().add_artist(hole)

## Add labels ##

# Molecular function
plt.text(0,0,"Molecular\nfunction",fontsize=title_font,ha='center',va='center',fontweight="bold")

# Placement positions #

Placement_right_1 = float(0.9)
Placement_right_2 = float(-1.4)


line_splace = np.linspace(-1.4, 0.9, num=8)

Placement_left_1 = float(line_splace[0])
Placement_left_2 = float(line_splace[1])
Placement_left_3 = float(line_splace[2])
Placement_left_4 = float(line_splace[3])
Placement_left_5 = float(line_splace[4])
Placement_left_6 = float(line_splace[5])
Placement_left_7 = float(line_splace[6])
Placement_left_8 = float(line_splace[7])


# binding #
plt.text(1.5,Placement_right_1,"{} ({}%)".format(labels_all[0],percent_list[0]),fontsize=font_size_text,ha='center',va='center')
# catalytic activity #
plt.text(1.1,Placement_right_2,"{} ({}%)".format(labels_all[1],percent_list[1]),fontsize=font_size_text,ha='center',va='center')
# Molecular function regulator #
plt.text(length_left,Placement_left_1,"{} ({}%)".format(labels_all[2],percent_list[2]),fontsize=font_size_text,ha='right',va='center')
# Transporter activitiy #
plt.text(length_left,Placement_left_2,"{} ({}%)".format(labels_all[3],percent_list[3]),fontsize=font_size_text,ha='right',va='center')
# Molecular transducer #
plt.text(length_left,Placement_left_3,"{} ({}%)".format(labels_all[4],percent_list[4]),fontsize=font_size_text,ha='right',va='center')
# Transcription regulator #
plt.text(length_left,Placement_left_4,"{} ({}%)".format(labels_all[5],percent_list[5]),fontsize=font_size_text,ha='right',va='center')
# Structural molecule #
plt.text(length_left,Placement_left_5,"{} ({}%)".format(labels_all[6],percent_list[6]),fontsize=font_size_text,ha='right',va='center')
# ATP-dependent #
plt.text(length_left,Placement_left_6,"{} ({}%)".format(labels_all[7],percent_list[7]),fontsize=font_size_text,ha='right',va='center')

# Small group collection #
plt.text(length_left,Placement_left_7,"Small group collection ({}%)".format(others_percentage),fontsize=font_size_text,ha='right',va='center')

# Unclassified #
plt.text(length_left+0.05,Placement_left_8,"Unclassified ({}%)".format(percent_list[-1]),fontsize=font_size_text,ha='right',va='center')


numpy_array_list = []
for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    numpy_array_list.append([ang,x,y])
      
#plt.plot(numpy_array_list[0][1], numpy_array_list[0][2], 'go--', linewidth=line_width, markersize=12)
#length = 1


## draw lines to labels ##

# binding #
length = 0.375
x_start = numpy_array_list[0][1]
y_start = numpy_array_list[0][2]
ang = numpy_array_list[0][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,x_end], [y_end,Placement_right_1-0.2], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

# catalytic activity #

length = 0.15
x_start = numpy_array_list[1][1]
y_start = numpy_array_list[1][2]
ang = numpy_array_list[1][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,x_end-(Placement_right_2/2)-0.1], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end-(Placement_right_2/2)-0.1,x_end-(Placement_right_2/2)-0.1], [y_end,y_end-0.22], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

### Break left ###

# Transporter activitiy #

length = 0.24
x_start = numpy_array_list[2][1]
y_start = numpy_array_list[2][2]
ang = numpy_array_list[2][0]
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

plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+0.025], [Placement_left_1,Placement_left_1], ls='dashed',color="grey", linewidth=line_width,clip_on = False)



# Molecular function regulator #

length = 0.15
x_start = numpy_array_list[3][1]
y_start = numpy_array_list[3][2]
ang = numpy_array_list[3][0]
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

plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+0.025], [Placement_left_2,Placement_left_2], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

# Molecular transducer #

length = 0.1
x_start = numpy_array_list[4][1]
y_start = numpy_array_list[4][2]
ang = numpy_array_list[4][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# four point #
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+0.1], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+0.1,length_left+length_middle], [y_end,Placement_left_3], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_3,Placement_left_3], ls='dashed',color="grey", linewidth=line_width,clip_on = False)


# Transcription regulator #

length = 0.1
x_start = numpy_array_list[5][1]
y_start = numpy_array_list[5][2]
ang = numpy_array_list[5][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# four point #
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+0.15], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+0.15,length_left+length_middle], [y_end,Placement_left_4], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_4,Placement_left_4], ls='dashed',color="grey", linewidth=line_width,clip_on = False)


# Structural molecule #

length = 0.1
x_start = numpy_array_list[6][1]
y_start = numpy_array_list[6][2]
ang = numpy_array_list[6][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# four point #
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+0.2], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+0.2,length_left+length_middle], [y_end,Placement_left_5], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_5,Placement_left_5], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

# ATP-dependent #

length = 0.1
x_start = numpy_array_list[7][1]
y_start = numpy_array_list[7][2]
ang = numpy_array_list[7][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# four point #
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+0.25], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+0.25,length_left+length_middle], [y_end,Placement_left_6], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_6,Placement_left_6], ls='dashed',color="grey", linewidth=line_width,clip_on = False)


# Small group collection #

small_group_length = 0.10
small_group_collection_point = 0.20
placement_justed = 0.45

length = small_group_length
x_start = numpy_array_list[8][1]
y_start = numpy_array_list[8][2]
ang = numpy_array_list[8][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+placement_justed-0.1+0.01], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+placement_justed-0.1+0.01,length_left+length_middle], [y_end,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_7,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

length = small_group_length
x_start = numpy_array_list[9][1]
y_start = numpy_array_list[9][2]
ang = numpy_array_list[9][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+placement_justed-0.05], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+placement_justed-0.05,length_left+length_middle], [y_end,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_7,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

length = small_group_length
x_start = numpy_array_list[10][1]
y_start = numpy_array_list[10][2]
ang = numpy_array_list[10][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+length_middle+placement_justed], [y_end,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle+placement_justed,length_left+length_middle], [y_end,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([length_left+length_middle,length_left+0.025], [Placement_left_7,Placement_left_7], ls='dashed',color="grey", linewidth=line_width,clip_on = False)


# Unclassified #

length = 0.04
x_start = numpy_array_list[11][1]
y_start = numpy_array_list[11][2]
ang = numpy_array_list[11][0]
x_end = x_start + length * np.cos(np.deg2rad(ang))
y_end = y_start + length * np.sin(np.deg2rad(ang))
# Adjustment of meeting point of 2 lines with fixed placement point of text and second line #
while_placement = Placement_left_8
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
plt.plot([x_start,x_end], [y_start,y_end], ls='dashed',color="grey", linewidth=line_width,clip_on = False)
plt.plot([x_end,length_left+0.1], [Placement_left_8,Placement_left_8], ls='dashed',color="grey", linewidth=line_width,clip_on = False)

if save_file == 'yes':
    print("Saving plot...")
    plt.savefig(os.path.join(Folder2,File2),dpi=600,bbox_inches='tight')
    print("Done.")
else:
    plt.show()
