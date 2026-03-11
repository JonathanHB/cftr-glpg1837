from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import os
from os import listdir
from os.path import isfile, join

import numpy as np
from operator import itemgetter

################################################################################
#--------------make an image panel from a directory full of images-------------#
################################################################################

inpath = "/home/jonathan/Documents/grabelab/cftr/revisions/abbv-974-1-figures"
projection_structure = "eq"
contact_types = ["prot_lig", "prot_lip", "prot_wat"]
bins = [1,10,20,30,40]

files = [[f"{inpath}/{projection_structure}_{ct}_{bin}" for ct in contact_types] for bin in bins]

#panel dimensions in number of images
n_horizontal_panels = len(contact_types) #3 images wide
n_vertical_panels = len(bins) #5 images high

#raw image size in pixels as saved from PyMOL; this value may be slightly different on the new workstation
input_panel_width = 1400
input_panel_height = 1000

#size in pixels for each image in this panel
output_panel_width = 600
output_panel_height = int(round(output_panel_width*input_panel_height/input_panel_width))

#define font for labels
pdbid_font = ImageFont.truetype(f'/home/jonathan/Documents/grabelab/cftr/cftr-glpg1837/revisions/figure_scripts/calibri-regular.ttf', 50)

vertical_spacing = 1.15

plot_array_base = Image.new('RGBA', (n_horizontal_panels*output_panel_width, int(n_vertical_panels*output_panel_height*vertical_spacing)))

for k in range(-1,9):

    #create empty image object with a transparent RGBA background
    plot_array_overlay = Image.new('RGBA', (n_horizontal_panels*output_panel_width, int(n_vertical_panels*output_panel_height*vertical_spacing)))

    #put images into array in arbitrary order
    index = 0
    for i in range(n_horizontal_panels):
        for j in range(n_vertical_panels):
            if k == -1:
                imgpath = f"{files[j][i]}.png"
            else:
                imgpath = f"{files[j][i]}_outline_{k}.png"
            
            if os.path.exists(imgpath):
                #print(i,j)
                im = Image.open(imgpath)
                im.thumbnail((output_panel_width, output_panel_height))
                plot_array_overlay.paste(im, (i*output_panel_width, int(j*output_panel_height*vertical_spacing - output_panel_height*(1-vertical_spacing))))

                I1 = ImageDraw.Draw(plot_array_overlay)
                #label panels
                I1.text((i*output_panel_width, j*output_panel_height*vertical_spacing+0.02*output_panel_height), f"{contact_types[i]}; ({(bins[j]-1)/10:.1f}-{(bins[j])/10:.1f}) nm", font=pdbid_font, fill=(0, 0, 0))

    plot_array_base = Image.alpha_composite(plot_array_base, plot_array_overlay)

serial_out = 2
#save the image
outpath = "/home/jonathan/Documents/grabelab/cftr/revisions/abbv-974-1-figures"
plot_array_base.save(f"{outpath}/protein_contacts_v{serial_out}.png")
