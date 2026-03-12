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

run_num = 3
runs = ["abbv-974-1",
        "abbv-974-2",
        "cftri-c10-1",
        "cftri-c10-2"]
run = runs[run_num]

water_figure_paths = {
    "abbv-974-1": ['/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_1-001913-000187-ancestors-2.5A-20A-round-473-bin1',   '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_1-001913-000187-ancestors-2.5A-20A-round-1286-bin10',    '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_1-001913-000187-ancestors-2.5A-20A-round-1678-bin20',    '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_1-001913-000187-ancestors-2.5A-20A-round-1788-bin30',    '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_1-001913-000187-ancestors-2.5A-20A-round-1856-bin40'],
    "abbv-974-2": ['/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_2-000691-000198-ancestors-2.5A-20A-round-82-bin1',    '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_2-000691-000198-ancestors-2.5A-20A-round-269-bin10',     '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_2-000691-000198-ancestors-2.5A-20A-round-583-bin20',     '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_2-000691-000198-ancestors-2.5A-20A-round-647-bin30',     '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/nonlip_glpg_2-000691-000198-ancestors-2.5A-20A-round-672-bin40'],
    "cftri-c10-1":['/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_1-001998-000130-ancestors-2.5A-20A-round-539-bin1',      '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_1-001998-000130-ancestors-2.5A-20A-round-1569-bin10',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_1-001998-000130-ancestors-2.5A-20A-round-1819-bin20',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_1-001998-000130-ancestors-2.5A-20A-round-1937-bin30',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_1-001998-000130-ancestors-2.5A-20A-round-1985-bin40'],
    "cftri-c10-2":['/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_2-001986-000211-ancestors-2.5A-20A-round-78-bin1',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_2-001986-000211-ancestors-2.5A-20A-round-1379-bin10',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_2-001986-000211-ancestors-2.5A-20A-round-1756-bin20',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_2-001986-000211-ancestors-2.5A-20A-round-1909-bin30',       '/home/jonathan/Documents/grabelab/cftr/cftr-potentiator-SAR/water-membrane-permeation/figures/lip_glpg_2-001986-000211-ancestors-2.5A-20A-round-1957-bin40']
}

serial_in = 1
inpath = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
projection_structure = "eq"
contact_types = ["prot_lig", "prot_lip", "prot_wat", "waters"]
bins = [1,10,20,30,40]

files = [[f"{inpath}/{projection_structure}_{ct}_{bin}_v{serial_in}" for ct in contact_types] for bin in bins]

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

            if i == 3:
                imgpathbase = water_figure_paths[run][j]
            else:
                imgpathbase = files[j][i]

            if k == -1:
                imgpath = f"{imgpathbase}.png"
            else:
                imgpath = f"{imgpathbase}_outline_{k}.png"
            
            if os.path.exists(imgpath):
                #print(i,j)
                im = Image.open(imgpath)
                im.thumbnail((output_panel_width, output_panel_height))
                plot_array_overlay.paste(im, (i*output_panel_width, int(j*output_panel_height*vertical_spacing - output_panel_height*(1-vertical_spacing))))

                I1 = ImageDraw.Draw(plot_array_overlay)
                #label panels
                I1.text((i*output_panel_width, j*output_panel_height*vertical_spacing+0.02*output_panel_height), f"{contact_types[i]}; ({(bins[j]-1)/10:.1f}-{(bins[j])/10:.1f}) nm", font=pdbid_font, fill=(0, 0, 0))

    plot_array_base = Image.alpha_composite(plot_array_base, plot_array_overlay)

serial_out = 1
#save the image
outpath = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
plot_array_base.save(f"{outpath}/protein_contacts_{run}_v{serial_out}.png")
