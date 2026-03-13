from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from os import listdir
from os.path import isfile, join

import numpy as np
from operator import itemgetter

################################################################################
#--------------make an image panel from a directory full of images-------------#
################################################################################

#path information
# serial_in = 1
# path = f"/project/bowmanlab/borowsky.jonathan/FAST-cs/gvp-performance/test"

# #list the image files of interest
# path_in = f"{path}/figures_v{serial_in}/"
# files = [f for f in listdir(path_in) if isfile(join(path_in, f)) and f[-4:] == ".png"]

# #PDB IDs, apo where applicable
# prot_names = [f[0:4].lower() for f in files if f[-4:] == ".png"]
# #print([[i,j] for i,j in zip(prot_names, files)])

# #load or define pocket direction and rmsd information to organize figure
# bricks = [['1ofv','A','1ofv','A','None'],['1qys','A','1qys','A','None'],['5bvl','A','5bvl','A','None'],
# ['2fd7','A','2fd7','A','None'],['4hjk','A','4hjk','A','None'],['1igd','A','1igd','A','None'],
# ['4tql','A','4tql','A','None'],['2alp','A','2alp','A','None'],['1amm','A','1amm','A','None']]
# bricks = [b[0] for b in bricks]

# moad_negatives = ['3GZ0','2ZHV','2OSS','5UOJ','1HCL','1BTP','5ZN5','4APE','5X8U','5E4P']

# pocket_dir_info = np.load(f"/project/bowmanlab/borowsky.jonathan/FAST-cs/pocket-tracking/pocket_vols_2/pocket_directions_rmsds_v6.npy", allow_pickle=True)
# distance_deltas = {e[0].lower():e[5][0]-e[5][1] for e in pocket_dir_info}

# #organize panels based on the above information
# files_brick = []
# files_moad = []
# files_pocket = []

# for x, pn in enumerate(prot_names):
#     if pn.lower() in bricks:
#         files_brick.append(files[x])
#     elif pn.upper() in moad_negatives:
#         files_moad.append(files[x])
#     else:
#         files_pocket.append([files[x], distance_deltas[pn.lower()]])

# #sort to put most forward pockets first
# files_pocket = sorted(files_pocket, key=itemgetter(1), reverse=True)
# fns_pocket = [i[0] for i in files_pocket]

# print(files_pocket)

# #move 2fjy next to 1kx9 since the direction calculation for the former is not
# # representative of the true change due to unresolved holo residues and to limitations
# # of the current definition of pocket direction
# kx9_ind = -1
# fjy_ind = -1
# for x, i in enumerate(fns_pocket):
#     if i[0:4].lower() == "1kx9":
#         kx9_ind = x
#     elif i[0:4].lower() == "2fjy":
#         fjy_ind = x

# fjy_obj = fns_pocket.pop(fjy_ind)
# fns_pocket.insert(kx9_ind+1, fjy_obj)

# #assemble indices
# fns_all = fns_pocket + files_brick+files_moad

run_num = 1
runs = ["abbv-974-1",
        "abbv-974-2",
        "cftri-c10-1",
        "cftri-c10-2"]
run = runs[run_num]

inpath = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
projection_structure = "eq"
contact_types = ["lig_prot", "lig_lip", "lig_wat"]
bins = [1,10,20,30,40]

files = [[f"{inpath}/{projection_structure}_{ct}_{bin}.png" for ct in contact_types] for bin in bins]

#panel dimensions in number of images
n_horizontal_panels = len(contact_types) #3 images wide
n_vertical_panels = len(bins) #5 images high

#raw image size in pixels as saved from PyMOL; this value may be slightly different on the new workstation
input_panel_width = 1500
input_panel_height = 1000

#size in pixels for each image in this panel
output_panel_width = 600
output_panel_height = int(round(output_panel_width*input_panel_height/input_panel_width))

# #pads to account for the fact that the screen's aspect ratio differs from that of the average protein
# #these are factors of w_img so that w_img can be used to adjust the image
# #resolution while leaving the ratios of all internal dimensions constant
# w_pad = -int(round(output_panel_width/30))
# h_pad = int(round(output_panel_width/30))

# #overall size of each pane of the panel
# #expressions could be factored but this seems pointless given that the computational expense is miniscule
# w = output_panel_width + w_pad
# h = output_panel_height + h_pad

#define font for labels
pdbid_font = ImageFont.truetype(f'/home/jonathan/Documents/grabelab/cftr/cftr-glpg1837/revisions/figure_scripts/calibri-regular.ttf', 50)

#create empty image object with a transparent RGBA background
plot_array = Image.new('RGBA', (n_horizontal_panels*output_panel_width, n_vertical_panels*output_panel_height))

#put images into array in arbitrary order
#index = 0
for i in range(n_horizontal_panels):
    for j in range(n_vertical_panels):
        #if index < len(files):
        #print(inpath[i][j])
        print(j,i)
        im = Image.open(files[j][i])
        im.thumbnail((output_panel_width, output_panel_height))
        #DO NOT DO THIS IT INVERTS CHIRALITY
        # if run == "abbv-974-2":
        #     im = flipped_img = im.transpose(Image.FLIP_TOP_BOTTOM)

        plot_array.paste(im, (i*output_panel_width, j*output_panel_height))

        I1 = ImageDraw.Draw(plot_array)
        I1.text((i*output_panel_width, j*output_panel_height), f"{contact_types[i]}; ({(bins[j]-1)/10:.1f}-{(bins[j])/10:.1f}) nm", font=pdbid_font, fill=(0, 0, 0))
        #index += 1

serial_out = 1
#save the image
outpath = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
#see if you can export as svg <-- unfortunately this does not work
#fix panel names to be actual words
plot_array.save(f"{outpath}/ligand_contacts_{run}_v{serial_out}.png")
