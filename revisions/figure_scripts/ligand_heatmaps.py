#Jonathan Borowsky
#Grabe lab
#3/9/26

################################################################################
#-----------------------------------imports------------------------------------#
################################################################################

import os

################################################################################
#------------------------------graphics settings-------------------------------#
################################################################################

def set_graphics():
    cmd.do("bg white")
    cmd.do("set ray_shadow=off")
    cmd.do("set cartoon_fancy_helices, 1")
    cmd.do("set specular, off")
    cmd.do("set orthoscopic,1")
    cmd.do("set depth_cue,0") #<--disable fog
    cmd.do("set auto_zoom, off")
    cmd.do("set sphere_quality, 5")
    cmd.do("set opaque_background, off")
    cmd.do("set ray_opaque_background, off")
    cmd.do("set antialias, 2")
    cmd.do("set ray_trace_mode, 1")
    cmd.do("set ray_trace_color, black")

set_graphics()


def save_png(inputfile, outputfile):
    cmd.delete("all")
    cmd.load(inputfile, "a")
    cmd.hide("everything")
    cmd.create("ligand", "resn LJP and not elem H")
    cmd.show("sticks", "object ligand")

    #fix bond orders
    cmd.bond("name C5","name C6",4)
    cmd.bond("name C5","name C15",4)
    cmd.bond("name C13","name C15",4)
    cmd.bond("name C13","name S",4)
    cmd.bond("name C6","name S",4)

    cmd.bond("name C11","name N3",4)
    cmd.bond("name N2","name N3",4)
    cmd.bond("name N2","name C10",4)
    cmd.bond("name C9","name C10",4)
    cmd.bond("name C9","name C11",4)

    cmd.unbond("name C12","name O2")
    cmd.bond("name C12","name O2",2)
    cmd.unbond("name C16","name O3")
    cmd.bond("name C16","name O3",2)


    cmd.orient("object ligand")
    cmd.spectrum("b", "white_blue")
    cmd.ray()
    cmd.png(outputfile, width=1500, height=1000, dpi=600)


################################################################################
#-------------------loop over different WE bins and molecules------------------#
################################################################################

pdb_dir = "/home/jonathan/Documents/grabelab/cftr/revisions/abbv-974-1-visualization"
projection_structure = "eq"
contact_types = ["lig_prot", "lig_lip", "lig_wat"]
bins = [1,10,20,30,40]

for contact_type in contact_types:
    for bin in bins:
        inputfile = f"{pdb_dir}/input_{contact_type}_{bin}.pdb"
        outputfolder = f"/home/jonathan/Documents/grabelab/cftr/revisions/abbv-974-1-figures"
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        outputfile=f"{outputfolder}/{projection_structure}_{contact_type}_{bin}.png"

        save_png(inputfile, outputfile)
        #break
    #break

