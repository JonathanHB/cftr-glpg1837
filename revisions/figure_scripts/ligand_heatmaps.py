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


def save_png(inputfile, outputfile, run):
    cmd.delete("all")
    cmd.load(inputfile, "a")
    cmd.hide("everything")
    cmd.create("ligand", "resn LJP and not elem H")
    cmd.show("sticks", "object ligand")

    #fix bond orders
    if run == "abbv-974-1" or run == "abbv-974-2":
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
        if run == "abbv-974-2":
            cmd.rotate("x",180)


    elif run == "cftri-c10-1" or run == "cftri-c10-2":

        cmd.bond("name C23","name C21",4)
        cmd.bond("name C23","name C15",4)
        cmd.bond("name C16","name C15",4)
        cmd.bond("name C16","name S1",4)
        cmd.bond("name C21","name S1",4)

        cmd.unbond("name C24","name O4")
        cmd.bond("name C24","name O4",2)
        cmd.unbond("name C20","name O3")
        cmd.bond("name C20","name O3",2)

        #orient has four local minima depending on the starting orientation
        # cmd.set_view((\
        #     -0.540823638,    0.075270735,    0.837761343,\
        #     0.709082484,    0.576545119,    0.405952841,\
        #     -0.452450871,    0.813590825,   -0.365182400,\
        #     0.000000000,    0.000000000,  -57.965442657,\
        #     68.898712158,   46.204513550,  114.862274170,\
        #     45.700397491,   70.230491638,   20.000000000 ))

        # cmd.set_view((\
        #     -0.440988332,    0.322001606,    0.837761581,\
        #     0.897402883,    0.172834381,    0.405952185,\
        #     -0.014076455,    0.930829763,   -0.365182847,\
        #     0.000051909,    0.000119299,  -49.102256775,\
        #     70.075874329,   43.904624939,  115.004043579,\
        #     37.985584259,   60.220439911,   20.000000000 ))
        cmd.set_view((\
            -0.440988332,    0.322001606,    0.837761581,\
            0.897402883,    0.172834381,    0.405952185,\
            -0.014076455,    0.930829763,   -0.365182847,\
            0.000051640,    0.000110404,  -49.102142334,\
            70.019950867,   43.736694336,  114.688774109,\
            37.985584259,   60.220439911,   20.000000000 ))
        cmd.rotate("x", 180)

    cmd.spectrum("b", "white_blue")
    cmd.ray()
    cmd.png(outputfile, width=1500, height=1000, dpi=600)


################################################################################
#-------------------loop over different WE bins and molecules------------------#
################################################################################

run_num = 3
runs = ["abbv-974-1",
        "abbv-974-2",
        "cftri-c10-1",
        "cftri-c10-2"]
run = runs[run_num]

pdb_dir = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-visualization"
projection_structure = "eq"
contact_types = ["lig_prot", "lig_lip", "lig_wat"]
bins = [1,10,20,30,40]

for contact_type in contact_types:
    for bin in bins:
        inputfile = f"{pdb_dir}/input_{contact_type}_{bin}.pdb"
        outputfolder = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        outputfile=f"{outputfolder}/{projection_structure}_{contact_type}_{bin}.png"

        save_png(inputfile, outputfile, run)
        #break
    #break

