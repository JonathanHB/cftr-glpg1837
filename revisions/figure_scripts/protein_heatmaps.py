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

#these should probably use cmd.set but maybe there's a reason not to
def set_graphics():
    cmd.do("bg white")
    cmd.do("set ray_shadow=off")
    cmd.do("set cartoon_fancy_helices, 1")
    cmd.do("set specular, off")
    cmd.do("set orthoscopic, 1")
    cmd.do("set depth_cue, 1") #<--enable fog
    cmd.do("set auto_zoom, off")
    cmd.do("set sphere_quality, 5")
    cmd.do("set opaque_background, off")
    cmd.do("set ray_opaque_background, off")
    cmd.do("set antialias, 2")
    cmd.do("set ray_trace_mode, 1")
    cmd.do("set ray_trace_color, black")
    cmd.do("set surface_quality, 1")
    #cmd.do("set ray_trace_color, black")
    cmd.do("set ray_opaque_background, 0")

set_graphics()

def tmd_query_pymol():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    return " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis])


def save_png(inputfile, outputfile, colorscale, reffile):
    
    cmd.delete("all")
    cmd.load(inputfile, "a")

    if run == "cftri-c10-1" or run == "cftri-c10-2":
        cmd.load(reffile, "ref")
        cmd.align("a", f"ref and ({tmd_query_pymol()})")

    cmd.hide("everything")
    cmd.create("pr", "a and poly and not elem H")
    cmd.show("surface", "pr")

    cmd.set_view((\
     0.933495879,    0.014861241,   -0.358267784,\
    -0.357729882,    0.107290752,   -0.927639842,\
     0.024650594,    0.994103849,    0.105471313,\
    -0.000507705,   -0.003106263, -105.927444458,\
    39.996837616,   38.508853912,  114.650337219,\
    84.638229370,  127.218566895,   20.000000000 ))

    # cmd.set_view((\
    #     0.931283236,    0.058830075,   -0.359503239,\
    #     -0.364245355,    0.135353819,   -0.921412706,\
    #     -0.005548787,    0.989036500,    0.147480413,\
    #     -0.000507705,   -0.003106263, -143.839187622,\
    #     39.996837616,   38.508853912,  114.650337219,\
    #     118.203819275,  169.476409912,   20.000000000 ))

    # set_view (\
    #      0.931753218,    0.209694207,   -0.296406806,\
    #     -0.362906545,    0.512254834,   -0.778381467,\
    #     -0.011384944,    0.832823277,    0.553390682,\
    #      0.000102840,   -0.004097723, -141.684692383,\
    #     43.730960846,   44.042785645,  115.130500793,\
    #    121.004135132,  161.892837524,   20.000000000 )

    cmd.spectrum("b", colorscale)
    cmd.set("ray_trace_mode", 1)
    cmd.do("set ray_trace_color, black")
    cmd.set("ray_trace_gain", 0.05)

    cmd.png(outputfile, width=1400, height=1000, dpi=600, ray=1)

    outline_resis = []#["873 and not name OE1"] #[931, 932, "873 and not name OE1", 312, "308 and not name CB and not name CA"]

    for ori, oresi in enumerate(outline_resis):
        
        cmd.set("ray_trace_gain", 0.5)
        cmd.hide("everything")
        cmd.show("surface", f"pr and resi {oresi}")
        cmd.set("ray_trace_mode", 2)
        cmd.do("set ray_trace_color, red")
        cmd.png(f"{outputfile[:-4]}_outline_{ori}.png", width=1200, height=1000, dpi=600, ray=1)



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
contact_types = ["prot_lig", "prot_lip", "prot_wat"]
colorscales = ["white_blue","white_green","white_red"]
bins = [1,10,20,30,40]

serial = 1

for i, contact_type in enumerate(contact_types):
    for bin in bins:
        inputfile = f"{pdb_dir}/input_{contact_type}_{bin}.pdb"
        outputfolder = f"/home/jonathan/Documents/grabelab/cftr/revisions/{run}-figures"
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        outputfile=f"{outputfolder}/{projection_structure}_{contact_type}_{bin}_v{serial}.png"

        reffile = f"/home/jonathan/Documents/grabelab/cftr/revisions/{runs[0]}-visualization/input_{contact_type}_{bin}.pdb"
        save_png(inputfile, outputfile, colorscales[i], reffile)
        #break
    #break

