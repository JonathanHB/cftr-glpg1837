import pymol
from pymol import cmd

import os

upperpath = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation"


colors = ["g", "c", "m"]
folders = [["nonlip_glpg_1", "001655-000156-ancestors-2.5A-20A"], ["nonlip_glpg_1", "001913-000187-ancestors-2.5A-20A"], ["nonlip_glpg_2", "000691-000198-ancestors-2.5A-20A"]]
frame_data = [["bound", [0,0,0]],
              ["sideways", [245, 246, 69]],
              ["h bond tangential", [256, 356, 79]],
              ["h bond radial", [338, 462, 128]],
              ["unbound", [-1, 573, 207]]]

state_ind = 3

#------------------------------------load data--------------------------------------------

cmd.delete("all")

cmd.load(f"{upperpath}/{folders[0][0]}/topology/input.gro", "ref")

for folder, frame, color in zip(folders, frame_data[state_ind][1], colors):
    if frame == -1:
        continue

    cmd.load(f"{upperpath}/{folder[0]}/topology/input.gro", f"{folder[0]}-{folder[1]}")

    we_round = int(round(frame*10/3))
    subframe = int(round(3*((frame*10/3) % 1)))

    #print([f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if f[-14:] == "-traj_comp.xtc"])
    trj_files = [f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if (f.endswith("-traj_comp.xtc") and int(f[0:6]) == we_round+1)]

    cmd.load_traj(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/{trj_files[0]}", f"{folder[0]}-{folder[1]}", start=subframe+1, stop=subframe+1)

    cmd.align(f"{folder[0]}-{folder[1]}", "ref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")

    #the linter is wrong about util
    if color == "g":
        util.cbag(f"{folder[0]}-{folder[1]}")
    elif color == "c":
        util.cbac(f"{folder[0]}-{folder[1]}")
    elif color == "m":
        util.cbam(f"{folder[0]}-{folder[1]}")
        
#------------------------------------graphics--------------------------------------------
cmd.hide("sticks", "resn PA+PC+OL")

cmd.show("sticks", "resi 873+933 and not elem H")

cmd.center("ref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")
cmd.hide("everything", "ref")

cmd.hide("sticks", "resn ACE+NME")

cmd.hide("nb_spheres")
cmd.hide("spheres")

cmd.hide("sticks", "resn CLR")
cmd.hide("sticks", "element H")