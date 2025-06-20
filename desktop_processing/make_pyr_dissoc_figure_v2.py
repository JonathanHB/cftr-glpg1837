import pymol
from pymol import cmd
#the linter is wrong about the imports here

import os

upperpath = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation"


colors = ["b", "c", "g", "y", "o"]
folders = [["nonlip_glpg_1", "001655-000156-ancestors-2.5A-20A"], ["nonlip_glpg_1", "001913-000187-ancestors-2.5A-20A"], ["nonlip_glpg_2", "000691-000198-ancestors-2.5A-20A"]]
frame_data = [["bound", [0,0,0]],
              ["sideways", [245, 246, 69]],
              ["h_bond_tangential", [256, 356, 79]],
              ["h_bond_radial", [338, 462, 128]],
              ["unbound", [-1, 573, 207]]]

trj_ind = 1
folder = folders[trj_ind]
#state_ind = 3


#------------------------------------load data--------------------------------------------
cmd.delete("all")

cmd.load(f"{upperpath}/{folders[0][0]}/topology/input.gro", "ref")


#this loop is orthogonal to the one in the original script
for frames, color in zip(frame_data, colors):
    
    frame = frames[1][trj_ind]
    if frame == -1:
        continue

    objname = f"{folder[0]}-{folder[1]}-{frames[0]}"

    cmd.load(f"{upperpath}/{folder[0]}/topology/input.gro", objname)

    we_round = int(round(frame*10/3))
    subframe = int(round(3*((frame*10/3) % 1)))

    #print([f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if f[-14:] == "-traj_comp.xtc"])
    trj_files = [f for f in os.listdir(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/") if (f.endswith("-traj_comp.xtc") and int(f[0:6]) == we_round+1)]

    cmd.load_traj(f"{upperpath}/{folder[0]}/2.5A-20A/{folder[1]}/{trj_files[0]}", objname, start=subframe+1, stop=subframe+1)

    cmd.align(objname, "ref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")

    #the linter is wrong about 'util' being syntactically incorrect here
    if color == "g":
        util.cbag(objname)
    elif color == "c":
        util.cbac(objname)
    elif color == "b":
        util.cbab(objname)
    elif color == "y":
        util.cbay(objname)
    elif color == "m":
        util.cbam(objname)
    elif color == "o":
        util.cbao(objname)


#------------------------------------graphics--------------------------------------------

#secondary structure assignment
cmd.dss("ref")

#hide/show
cmd.hide("everything")

cmd.show("cart", "ref and poly")

cmd.show("sticks", "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and not name C+N+O")
cmd.show("spheres", "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and name CA")
cmd.set("sphere_scale", 0.3, "ref and resi 229+233+236+304+305+308+309+312+313+316+928+930+931+932 and name CA")

cmd.show("spheres", "ref and resi 873+933 and not name C+N+O") #926, 931, and 932 form interactions in some cases, but not in others

cmd.show("sticks", "resname LJP and not ref")

cmd.hide("sticks", "elem H and not (resname LJP and name H12+H13+H14)")
cmd.hide("spheres", "elem H")

#cmd.show("spheres", "name P31")

# cmd.hide("sticks", "resn PA+PC+OL")
# cmd.hide("sticks", "resn ACE+NME")

# cmd.hide("nb_spheres")
# cmd.hide("spheres")

# cmd.hide("sticks", "resn CLR")
# cmd.hide("sticks", "element H")


#coloring
util.cbaw("poly")

##mark edge of membrane
#cmd.color("black", "resi 77+149 or resi 192+245 or resi 298+362 or resi 988+1034 or resi 857 or resi 942 or resi 1094+1154") 


#set view

cmd.set_view((\
    -0.912021697,   -0.409308732,   -0.025543697,\
     0.393961310,   -0.857123911,   -0.331812501,\
     0.113928899,   -0.312682658,    0.942980587,\
     0.000874704,    0.003530707,  -93.388008118,\
    41.787548065,   41.829425812,  109.479858398,\
    65.212921143,  120.734275818,  -20.000000000 ))

cmd.png(f"pyr_dissoc_{folders[trj_ind][0]}_{folders[trj_ind][1]}.png", width=2400, height=1800, ray=True)
#cmd.center("ref and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")
#cmd.center("resi 305+236+232+233+229+309+308+311+312+313+316+315+304+925+926+927+928+929+930+931+932+933+934+935+936+861+862+863+864+865+866+867+868+869+870+871+872+873+874+875+876+877+996")

# cmd.set_view((\
#     -0.911331832,   -0.411233664,   -0.018152758,\
#      0.394373268,   -0.859640598,   -0.324736625,\
#      0.117947072,   -0.303101003,    0.945611656,\
#      0.000000000,    0.000000000, -107.594596863,\
#     38.358119965,   44.360671997,  110.733192444,\
#     75.159851074,  140.029495239,  -20.000000000 ))