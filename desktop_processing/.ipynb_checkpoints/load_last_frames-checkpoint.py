import os

import pymol
from pymol import cmd
#the linter is wrong about the imports here

upperpath = "/home/jonathan/Documents/grabelab/cftr/last_rounds/und-1/002019"
toppath = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/lip_glpg_1/topology/input.gro"

#upperpath = "/home/jonathan/Documents/grabelab/cftr/last_rounds/pyr-2/001088"
#toppath = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/topology/input.gro"
#upperpath = "/Users/jonathanborowsky/Documents/grabelab/cftr-glpg/last_we_rounds/pyr-1/002010"
#toppath = "/Users/jonathanborowsky/Documents/grabelab/cftr-glpg/independent-partial-dissociation/nonlip_glpg_1/topology"

cmd.delete("all")

cmd.load(toppath)

modphase = 1

for folder in os.listdir(upperpath):
    if int(folder) % 2 == modphase:
        cmd.load(f"{upperpath}/{folder}/seg.gro", f"seg-{folder}")
        cmd.align(f"seg-{folder}", "input and (resi 77-149 or resi 192-245 or resi 298-362 or resi 988-1034 or resi 857-889 or resi 900-942 or resi 1094-1154)")

        cmd.hide("everything")

cmd.show("ribbon", "poly")
cmd.show("sticks", "resn LJP and not elem H")
cmd.show("spheres", "name P31")
cmd.set("sphere_scale", 0.1)
cmd.set("sphere_scale", 0.3, "input")
#cmd.hide("sticks", "not resn LJP")

cmd.color("blue", "index 56904")
cmd.color("blue", "index 57574") #57599 is the index of the last atom in the same molecule
cmd.color("blue", "index 19786") #phosphate index = end index -25

