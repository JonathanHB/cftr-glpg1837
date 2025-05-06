#Jonathan Borowsky
#Grabe lab
#9/20/24

#load protein and ligand configurations in pymol, align them, and color the ligand by the dissociation state

import pymol
from pymol import cmd
import os

bin_stride = 4

upperpaths = ["/home/jonathan/Documents/grabelab/cftr/chloe-data/wstp_cftr_1_degrabo",
              "/home/jonathan/grabelab/cftr/chloe-data/wstp_lip_glpg_1"]
upperpath = upperpaths[0]


cmd.delete("all")

cmd.load(f"{upperpath}/topology/input.gro")
cmd.create("input-ca", "input and name CA")

cmd.disable("input")

for bin in sorted(os.listdir(f"{upperpath}/bin_samples")):
    bin_num = int(bin[4:])
    if bin_num%bin_stride == 0:
        for f in os.listdir(f"{upperpath}/bin_samples/{bin}"):
            if f[-3:] == "pdb" and f[-6:-4] == "lo":
                #load alpha carbons and ligand, align both to the reference structure using the alpha carbons, and color them by bin
                cmd.load(f"{upperpath}/bin_samples/{bin}/{f}", f"{bin}-{f}")
                cmd.align(f"{bin}-{f[:-4]}", "object input-ca and poly")
                cmd.color(f"grey{str(bin_num*2).zfill(2)}", f"{bin}-{f[:-4]} and elem C")

cmd.hide("sticks", "elem H or resn ACE+NME")
