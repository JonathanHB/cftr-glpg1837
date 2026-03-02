import pymol
from pymol import cmd
#the linter is wrong about the imports here

upperpath = "/Users/jonathanborowsky/Documents/grabelab/cftr-glpg/last_we_rounds/pyr-1/002010"
toppath = "/Users/jonathanborowsky/Documents/grabelab/cftr-glpg/independent-partial-dissociation/nonlip_glpg_1/topology"
mod = 0

for folder in os.listdir(upperpath):
    if int(folder) % 10 == mod:
        cmd.load(f"{upperpath}/{folder}/seg.gro", f"seg-{folder}")