import scipy.stats
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances, align
import MDAnalysis.transformations as trans
import sys
import numpy as np
import os


#----------select reference components----------
# using the structure right before minimization

pre_min = mda.Universe('/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo/bstates/input/min.gro')
#pre_min = mda.Universe('/wynton/home/grabe/csheen/cftr-project/wstp_cftr_1_degrabo/bstates/input/min.gro')
pre_min_protein = pre_min.select_atoms('protein')
pre_min_not_protein = pre_min.select_atoms('not protein')

ref = mda.Universe('/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo/bstates/input/input.gro')
#ref = mda.Universe('/wynton/home/grabe/csheen/cftr-project/wstp_cftr_1_degrabo/bstates/input/input.gro')
ref_ligand = ref.select_atoms('resname LJP').center_of_mass('group')
#print(ref_ligand)
#print(ref_ligand)

ref_protein = ref.select_atoms('protein')
ref_not_protein = ref.select_atoms('not protein')

"""
transforms = [trans.unwrap(ref_protein),
              trans.center_in_box(ref_protein, wrap=True),
              trans.wrap(ref_not_protein, compound)]

ref.trajectory.add_transformations(*transforms)

transforms2 = [trans.unwrap(step_5_input_protein),
              trans.center_in_box(step_5_input_protein, wrap=True),
              trans.wrap(step_5_input_not_protein)]

step_5_input.trajectory.add_transformations(*transforms2)
#----------calculate positions and distances----------
# take the last frame of the trajectory and get the distance from the position of the GLPG
# to the starting position of the GLPG in the CFTR channel
# fit current protein to reference protein position

"""
align.alignto(pre_min_protein, ref_protein)

# # grab the last frame of the trajectory
pre_min.trajectory[-1]
# calculate COM of current GLPG position
pre_min_ligand = pre_min.select_atoms('resname LJP').center_of_mass('group')
# print(pre_min_ligand)

dist = distances.distance_array(pre_min_ligand, ref_ligand, box = ref.dimensions)

# print(*dist[0])
