import scipy.stats
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances, align
import MDAnalysis.transformations as trans
import sys
import os
from ref_coordinates import pre_min, pre_min_ligand, pre_min_protein

#----------load trajectory----------
currdir = os.getcwd()
# load in end of eq19 for testing
#curr_trj = mda.Universe('000080/input.gro', '000080/seg.trr')
#curr_trj = mda.Universe(currdir+'/input.gro', currdir+'/seg.trr')
curr_trj = mda.Universe('wstp_lip_glpg_1/topology/input.gro', 'wstp_lip_glpg_1/001980-000132-ancestors/traj_segs/001980/000132/traj_comp.xtc')
#wstp_lip_glpg_1/001980-000132-ancestors/traj_segs/001980/000132/traj_comp.xtc
#curr_trj = mda.Universe('wstp_cftr_1_degrabo/topology/input.gro', 'wstp_cftr_1_degrabo/001456-000152-ancestors/001456/000152/traj_comp.xtc')
curr_protein = curr_trj.select_atoms('protein')
not_protein = curr_trj.select_atoms('not protein')

#remove this for now
# transforms = [trans.unwrap(curr_protein),
#               trans.center_in_box(curr_protein, wrap=True),
#               trans.wrap(not_protein)]

# curr_trj.trajectory.add_transformations(*transforms)
#----------calculate positions and distances----------
# take the last frame of the trajectory and get the distance from the position of the GLPG
# to the starting position of the GLPG in the CFTR channel
# fit current protein to reference protein position


# grab the last frame of the trajectory
curr_trj.trajectory[-1]
align.alignto(curr_protein, pre_min_protein)

# calculate COM of current GLPG position
curr_ligand = curr_trj.select_atoms('resname LJP').center_of_mass('group')
dist = distances.distance_array(curr_ligand, pre_min_ligand, box = pre_min.dimensions)

print(*dist[0])
