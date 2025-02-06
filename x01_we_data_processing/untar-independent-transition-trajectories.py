import numpy as np
import h5py
import os
import sys

#put westpa_msm_functions.py in the same folder I believe
from get_unique_transitions import get_unique_transitions

#-----------------------------------------------methods-----------------------------------------------

#------------------------
# Parameters: 
# pc: a progress coordinate (currently a 1d array [or list] of floats, 
# but could be changed to a 2d array if the PC had multiple time points 
# and the indexing in get_unique_transitions() was adjusted accordingly)
#------------------------
# Returns: 
# macrostate: An integer. The possible values should be consecutive integers counting up from 0, and -1.
#             takes values in the range [-1,4] in this implementation
# note that this method must return a number other than -1 for the macrostate of the initial structure

def pc_2_macrostate(pc):
    pc0_min = 4 #11.5
    pc0_max = 16 #16 
    
    if pc[0] < pc0_min: 
        macrostate = 0
    elif pc[0] > pc0_max:
        macrostate = 1
    else:
        return -1
        
    return macrostate


def get_transition_representatives(h5path, macrostate_classifier, n_macrostates):

    transitions, wbf, itt, pcs_by_tss = get_unique_transitions(h5path, macrostate_classifier, n_macrostates, minround=0, maxround=-1)
    tcounts = np.array([[len(j) for j in i] for i in transitions])
    print(tcounts)
    
    #most_advanced_descendant_of_each_independent_transition
    madoeit = []

    for k in pcs_by_tss.keys():
        #print(k)
        pcs = [i[2] for i in pcs_by_tss[k][1]]

        #print(len(pcs))
        #print(max(pcs))
        #print(min(pcs))

        #print(pcs_by_tss[k][1][np.argmax(pcs)])
        madoeit.append(pcs_by_tss[k][1][np.argmax(pcs)])
        #print()
    #transition_rep_sets = [tr_0_i for tr_0_i in transitions[0][1:] if len(tr_0_i) > 0]

    return [madoeit]

n_macrostates=2


#---------------------------------paths and target walker information---------------------------------

we_data_path_1 = sys.argv[1] #data from Chloe's initial simulation
we_data_path_2 = sys.argv[2] #data from my extension

# upperpath = sys.argv[1]
# runpath = sys.argv[2]

# os.mkdir(runpath)
# os.chdir(runpath)

#os.getcwd() #"/media/X01Raid01/Data_Backup/home/jborowsky/westpa/westpa-28d8"
#print(f"collecting data from: {upperpath}/{runpath}")

#if len(sys.argv) == 1:
h5path = f'{we_data_path_2}/west.h5'
#'/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data/wstp_lip_glpg_2_continued/west.h5'
#f'{upperpath}/{runpath}/west.h5'
# else: 
#     h5path = f'{upperpath}/{sys.argv[1]}.h5'
print(f"reading {h5path}")

#---------------------------copy topology files for trajectory processing-----------------------------

if not os.path.exists("topology"):
    os.mkdir(f"topology")
    os.system(f"cp {we_data_path_1}/bstates/input/input.gro topology")

#os.system(f"cp {upperpath}/{runpath}/bstates/input/index.ndx topology")
#os.system(f"cp {upperpath}/{runpath}/bstates/input/input.gro {runpath}/topology")

#-----------------collect trajectory ancestors for each transition representative---------------------

#make directory to store relevant traj segs
#trj_seg_dir = "representatives-transitions"
# if not os.path.exists(upperpath + "/" + trj_seg_dir):
#     os.mkdir(upperpath + "/" + trj_seg_dir)

#get transition representatives
transition_representatives = get_transition_representatives(h5path, pc_2_macrostate, n_macrostates)
print(transition_representatives)

#get ancestors of target walkers
# walker_ancestry = walker_ancestors(h5path, walker_round, walker_num)
# walker_ids = walker_ancestry[0]


#extract trajectory files of the transition representatives
for ms_ind, tr_set in enumerate(transition_representatives):
    #make a directory for each type of macrostate transition
    #if not os.path.exists(f"{upperpath}/{trj_seg_dir}/reps_{ms_ind}"):
        #os.mkdir(f"{upperpath}/{trj_seg_dir}/reps_{ms_ind}")

    for tr in tr_set:
        print(tr)
        #archive = f"round-{str(tr[0]+1).zfill(6)}-segs"
        #print(archive)

        os.system(f"python3 ../collect-trj-segs.py {we_data_path_1} {we_data_path_2} {tr[0]+1} {tr[1]}")

        #print(walker_ids[round])
        #os.system(f"tar -zxf {upperpath}/traj_segs/{archive}.tar.gz traj_segs/{str(tr[0]+1).zfill(6)}/{str(tr[1]).zfill(6)}/traj_comp.xtc; \
        #        mv {upperpath}/traj_segs/{str(tr[0]+1).zfill(6)}/{str(tr[1]).zfill(6)}/traj_comp.xtc {upperpath}/{trj_seg_dir}/reps_{ms_ind}/{str(tr[0]+1).zfill(6)}-{str(tr[1]).zfill(6)}-trj.xtc")
        #os.system(f"tar -xvf traj_segs/{archive}.tar.gz -C {trj_seg_dir}/")
        #os.system(f"cp {trj_seg_dir}/traj_segs/{str(round+1).zfill(6)}/{str(walker_ids[round]).zfill(6)}/traj_comp.xtc {trj_seg_dir}/{str(round+1).zfill(6)}-trj.xtc;\
        #          rm -r {trj_seg_dir}/traj_segs/{str(round+1).zfill(6)}")
    
#concatenate trajectories
# trjcat_commands = [
#     "module load mpi",
#     "module load Sali",
#     "module load cuda/10.0.130",
#     f"/wynton/home/grabe/shared/gromacs/gromacs-2020.6_CUDA10_SSE4/bin/gmx trjcat -f {trj_seg_dir}/*-trj.xtc -o {trj_seg_dir}/{str(walker_round).zfill(6)}-{str(walker_num).zfill(6)}-full-trj.xtc -cat"
# ]

#os.system("; ".join(trjcat_commands))
