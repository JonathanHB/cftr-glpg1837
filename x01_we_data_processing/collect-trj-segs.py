#import numpy as np
#import matplotlib.pyplot as plt
import h5py
#import sklearn
#from sklearn.preprocessing import normalize
import os
import sys

def walker_ancestors(h5path, walker_round, walker_num):
    
    #load h5 file
    with h5py.File(h5path, 'r') as f:
        
        pcoords = []
        walker_ids = []
        
        #determine names of westpa rounds
        #note that round 1's name is at index 0
        iterations = [iter for iter in f["iterations"]]
        
        for iter_ind in range(walker_round-1, -1, -1):
        
            #this extracts only the iteration name, not the iteration data
            iter_name = iterations[iter_ind] 
            #using the iteration name to extract the data
            
            #log walker ID and progress coordinate
            #zeros are for trimming excess nested array layers
            pcoords.append(f["iterations"][iter_name]["pcoord"][walker_num][0][0])
            walker_ids.append(walker_num)
            
            #update round and walker numbers
            #current_round -= 1
            if iter_ind == 0:
                break
            walker_num = f["iterations"][iter_name]["seg_index"][walker_num][1]

    walker_ids.reverse()
    pcoords.reverse()
    return [walker_ids, pcoords]
            
#---------------------------------paths and target walker information---------------------------------

#specify input file and walker
#abspath = "/wynton/home/grabe/jborowsky/aac1/westpa-18"
#"/Users/jonathanborowsky/Documents/grabelab/aac1-cycle/wynton-trj/we18"
#date_time = "111923-2221" #"102423-0931"

#upperpath = sys.argv[1] #"/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo"
we_data_path_1 = sys.argv[1] #data from Chloe's initial simulation
we_data_path_2 = sys.argv[2] #data from my extension

h5path = f'{we_data_path_2}/west.h5'
#'/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data/wstp_lip_glpg_2_continued/west.h5'
#f'{upperpath}/{runpath}/west.h5'

walker_round = int(sys.argv[3])
walker_num = int(sys.argv[4])

#-----------------------------------------------------------------------------------------------------

#make directory to store relevant traj segs
trj_seg_dir = f"{str(walker_round).zfill(6)}-{str(walker_num).zfill(6)}-ancestors"
if not os.path.exists(trj_seg_dir):
    os.mkdir(trj_seg_dir)

#get ancestors of target walkers
walker_ancestry = walker_ancestors(h5path, walker_round, walker_num)
walker_ids = walker_ancestry[0]

print(walker_ancestry)

#extract trajectory files of the ancestors of the target trajectory
for round in range(1, walker_round+1):
    archive = f"round-{str(round).zfill(6)}-segs"
    #print(archive)
    #print(walker_ids[round])

    #wedp = we_data_path_1

    use_path_2 = False
    if os.path.exists(f"{we_data_path_1}/traj_segs/{archive}.tar.gz"):
        os.system(f"tar -zxf {we_data_path_1}/traj_segs/{archive}.tar.gz traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc")
        if not os.path.exists(f"{we_data_path_1}/traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc"):
            use_path_2 = True
    else:
        use_path_2 = True

    if use_path_2:
        os.system(f"tar -zxf {we_data_path_2}/traj_segs/{archive}.tar.gz traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc")
    #    wedp = we_data_path_2

    #the untarred file comes out named 'traj_segs' and I can't figure out how to untar only one file directly from a .tar.gz to a nonstandard new location
    os.system(f"mv traj_segs/{str(round).zfill(6)}/{str(walker_ids[round-1]).zfill(6)}/traj_comp.xtc {trj_seg_dir}/{str(round).zfill(6)}-{str(walker_ids[round-1]).zfill(6)}-traj_comp.xtc") 





    #mv traj_segs/{str(round+1).zfill(6)}/{str(walker_ids[round]).zfill(6)}/traj_comp.xtc {trj_seg_dir}/{str(round+1).zfill(6)}-{str(walker_ids[round]).zfill(6)}-trj.xtc")
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
