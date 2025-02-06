#Jonathan Borowsky
#Grabe lab
#091924

#a script to extract random walkers in each progress coordinate interval in order to 
# gain some intuition about what different progress coordinate values mean

import numpy as np
import h5py
import os
#import matplotlib.pyplot as plt
#import matplotlib as mpl

#load walker weights and progress coordinates from west.h5 file

def load_h5_pc_data(h5_path, miniter=0, maxiter=-1, include_last_round=False):

    #---------------------------------------------------
    #load westpa output file

    with h5py.File(h5_path, 'r') as f:

        iter_data = [iter for iter in f["iterations"]]

        if maxiter == -1:
            if include_last_round:
                maxiter = len(iter_data)
            else:
                maxiter = len(iter_data) - 1

        #note that the number of containing arrays is sensitive to 
        #whether the data is extracted by indexing or list comprehension
        pcoord_ndim = len(f["iterations"][f'iter_{str(miniter+1).zfill(8)}']['pcoord'][0][0])

        #this could be done by another layer of list comprehension but this is more readable
        pcoords = []
        for pcd in range(pcoord_ndim):
            pcoords.append([[i[0][pcd] for i in f["iterations"][f'iter_{str(i+1).zfill(8)}']['pcoord']] for i in range(miniter, maxiter)])

        iters = [[i for j in range(len(pcoords[0][i-miniter]))] for i in range(miniter, maxiter)]
        walkers = [[j for j in range(len(pcoords[0][i-miniter]))] for i in range(miniter, maxiter)]
        weights = [[i[0] for i in f["iterations"][f'iter_{str(i+1).zfill(8)}']["seg_index"]] for i in range(miniter, maxiter)]

    print(f"loading westpa iterations {miniter+1} to {maxiter}")

    #---------------------------------------------------
    #make flattened list versions

    pcoords_flat = [[j for i in pcoords[pcd] for j in i ] for pcd in range(pcoord_ndim)]
    iters_flat   = [j for i in iters for j in i ]
    walkers_flat = [j for i in walkers for j in i ]
    #weights_flat = [j for i in weights for j in i ]    
    
    return pcoords_flat, iters_flat, walkers_flat
    
    #weights_flat, maxiter, pcoords, weights


#load data

upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data/wstp_lip_glpg_2_continued"
#"/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_lip_glpg_1"
#"/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/cftr-project/wstp_lip_glpg_2"

h5_path = "/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data/wstp_lip_glpg_2_continued/west.h5"
#upperpath+"/west.h5"
#"/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data/wstp_lip_glpg_2_continued/west.h5"

miniter=1500
maxiter=-1
include_last_round=False

pcoords_flat, iters_flat, walkers_flat = load_h5_pc_data(h5_path, miniter, maxiter, include_last_round)


import random
bins = np.arange(0,50,1)
print(f"bins: {bins}")
data_bins = np.digitize(pcoords_flat, bins)[0]
binned_data = [[i for i,b in enumerate(data_bins) if b == bb+1] for bb in bins]
print(f"number of walkers: {[len(i) for i in binned_data]}")

n_samples = 10
sampled_inds = []
sampled_rounds = []
sampled_walkers = []

for b in binned_data:
    n_walkers_ = min(n_samples, len(b))

    if n_walkers_ != 0:
        sampled_ind = random.choices(b, k=n_walkers_)
        sampled_inds.append(sampled_ind)
        sampled_rounds.append([iters_flat[i] for i in sampled_ind])
        sampled_walkers.append([walkers_flat[i] for i in sampled_ind])

    else:
        sampled_inds.append([])
        sampled_rounds.append([])
        sampled_walkers.append([])


print(sampled_rounds)
print(sampled_walkers)

#TODO draw extra samples and use them in case of corrupted walkers

x = 0
for rr, ww, in zip(sampled_rounds, sampled_walkers):
    
    for r, w, in zip(rr, ww):
        
        archive = f"round-{str(r+1).zfill(6)}-segs"
        if os.path.exists(f"{upperpath}/traj_segs/{archive}.tar.gz"):
            os.system(f"tar -zxf {upperpath}/traj_segs/{archive}.tar.gz traj_segs/{str(r+1).zfill(6)}/{str(w).zfill(6)}/traj_comp.xtc")

    os.system(f"mv traj_segs/ bin-{str(x).zfill(3)}/")

    x+=1

#files triggering errors (rounds 0 indexed)
#948 101
