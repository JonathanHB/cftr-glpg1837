import os
import sys
import numpy as np
import time

import MDAnalysis as mda

#custom PC calculation method
#from protein_ligand_water_contacts import main, get_n_observables
from calc_lig_unaligned_z import main, get_n_observables, observable_names


#####################################################################################################
#                                            USER INPUT
#####################################################################################################

init_round = int(sys.argv[1])
final_round = int(sys.argv[2])+1

#####################################################################################################
#                                            FILE PATHS
#####################################################################################################

abspath = os.getcwd()
run_name_ = abspath.split("/")[-1]

name_mapping = {"abbv-974-1":"nonlip_glpg_1",
               "abbv-974-2":"nonlip_glpg_2",
               "cftri-c10-1":"lip_glpg_1",
               "cftri-c10-2":"lip_glpg_2"
}

run_name = name_mapping[run_name_]

csheen_folder = "/media/X01Raid01/Data_Backup/home/csheen/cftr-project"
jborowsky_folder = "/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data"

#folder containing traj_segs/ directory
pathdict = {"nonlip_glpg_1": [f"{csheen_folder}/wstp_cftr_1_degrabo/", f"{jborowsky_folder}/wstp_nonlip_glpg_1_continued/"],
            "nonlip_glpg_2": [f"{csheen_folder}/wstp_cftr_2_wynton/", f"{jborowsky_folder}/wstp_nonlip_glpg_2_continued/"],
            "lip_glpg_1": [f"{csheen_folder}/wstp_lip_glpg_1/", f"{jborowsky_folder}/wstp_lip_glpg_1_continued/"],
            "lip_glpg_2": [f"/media/X01Raid01/Data_Backup/shared/CFTR_potentiators/westpa/CFTRi-C10-we2/CFTRi-C10-we2-part1", f"{jborowsky_folder}/wstp_lip_glpg_2_continued/"]}

pyr_ref_path = "/media/X01Raid01/Data_Backup/shared/CFTR_potentiators/equilibration/ABBV-974-equil/run05/output/eq19.gro"
c10_ref_path = "/media/X01Raid01/Data_Backup/shared/CFTR_potentiators/equilibration/CFTRi-C10-equil/run01/output/eq19.gro"
refpathdict = {"nonlip_glpg_1": pyr_ref_path,
               "nonlip_glpg_2": pyr_ref_path,
               "lip_glpg_1": c10_ref_path,
               "lip_glpg_2": c10_ref_path
}

refpath = refpathdict[run_name]

we_data_paths = pathdict[run_name]
we_data_path_1 = we_data_paths[0]
we_data_path_2 = we_data_paths[1]

#####################################################################################################
#                                            FILE PATHS
#####################################################################################################

n_observables = get_n_observables()

#for debugging
check_existing = False

#output file versioning in case this needs to be rerun
serial = 1 

increment = 10

#ref = mda.Universe(refpath)

#for all westpa rounds
for xround in range(init_round, final_round, increment):

    print(f"round {xround}")
    t0 = time.time()
    #-----------------------------------------------------------------------
    #figure out where the archive for the current round is and untar it

    tar1 = f"{we_data_path_1}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"
    tar2 = f"{we_data_path_2}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"

    if os.path.exists(tar1):
        archive = tar1
    elif os.path.exists(tar2):
        archive = tar2
    else:
        print(f"archive not found for round {xround}; skipping")
        continue

    os.system(f"tar zxf {archive} -C {abspath}")
    
    t1 = time.time()
    print(f"untarred walkers in {t1-t0} s")
    
    observables_allwalkers = []

    nwalkers = len(os.listdir(f'{abspath}/traj_segs/{str(xround).zfill(6)}'))
    print(f"processing {nwalkers} walkers")

    for xwalker in range(9999):

        #if xwalker % 10 == 0:
        #print(f"walker {xwalker}")
        #print(os.listdir())

        wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
        if os.path.exists(wwfolder):
            try:
                observables = main(refpath, f"{wwfolder}/traj_comp.xtc", run_name_, frame=-1)
                #main(refpath, refpath, f"{wwfolder}/traj_comp.xtc")
                
                #this is some kind of debugging code
                if check_existing:
                    extantdata = np.load(f"pc_data_{xround}_v1.npy")
                    for ed in extantdata:
                        if int(ed[0]) == xround and int(ed[1]) == xwalker:
                            print(ed[2])
                    sys.exit(0)


                observables_allwalkers.append([xwalker] + observables)

            #handle the tiny fraction of files which are corrupted without crashing
            except Exception as e:
                print(f"exception: {e}")
                
                observables_allwalkers.append([None] + [None for a in range(n_observables)])

        else:
            break

    t2 = time.time()
    print(f"processed walkers at {(t2-t1)/nwalkers} seconds per walker")

    #remove unpacked archive; if statement is used in case something happened to the archive while the folder was being processed
    if os.path.exists(archive):
        os.system(f"rm -r {abspath}/traj_segs/{str(xround).zfill(6)}/")


    obs_names = ["walker_numbers"] + observable_names()

    for i_obs in range(n_observables+1):
        #print(i_obs)
        #print(observables_allwalkers[0])
        #print(len(observables_allwalkers))
        #print([o[i_obs] for o in observables_allwalkers])
        values = [o[i_obs] for o in observables_allwalkers if o[i_obs] is not None]

        #THIS CODE IS NOT GENERAL; IT DEPENDS ON THE DETAILS OF THE IMPORTED main() METHOD    
        
        #handle lists of integers
        if i_obs in [0]:#,8,9,11]:
            output = values

        #handle variable number of water coordinates
        #elif i_obs == 1:
            #print(values)   
        #    maxwaters = max([v.shape[0] for v in values])             
            #print(maxwaters)
        #    output = np.zeros((len(values),maxwaters,3))
        #    for i,v in enumerate(values):
        #        output[i,0:len(v)] = v

        #stack fixed-size arrays 
        else:
            output = np.stack(values)


        np.save(f"{abspath}/pc_data_round_{xround}_{obs_names[i_obs]}_v{serial}", output)
