import os
import sys
import numpy as np

#custom PC calculation method
from protein_ligand_water_contacts import main, get_n_observables

#####################################################################################################
#                                            USER INPUT
#####################################################################################################

init_round = int(sys.argv[1])
final_round = int(sys.argv[2])+1

#####################################################################################################
#                                            FILE PATHS
#####################################################################################################

abspath = os.getcwd()

csheen_folder = "/media/X01Raid01/Data_Backup/home/csheen/cftr-project"
jborowsky_folder = "/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/jhb-simulation-data"

pathdict = {"nonlip_glpg_1": [f"{csheen_folder}/wstp_cftr_1_degrabo/", f"{jborowsky_folder}/wstp_nonlip_glpg_1_continued/"],
            "nonlip_glpg_2": [f"{csheen_folder}/wstp_cftr_2_wynton/", f"{jborowsky_folder}/wstp_nonlip_glpg_2_continued/"],
            "lip_glpg_1": [f"{csheen_folder}/wstp_lip_glpg_1/", f"{jborowsky_folder}/wstp_lip_glpg_1_continued/"],
            "lip_glpg_2": [f"{jborowsky_folder}/chloe-cftr-project/wstp_lip_glpg_2/", f"{jborowsky_folder}/wstp_lip_glpg_2_continued/"]}

pyr_ref_path = "/media/X01Raid01/Data_Backup/shared/CFTR_potentiators/equilibration/ABBV-974-equil/run05/output/eq19.gro"
c10_ref_path = "/media/X01Raid01/Data_Backup/shared/CFTR_potentiators/equilibration/CFTRi-C10-equil/run01/output/eq19.gro"
refpathdict = {"nonlip_glpg_1": pyr_ref_path,
               "nonlip_glpg_2": pyr_ref_path,
               "lip_glpg_1": c10_ref_path,
               "lip_glpg_2": c10_ref_path
}

run_name = abspath.split("/")[-1]

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

#for all westpa rounds
for xround in range(init_round, final_round):

    print(f"round {xround}")

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

    observables_allwalkers = []

    for xwalker in range(9999):

        #if xwalker % 10 == 0:
        #print(f"walker {xwalker}")
        #print(os.listdir())

        wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
        if os.path.exists(wwfolder):
            try:
                observables = main(refpath, f"{abspath}/topology/input.gro", f"{wwfolder}/traj_comp.xtc")
                
                #this is some kind of debugging code
                if check_existing:
                    extantdata = np.load(f"pc_data_{xround}_v1.npy")
                    for ed in extantdata:
                        if int(ed[0]) == xround and int(ed[1]) == xwalker:
                            print(ed[2])
                    sys.exit(0)


                observables_allwalkers.append((xwalker) + observables)

            #handle the tiny fraction of files which are corrupted without crashing
            except Exception as e:
                print(e)
                observables_allwalkers.append((xwalker, [None for a in range(n_observables)]))

        else:
            break

    #remove unpacked archive; if statement is used in case something happened to the archive while the folder was being processed
    if os.path.exists(archive):
        os.system(f"rm -r {abspath}/traj_segs/{str(xround).zfill(6)}/")


    for i_obs in range(n_observables):
        values = [o[i_obs] for o in observables_allwalkers if o is not None]

        #THIS CODE IS NOT GENERAL; IT DEPENDS ON THE DETAILS OF THE IMPORTED main() METHOD    
        if i_obs == 0 or i_obs == 1:
            output = np.concatenate(values)
        else:
            output = np.stack(values)

        np.save(f"{abspath}/pc_data_round_{xround}_obs_{i_obs}_v{serial}", output)
