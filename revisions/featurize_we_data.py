import os
import sys
import numpy as np

#custom PC calculation method
from protein_ligand_water_contacts import main

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

we_data_paths = pathdict[abspath.split("/")[-1]]
we_data_path_1 = we_data_paths[0]
we_data_path_2 = we_data_paths[1]

#####################################################################################################
#                                            FILE PATHS
#####################################################################################################

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


    for xwalker in range(9999):

        #if xwalker % 10 == 0:
        #print(f"walker {xwalker}")
        #print(os.listdir())

        wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
        if os.path.exists(wwfolder):
            try:
                frame = md.load(f"{wwfolder}/traj_comp.xtc", top = f"{abspath}/topology/input.gro")[-1]
                #pc = calc_cftr_pc(frame)
                
                if check_existing:
                    print(pc)
                    extantdata = np.load(f"pc_data_{xround}_v1.npy")
                    for ed in extantdata:
                        if int(ed[0]) == xround and int(ed[1]) == xwalker:
                            print(ed[2])
                    sys.exit(0)


                pcs_all.append([xround, xwalker, pc])            
            except Exception as e:
                #in practice this is to deal with the tiny fraction of files which are corrupted
                print(e)
                pcs_all.append([xround, xwalker, -1])            

            #wire_inds.append(inds)
        else:
            break

    #pcs_all.append(pcs)
    #wire_inds_all.append(wire_inds)

    #remove unpacked archive; if statement in case something happened to the archive while the folder was being processed
    if os.path.exists(archive):
        os.system(f"rm -r {abspath}/traj_segs/{str(xround).zfill(6)}/")

    # if xround % 10 == 0:
    np.save(f"{abspath}/pc_data_{xround}_v{serial}", pcs_all)

        # with open(f"{abspath}/wire_inds_{xround}_v{serial}", "w") as f:
        #     wr = csv.writer(f)
        #     wr.writerows(wire_inds_all)

    # else:
    #     break

#print(np.array(pcs_all))
#print(pcs_all)
#np.save(f"{abspath}/pc_data_{sys.argv[1]}_{sys.argv[2]}_v{serial}", pcs_all)

# with open(f"{abspath}/wire_inds_{sys.argv[1]}_{sys.argv[2]}_v{serial}", "w") as f:
#     wr = csv.writer(f)
#     wr.writerows(wire_inds_all)
        
