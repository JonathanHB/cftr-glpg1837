import os
import sys
import numpy as np
import mdtraj as md

#custom PC calculation method
from calc_cftr_pc import calc_cftr_pc

if len(sys.argv) == 1:
    print("arguments:\n    path to westpa folder for Chloe's simulation, path to westpa folder for my simulation, first round to process, last round to process, 'lip' or 'nonlip' depending upon whether the ligand in the simulation is lipidated\n\
each folder specified as an argument should contain a 'traj_segs' subfolder, and rounds are 1-indexed and inclusive\n\
for example, to process the first 10 rounds of a westpa simulation, enter the following as arguments:\n    path1, path2, 1, 10, if_lip")
    sys.exit(0)

abspath = os.getcwd()

pathdict = {"nonlip_glpg_1": ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo/", "../../jhb-simulation-data/wstp_nonlip_glpg_1_continued/"],
            "nonlip_glpg_2": ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_2_wynton/", "../../jhb-simulation-data/wstp_nonlip_glpg_2_continued/"],
            "lip_glpg_1": ["/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_lip_glpg_1/", "../../jhb-simulation-data/wstp_lip_glpg_1_continued/"],
            "lip_glpg_2": ["../../chloe-cftr-project/wstp_lip_glpg_2/", "../../jhb-simulation-data/wstp_lip_glpg_2_continued/"]}

we_data_paths = pathdict[abspath.split("/")[-1]]
we_data_path_1 = we_data_paths[0]
we_data_path_2 = we_data_paths[1]

# we_data_path_1 = sys.argv[1] #data from Chloe's initial simulation
# we_data_path_2 = sys.argv[2] #data from my extension

# if sys.argv[5] == 'nonlip':
#     lipidated = False
# elif sys.argv[5] == 'lip':
#     lipidated = True
# else:
#     sys.exit(0)

check_existing = False

serial = 1 #versioning

pcs_all = []
#wire_inds_all = []

#for all westpa rounds
for xround in range(int(sys.argv[1]), int(sys.argv[2])+1):

    print(f"round {xround}")

    #name of the tar file for this westpa round
    if os.path.exists(f"{we_data_path_1}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"):
        archive = f"{we_data_path_1}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"
    elif os.path.exists(f"{we_data_path_2}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"):
        archive = f"{we_data_path_2}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"
    else:
        break

    #if os.path.exists(archive):
    os.system(f"tar zxf {archive} -C {abspath}")

    #for all walkers in westpa round
    #pcs = []
    #wire_inds = []

    for xwalker in range(9999):

        #if xwalker % 10 == 0:
        print(f"walker {xwalker}")
        #print(os.listdir())

        wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
        if os.path.exists(wwfolder):
            try:
                frame = md.load(f"{wwfolder}/traj_comp.xtc", top = f"{abspath}/topology/input.gro")[-1]
                pc = calc_cftr_pc(frame)
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

    if xround % 10 == 0:
        np.save(f"{abspath}/pc_data_{xround}_v{serial}", pcs_all)

        # with open(f"{abspath}/wire_inds_{xround}_v{serial}", "w") as f:
        #     wr = csv.writer(f)
        #     wr.writerows(wire_inds_all)

    # else:
    #     break

#print(np.array(pcs_all))
#print(pcs_all)
np.save(f"{abspath}/pc_data_{sys.argv[1]}_{sys.argv[2]}_v{serial}", pcs_all)

# with open(f"{abspath}/wire_inds_{sys.argv[1]}_{sys.argv[2]}_v{serial}", "w") as f:
#     wr = csv.writer(f)
#     wr.writerows(wire_inds_all)
        
