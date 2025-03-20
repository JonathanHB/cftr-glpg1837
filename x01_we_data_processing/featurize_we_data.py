import os
import sys
import numpy as np
import mdtraj as md

#custom PC calculation method
from calc_cftr_pc import calc_cftr_pc

abspath = os.getcwd()
we_data_path_1 = sys.argv[1] #data from Chloe's initial simulation
we_data_path_2 = sys.argv[2] #data from my extension

if sys.argv[5] == 'nonlip':
    lipidated = False
elif sys.argv[5] == 'lip':
    lipidated = True
else:
    sys.exit(0)

serial = 1 #versioning

pcs_all = []
#wire_inds_all = []

#for all westpa rounds
for xround in range(int(sys.argv[3]), int(sys.argv[4])):

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

	    #kluge for corrupted file
        if sys.argv[6] == "nonlip_1_skip" and xround == 494 and xwalker == 46:
            pcs_all.append([494,46,-1])
            continue

        #if xwalker % 10 == 0:
        print(f"walker {xwalker}")
        #print(os.listdir())

        wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
        if os.path.exists(wwfolder):
            frame = md.load(f"{wwfolder}/traj_comp.xtc", top = f"{abspath}/topology/input.gro")[-1]
            pc = calc_cftr_pc(frame, lipidated)
            pcs_all.append([xround, xwalker, pc])
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
np.save(f"{abspath}/pc_data_{sys.argv[3]}_{sys.argv[4]}_v{serial}", pcs_all)

# with open(f"{abspath}/wire_inds_{sys.argv[1]}_{sys.argv[2]}_v{serial}", "w") as f:
#     wr = csv.writer(f)
#     wr.writerows(wire_inds_all)
        
