import os
import sys
import numpy as np
import mdtraj as md

#custom method
import calc_cftr_pc from calc_cftr_pc

abspath = ""
lipidated = False

serial = 1 #versioning

pcs_all = []
wire_inds_all = []

#for all westpa rounds
for xround in range(int(sys.argv[1]), int(sys.argv[2])):

    print(f"round {xround}")

    #name of the tar file for this westpa round
    archive = f"{abspath}/traj_segs/round-{str(xround).zfill(6)}-segs.tar.gz"

    if os.path.exists(archive):
        os.system(f"tar zxf {archive}")

        #for all walkers in westpa round
        #pcs = []
        wire_inds = []

        for xwalker in range(999):

            if xwalker % 10 == 0:
                print(f"walker {xwalker}")

            wwfolder = f"{abspath}/traj_segs/{str(xround).zfill(6)}/{str(xwalker).zfill(6)}"
            if os.path.exists(wwfolder):
                frame = md.load(f"{wwfolder}/traj_comp.xtc", top = f"{abspath}/gromacs_config/input.gro")[-1]
                pc = calc_cftr_pc(frame)
                pcs_all.append([xround, xwalker, pc])
                wire_inds.append(inds)
            else:
                break

        #pcs_all.append(pcs)
        wire_inds_all.append(wire_inds)

        #remove unpacked archive; redundant if statement in case something happened to the archive while the folder was being processed
        if os.path.exists(archive):
            os.system(f"rm -r {abspath}/traj_segs/{str(xround).zfill(6)}/")

        if xround % 10 == 0:
            np.save(f"{abspath}/pc_data_{xround}_v{serial}", pcs_all)

            with open(f"{abspath}/wire_inds_{xround}_v{serial}", "w") as f:
                wr = csv.writer(f)
                wr.writerows(wire_inds_all)

    else:
        break

#print(np.array(pcs_all))
#print(pcs_all)
np.save(f"{abspath}/pc_data_{sys.argv[1]}_{sys.argv[2]}_v{serial}", pcs_all)

with open(f"{abspath}/wire_inds_{sys.argv[1]}_{sys.argv[2]}_v{serial}", "w") as f:
    wr = csv.writer(f)
    wr.writerows(wire_inds_all)
        