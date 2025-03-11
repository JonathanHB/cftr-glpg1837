import numpy as np
import mdtraj as md

def tmd_query():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    #print("color blue, " + " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis]))

    segment_resis_all = [i for sr in segment_resis for i in range(sr[0], sr[1]+1)]
    query = " or ".join([f"resSeq {sr}" for sr in segment_resis_all])

    #indices = frame.top.select(f"protein and ({query})")

    #print_indices = frame.top.select(f"protein and name CA and ({query})")
    #print("+".join([str(i+1) for i in print_indices]))
    return f"protein and ({query})"


def ring_inds(lipidated):
    #see https://pymolwiki.org/index.php/Label
    #i.e. "label i. 1453, name"

    bicyclic_ring_atom_names_glpg1837 =  {'S':[0],\
                                          'C':[1, 2, 3, 4, 5, 6, 7, 8, 14,13,15,16],\
                                          'O':[1,3],\
                                          'N':[4,1],\
                                          'H':[1,2,3, 4,5,6, 9,10,11, 15,16,17, 7,8, 18,19, 20]}

    bicyclic_ring_atom_names_lipidated = {'S':[1],\
                                          'C':[11,12,13,14,15,16,17,18,22,21,23,24],\
                                          'O':[2,4],\
                                          'N':[2,1],\
                                          'H':[23,24,25, 27,28,29, 32,33,34, 36,37,38, 30,31, 39,40, 26]}


    if not lipidated:
        atomlist = bicyclic_ring_atom_names_glpg1837
    else:
        atomlist = bicyclic_ring_atom_names_lipidated


    atom_query = []
    for k in atomlist.keys():
        for ki in atomlist[k]:
            if ki == 0:
                atom_query.append(k)
            else:
                atom_query.append(k+str(ki))

    query = " or ".join(["name " + qi for qi in atom_query])

    #indices = frame.top.select(f"resname LJP and ({query})")

    #print("+".join([str(i+1) for i in indices]))
    return f"resname LJP and ({query})"


def calc_cftr_pc(frame, lipidated):

    #project coordinates into the xy plane
    frame.xyz[:,:,2] = 0

    #compute center of mass in the xy plane
    prot_xycom = md.compute_center_of_mass(frame, select=tmd_query())[-1]
    lig_xycom  = md.compute_center_of_mass(frame, select=ring_inds(lipidated))[-1]

    #set the first atom to the protein center of mass and the second to the ligand center of mass
    for i in [0,1,2]:
        
        frame.xyz[-1,0,i] = prot_xycom[i]
        frame.xyz[-1,1,i] = lig_xycom[i]

    #compute periodicity-corrected distances between the dummy atoms at indices 0 and 1
    dists = md.compute_distances(frame, [[0,1]], periodic=True, opt=True)

    return dists[0][0]


def test_methods(lipidated):

    if not lipidated:
        frame = md.load("/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/wstp_cftr_1_degrabo/topology/input.gro")[-1]
    else:
        frame = md.load("/media/X01Raid01/Data_Backup/home/jborowsky/cftr-analysis/wstp_lip_glpg_1/topology/input.gro")[-1]

    print(calc_cftr_pc(lipidated))
