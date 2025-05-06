import numpy as np
import mdtraj as md
import sys



#----------methods----------

def mean_z_coord(trj, j, atominds, edge_pad=0.231):
    
    boxz = trj.unitcell_lengths[j][2]
    
    zcoords = [trj.xyz[j][i][2] for i in atominds]

    #if the membrane lies on the box boundary, shift the bottom half of the box up past the top half
    #so that the membrane should no longer cross the boundary
    if max(zcoords) > boxz - edge_pad and min(zcoords) < edge_pad:
        zcoords = [i + boxz if i < boxz/2 else i for i in zcoords]
        
    #calculate mean z coordinate
    mean_z = np.average(zcoords, weights = [trj.top.atom(i).element.mass for i in atominds])

    #make sure mean z coordinate lies in the box
    #  the z coordinate can only lie outside of the box if the case above was activated and 
    #  the center of mass happened to lie near the (now-shifted) bottom half of the box
    #  which got moved to be just above the top
    #  so there is no need to check for mean_z < 0 as this will never occur
    if mean_z > boxz:
        mean_z -= boxz
        
    return mean_z


#check if going around the box the other way would get you to the membrane faster
#adjust z coordinates to represent the shortest membrane-ligand distance
#retaining sign information
def periodic_distance_correction(boxz, lrz):

    #this returns the magnitude of the shortest z distance but not its sign
    #ligand_rel_zcoord_centered = min([abs(lrz), abs(lrz - boxz), abs(lrz + boxz)])

    #this preserves the sign of the delta-z coordinate
    if lrz > boxz/2:
        ligand_rel_zcoord_centered = lrz - boxz
    elif lrz < -boxz/2:
        ligand_rel_zcoord_centered = lrz + boxz
    else:
        ligand_rel_zcoord_centered = lrz
        
    return ligand_rel_zcoord_centered


def ljp_ring_query(trj):

    n_ligand_atoms = len(trj.top.select("resname LJP"))
    if n_ligand_atoms == 71:
        lipidated = True
    elif n_ligand_atoms == 44:
        lipidated = False
    else:
        print(f"error: {n_ligand_atoms} atoms in ligand")
        sys.exit(f"error: {n_ligand_atoms} atoms in ligand")


    bicyclic_ring_atom_names_glpg1837 =  {  'S':[0],\
                                            'C':[1, 2, 3, 4, 5, 6, 7, 8, 14,13,15,16],\
                                            'O':[1,3],\
                                            'N':[4,1],\
                                            'H':[1,2,3, 4,5,6, 9,10,11, 15,16,17, 7,8, 18,19, 20]}

    bicyclic_ring_atom_names_lipidated = {  'S':[0],\
                                            'C':[i for i in range(1,13)],\
                                            'O':[1,2],\
                                            'N':[1,2],\
                                            'H':[i for i in range(23,41)]} 

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

    return query


#compute the z coordinate of the ligand center of mass relative to the membrane center of mass, 
#corrected for periodic box boundary crossing.
#inputs: an mdtraj trajectory, a frame index, and True or False to indicate whether the ligand is lipidated
#outputs: the ligand z coordinate relative to the membrane at frame j of trj
def __calc_pc(trj, j):

    #----------select components----------

    #select lipids
    lipidnames=["POPC"]

    lipidresns = ["resn "+i for i in lipidnames]
    lipidquery = " or ".join(lipidresns)
    lipid = trj.top.select(lipidquery)

    #select ligand
    ligandname="LJP"
    
    atomname_query = ljp_ring_query(trj)

    ligand = trj.top.select("resn %s" % ligandname)
    ligand_head = trj.top.select(f"resn %s and (%s)" % (ligandname, atomname_query))
    ligand_r = [i for i in ligand if i not in ligand_head]
    #ligand_r = trj.top.select(f"resn {ligandname} and not ({atomname_query})")

    #how close does an object need to be to both sides of the box to be assumed to cross the boundary
    #this value was selected to be 1.5 times the 1.54 A C-C single bond length, 
    #which is about the longest that normal C-C bonds get
    edge_pad = 0.231 #nm
    
    #----------calculate centers of mass for each group

    mlpz = mean_z_coord(trj, j, lipid, edge_pad=edge_pad)
    
    mlgz = mean_z_coord(trj, j, ligand, edge_pad=edge_pad)
    mlgz_head = mean_z_coord(trj, j, ligand_head, edge_pad=edge_pad)
    mlgz_r = mean_z_coord(trj, j, ligand_r, edge_pad=edge_pad)

    #----------ligand-membrane distance
    
    #ligand z coordinate relative to membrane center
    lrz = mlgz - mlpz
    lrz_head = mlgz_head - mlpz
    lrz_r = mlgz_r - mlpz

    #----------periodic distance correction

    boxz = trj.unitcell_lengths[j][2]

    lrz = periodic_distance_correction(boxz, lrz)
    lrz_head = periodic_distance_correction(boxz, lrz_head)
    lrz_r = periodic_distance_correction(boxz, lrz_r)

    return [lrz, lrz_head, lrz_r]

#----------calculate and return progress coordinates----------

def calc_pcoord(trj, top):

    trj = md.load(trj, top=top)

    return __calc_pc(trj, -1)
