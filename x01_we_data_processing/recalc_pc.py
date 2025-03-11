import numpy as np
import mdtraj as md
import itertools

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse.csgraph import shortest_path

import os
import sys
import h5py

import matplotlib.pyplot as plt
import time
import csv


# Determine the atom indices of all water oxygens in a cylinder 
# with its radial axis lying along z and passing through the protein's center of mass.
#---------------------------------
#Parameters
# frame: an MD trajectory frame
# queries: a list of mdtraj selections for selecting chemical groups through which the water wire may pass
#          these are put in parenthesis and joined by "or" 
#          These groups are all referred to as waters below for conciseness.
# r_dist: the radius in nm of the outer cylinder, the axis of which passes through the protein center of mass
# z_dist_n: the z distance in nm from the protein center of mass to the -z end of the outer cylinder
# z_dist_p: the z distance in nm from the protein center of mass to the +z end of the outer cylinder
# surface_pad: the minimum distance in nm of the inner cylinder from the outer cylinder. 
#              Note that larger pad values mean a smaller inner cylinder
#    The purpose of this is to restrict water molecules which can be used for path endpoints 
#    to those at least surface_pad nm from the surface of the cylinder. 
#    As an endpoint on the surface of the cylinder can only be approached from one side, 
#    the average best path to it is worse than the average best path through a water within the cylinder
#    this effect shows up noticeably on plots of water paths and affects the value of the progress coordinate 
#    once the protein channel is fully or mostly open.
#    we therefore calculate a smaller inner cylinder and require the endpoints to be inside of it.
# check_queries: whether to check if the chemical entities listed in the query can actually be found. 
#                Returns none if they cannot
#---------------------------------
# Returns
# cyl_o_inds: the atom indices of the oxygens of all waters in the outer cylinder
# note that these are 0-indexed, while pymol atom indices are 1-indexed
# n_water_ind: index of the water closest to the negative face of the inner cylinder in cyl_o_inds (not the atom index in the trajectory)
# p_water_ind: index of the water closest to the positive face of the inner cylinder in cyl_o_inds (not the atom index in the trajectory)
# 
# Note that if the structure wraps across the periodic boundary, 
# the water at n_water_ind may have a larger z coordinate than its counterpart at p_water_ind,
# but these atoms still maintain their positions with respect to the protein as if the structure did not wrap.

def get_water_cylinder_inds(frame, queries, r_dist, z_dist_n, z_dist_p, surface_pad, check_queries=False):

    prot_com = md.compute_center_of_mass(frame, select = "protein")[0] #the [0] is to get rid of a redundant layer of array structure

    #square outer radius for faster calculation
    r_sqdist = r_dist**2
    
    #inner cylinder dimensions
    z_dist_p_end = z_dist_p - surface_pad
    z_dist_n_end = z_dist_n - surface_pad
    r_dist_end = r_dist - surface_pad
    #square inner radius for faster calculation
    r_sqdist_end = r_dist_end**2

    #select protonatable groups; referred to as 'water' hereafter
    combined_queries = " or ".join([f"({q})" for q in queries])
    
    #check if all queries identify at least one atom
    if check_queries:
        #print(combined_queries)

        for q in queries:
            if len(frame.top.select(q)) == 0:
                print(f"error: no results for {q}")
                return 0

    #indices of protonatable atoms
    water_o = frame.top.select(combined_queries)

    #coordinates relative to protein center of mass, not PBC adjusted
    frame_xyz_prot = frame.xyz[0]-prot_com
    
    #get indices of protonatable groups in inner and outer cylinders
    #faster calculation as long as outer cylinder does not cross pbc:
    if  prot_com[0]-r_dist   > 0 and prot_com[0]+r_dist   < frame.unitcell_lengths[0][0] and \
        prot_com[1]-r_dist   > 0 and prot_com[1]+r_dist   < frame.unitcell_lengths[0][1] and \
        prot_com[2]-z_dist_n > 0 and prot_com[2]+z_dist_p < frame.unitcell_lengths[0][2]:
        
        #find water oxygens in outer cylinder
        cyl_o_inds = [ai for ai in water_o if (-z_dist_n < frame_xyz_prot[ai][2] < z_dist_p and \
                                     frame_xyz_prot[ai][0]**2 + frame_xyz_prot[ai][1]**2 < r_sqdist)]                                               

        #find water oxygens in inner cylinder for endpoints
        cyl_o_inds_end = [ai for ai in water_o if (-z_dist_n_end < frame_xyz_prot[ai][2] < z_dist_p_end and \
                                     frame_xyz_prot[ai][0]**2 + frame_xyz_prot[ai][1]**2 < r_sqdist_end)]  
        
    else:
        print("cylinder crosses PBC boundary; water wire will calculate fine but path classification will fail")
        
        
        #this is faster than running this operation only for indices ai
        image_dists = frame.unitcell_lengths[0]-abs(frame_xyz_prot)
        #The min() expression could probably be turned into a vector operation with careful array slicing. 
        #The list operations below could then be extracted from the if statement
        #alternatively if this can be made as fast as the non-pbc case just run it all the time
        
        #find water oxygens in outer cylinder
        cyl_o_inds = [ai for ai in water_o if (min(abs(frame_xyz_prot[ai][0]), image_dists[ai][0])**2 \
                                             + min(abs(frame_xyz_prot[ai][1]), image_dists[ai][1])**2 < r_sqdist \
                               and -z_dist_n < min(abs(frame_xyz_prot[ai][2]), image_dists[ai][2]) < z_dist_p)]

        #find water oxygens in inner cylinder for endpoints
        cyl_o_inds_end = [ai for ai in water_o if (min(abs(frame_xyz_prot[ai][0]), image_dists[ai][0])**2 \
                                                 + min(abs(frame_xyz_prot[ai][1]), image_dists[ai][1])**2 < r_sqdist_end \
                               and -z_dist_n_end < min(abs(frame_xyz_prot[ai][2]), image_dists[ai][2]) < z_dist_p_end)]                                               


    #get waters at either end of the smaller cylinder for path calculation
    cyl_o_z = [frame_xyz_prot[i][2] for i in cyl_o_inds_end]
    #cyl_o_r = [(frame.xyz[0][i][0]-prot_com[0])**2 + (frame.xyz[0][i][1]-prot_com[1])**2 for i in cyl_o_inds]

    #get indices of the terminal waters in the cyl_o_inds array (rather than the cyl_o_inds_end array)
    n_water_ind = cyl_o_inds.index(cyl_o_inds_end[np.argmin(cyl_o_z)])
    p_water_ind = cyl_o_inds.index(cyl_o_inds_end[np.argmax(cyl_o_z)])
        
    #print("hide spheres; show spheres, index "+"+".join([str(i+1) for i in cyl_o_inds]))
        
    return cyl_o_inds, n_water_ind, p_water_ind
    

def get_distance_matrix(frame, indices, periodic=True):

    # list of atom pairs
    index_pairs = list(itertools.product(indices, indices))
    
    #identify connections between atoms on the same residue to exclude them 
    #since hydrogen can't go directly from one glutamate O to the other
    #note that you can't just trim them out at this stage since 
    #it would mess up the shape of the distance matrix 
    intramolecular_pairs = []
    
    for i, ip in enumerate(index_pairs):
        if frame.top.atom(ip[0]).residue == frame.top.atom(ip[1]).residue:
            intramolecular_pairs.append([i//len(indices), i%len(indices)])
    
    # calculate pairwise distances
    # periodic = True is an issue for rcsb pdb structures with unit cell information if the residues are far enough apart
    distances = md.compute_distances(frame, index_pairs, periodic=periodic, opt=True).reshape(len(indices), len(indices))

    #set distances between atoms on the same residue to large values to exclude them from the MST
    for imp in intramolecular_pairs:
        distances[imp[0]][imp[1]] = 999
        distances[imp[1]][imp[0]] = 999

    return distances


def get_water_pc(frame, queries_1, r_dist, z_dist_n, z_dist_p, surface_pad, check_queries = False, printcommands = False):
    
    #get cylinder waters
    cyl_o_inds, c_water_ind, m_water_ind = get_water_cylinder_inds(frame, queries_1, r_dist, z_dist_n, z_dist_p, surface_pad, check_queries)

    #visually inspect spatial distribution of cylinder waters
    #water_dist = plt.hist2d(cyl_o_x, cyl_o_z)
    #plt.show()

    #calculate distances between the waters in the cylinder
    o_dist_mat = get_distance_matrix(frame, cyl_o_inds)
    #calculate the minimum spanning tree, which contains the minmax path between every pair of waters
    #mst is some kind of scipy sparse matrix or graph object
    mst = minimum_spanning_tree(o_dist_mat)
    #distance matrix containing the mst, with 0 used for nonexistent edges
    #note that this matrix is not symmetrical; each nonzero entry has a 0 entry at the transposed position
    #mst_arr = mst.toarray()

    #determine the shortest paths between the water with the smallest z coordinate and the water with the largest 
    dist_matrix, predecessors = shortest_path(mst, directed=False, indices=c_water_ind, return_predecessors=True)

    #traverse the predecessor list to find the shortest path
    nodelist = [m_water_ind]
    path_dists = []

    for t in range(len(predecessors)):
        if t != 0:
            path_dists.append(o_dist_mat[nodelist[-2]][nodelist[-1]])

        if predecessors[nodelist[-1]] < 0:
            break
        nodelist.append(predecessors[nodelist[-1]])
        
    #if printcommands: 
    #    print_pymol_commands(path_dists, nodelist, cyl_o_inds)
        
    path_o_inds = [str(cyl_o_inds[i]) for i in nodelist]
        
    return path_dists, nodelist, path_o_inds


#get the vertices of the hexagonal pyramid as a function of time
#--Parameters:
#  trj: an MDtraj trajectory
#--Returns:
#  a 7x3 array where each entry is the coordinates of a vertex, starting with the 'top' of the pyramid

def get_vertices(trj):
    
    #the order of these matters; if it is changed, 
    # the if statements enclosing the return statements in facecrossings() must be adjusted
    pro_resseqs = [32, 78, 132, 178, 231, 272] 
    val_thr_resseqs = [38, 138, 237]
    #could also use ARG 39, 139, and 238

    #for checking that the indices are right
    #print("+".join([str(i) for i in pro_inds+val_thr_inds]))

    val_thr_str = " or ".join([f"resSeq {r}" for r in val_thr_resseqs]) #resSeq is case sensitive
    val_thr_ca_inds = trj.top.select(f"({val_thr_str}) and name CA")

    #average valine and threonine CAs to get 'top' of pyramid
    top_vertex = np.mean(trj.xyz[-1][val_thr_ca_inds], axis=0)
    
    pro_str = " or ".join([f"resSeq {r}" for r in pro_resseqs]) #resSeq is case sensitive
    pro_ca_inds = trj.top.select(f"({pro_str}) and name CA")
    
    base_vertices = trj.xyz[-1][pro_ca_inds]
    
    return np.vstack([top_vertex, base_vertices])


def facecrossings(vertices, d, e):
        
    for i in range(6):
        
        a = vertices[0]
        b = vertices[i+1]
        c = vertices[(i+1)%6 + 1]
        
        ba = b-a
        ca = c-a
        da = d-a
        ea = e-a
        
        normvec = np.cross(ba, ca)
        on_opp_side = np.dot(normvec,da)*np.dot(normvec,ea)
    
        if on_opp_side == 0:
            #statistically should barely ever happen
            print("warning: 0 dot product; on or both points lie in plane")
        
        if on_opp_side <= 0:
           
            da_perp = np.dot(normvec, da)/np.dot(normvec,normvec)*normvec
            da_parr = da - da_perp
            ea_perp = np.dot(normvec, ea)/np.dot(normvec,normvec)*normvec
            ea_parr = ea - ea_perp

            #the point where line de intersects triangle abc
            intersection_point = (da_parr*np.sqrt(ea_perp**2)+ea_parr*np.sqrt(da_perp**2))/(np.sqrt(da_perp**2)+np.sqrt(ea_perp**2))

            #negative if b and c are on opposite sides of the line between a and the intersection
            sgnbc = np.dot(np.cross(intersection_point, ba), np.cross(intersection_point, ca))
            #negative if a and c are on opposite sides of the line between b and the intersection
            sgnca = np.dot(np.cross(intersection_point-ba, -ba), np.cross(intersection_point-ba, c-b))
        
        
            if sgnbc<=0 and sgnca<=0: #see notebook
                if i == 0 or i == 1:
                    return 0
                elif i == 2 or i == 3:
                    return 1
                else:
                    return 2
    
    return -1


def classify_path(trj, path_o_inds, vertices):
    fcs = []
    for i in range(len(path_o_inds)-1):
        #print(i)
        fc = facecrossings(vertices, trj.xyz[0][int(path_o_inds[i])], trj.xyz[0][int(path_o_inds[i+1])])
        if fc != -1:
            fcs.append(fc)
    
    if len(fcs) == 1:
        return fcs[0]
    elif len(fcs)>1:
        print("multiple crossings detected")
        print(fcs)
        return 4
    else:
        return 3
    

def calc_pc_inds(frame):
    #water oxygens and protonatable groups and OH groups on proteins
    queries_1 = [ "resname HOH and element O",
                "resname DNF and name O1"]

    queries_2 = [ "resname HOH and element O",
                "resname DNF and name O1",
                "resname GLU and (name OE1 or name OE2)",
                "resname ASP and (name OD1 or name OD2)",
                "resname SER and name OG",
                "resname THR and name OG1",
                "resname HIS and (name ND1 or name NE2)",
                "resname TYR and name OH"]

    #cylinder size
    r_dist = 2
    z_dist_p = 3.5
    z_dist_n = 4
    surface_pad = 0.4

    # for fi, frame in enumerate(trj[-1]):
    #     print("---------------------------")
    #     print(fi)
    #calculate water wire
    path_dists, nodelist, path_o_inds = get_water_pc(frame, queries_1, r_dist, z_dist_n, z_dist_p, surface_pad, check_queries = False, printcommands = False)
    
    v = get_vertices(frame)
    pt = classify_path(frame, path_o_inds, v)
    
    return [max(path_dists), pt], path_o_inds

abspath = "/wynton/home/grabe/jborowsky/aac1/westpa-26"

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
                pc, inds = calc_pc_inds(frame)
                pcs_all.append([xround, xwalker] + pc)
                wire_inds.append(inds)
            else:
                break

        #pcs_all.append(pcs)
        wire_inds_all.append(wire_inds)

        #remove unpacked archive; redundant if in case something happened to the archive while the folder was being processed
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
        
