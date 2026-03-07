import sys
import os
import numpy as np
import time
import MDAnalysis as mda
from MDAnalysis.analysis import align

import matplotlib.pyplot as plt

#import utility

from scipy.spatial.distance import cdist

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


def tmd_query():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    #print("color blue, " + " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis]))

    segment_resis_all = [i for sr in segment_resis for i in range(sr[0], sr[1]+1)]
    query = " or ".join([f"resid {sr}" for sr in segment_resis_all])

    # print_indices = frame.top.select(f"protein and name CA and ({query})")
    # print("+".join([str(i+1) for i in print_indices]))
    return f"protein and ({query})"


def load_frame(top_path, xtc_path=None, frame_index = -1):

    if xtc_path is None:
        u = mda.Universe(top_path)
    else:
        u = mda.Universe(top_path, xtc_path) #, topology_format='ITP'

    u.trajectory[frame_index]

    # create a new universe from that frame
    u_last = mda.Merge(u.atoms)

    # copy box, which for some reason merge does not do
    u_last.dimensions = u.dimensions

    return u_last


#these observables are all best calculated from a non-aligned trajectory since they depend on z
def get_waters_phosphates_boxvectors(u):

    # ============================
    # MEMBRANE PHOSPHATES
    # ============================

    phosphates = u.select_atoms("resname PC and name P31")

    if len(phosphates) == 0:
        raise ValueError("No phosphate atoms found. Check lipid resnames and phosphate atom name.")

    #print(phosphates.positions.shape)

    z_positions = phosphates.positions[:, 2]
    z_mean = np.mean(z_positions)

    # Separate into upper and lower leaflets
    upper_leaflet = phosphates[z_positions > z_mean]
    lower_leaflet = phosphates[z_positions <= z_mean]

    if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
        raise ValueError("Leaflet separation failed.")

    z_upper_avg = np.mean(upper_leaflet.positions[:, 2])
    z_lower_avg = np.mean(lower_leaflet.positions[:, 2])

    membrane_thickness = z_upper_avg-z_lower_avg

    # ============================
    # WATERS WITHIN MEMBRANE
    # ============================

    waters = u.select_atoms("resname TP3 and name O")

    water_z = waters.positions[:, 2]
    z_pad = 0

    inside_mask = (water_z > z_lower_avg-z_pad) & (water_z < z_upper_avg+z_pad)
    waters_inside = waters[inside_mask]

    return u.dimensions, membrane_thickness, waters_inside


def contacts_bin(group_1, group_2, cutoff=5.0, flat=False):

    sqdistances = cdist(group_1.positions, group_2.positions, 'sqeuclidean')

    #identify contacts
    contacts = np.where(sqdistances < cutoff**2, 1, 0)

    #calculate whether each atom of group 1 contacts any atom of group 2
    if flat:
        contacts = np.where(np.sum(contacts, axis=1) > 0, 1, 0)

    return contacts


#align trajectory and calculate aligned water positions and various contacts
def get_water_contacts(u, ref, water_sel):

    align.AlignTraj(u, ref, select=f"{tmd_query()} and name CA", in_memory=True).run()

    #protein, ligand, and lipid selections
    protein_sel = u.select_atoms(f"{tmd_query()} and not name H*")
    ligand_sel  = u.select_atoms(f"resname LJP and not name H*")
    lipid_sel   = u.select_atoms(f"(resname PA or resname OL) and not name H*")
    pyr_hb_donor = u.select_atoms(f"resname LJP and name N2")
    pyr_hb_acceptor = u.select_atoms(f"resname LJP and name N3")
    R933_hb_donor = u.select_atoms(f"resname ARG and resid 933 and name NE NH1 NH2")
    E873_hb_acceptor = u.select_atoms(f"resname GLU and resid 873 and name OE1 OE2")
    
    if False:
        print(f"{len(protein_sel)} protein atoms")
        print(f"{len(ligand_sel)} ligand atoms")
        print(f"{len(lipid_sel)} lipid atoms")
        print(f"{len(water_sel)} water atoms")

    prot_lig = contacts_bin(protein_sel, ligand_sel, cutoff=5.0, flat=False)
    prot_lip = contacts_bin(protein_sel, lipid_sel, cutoff=5.0, flat=True)
    prot_wat = contacts_bin(protein_sel, water_sel, cutoff=5.0, flat=True)
    lig_lip = contacts_bin(ligand_sel, lipid_sel, cutoff=5.0, flat=True)
    lig_wat = contacts_bin(ligand_sel, water_sel, cutoff=5.0, flat=True)

    hb_ljp_donor = contacts_bin(pyr_hb_donor, E873_hb_acceptor, cutoff=3.3, flat=True)[0]
    hb_ljp_acceptor = contacts_bin(pyr_hb_acceptor, R933_hb_donor, cutoff=3.3, flat=True)[0]
    lig_saltb_hbonds = np.array([hb_ljp_donor,hb_ljp_acceptor])

    return protein_sel, ligand_sel, water_sel.positions, ligand_sel.positions, prot_lig, prot_lip, prot_wat, lig_lip, lig_wat, lig_saltb_hbonds


def get_n_observables():
    return 10


def main(ref_path, gro_path, xtc_path, frame=-1):
    u = load_frame(gro_path, xtc_path, frame)
    boxdims, thickness, waters_inside = get_waters_phosphates_boxvectors(u)
    ref = mda.Universe(ref_path)
    protein_sel, ligand_sel, water_pos, lig_pos, prot_lig, prot_lip, prot_wat, lig_lip, lig_wat, lig_saltb_hbonds = get_water_contacts(u, ref, waters_inside)

    output = [water_pos, lig_pos, prot_lig, prot_lip, prot_wat, lig_lip, lig_wat, lig_saltb_hbonds, boxdims, thickness]

    if len(output) != get_n_observables():
        print(f"error: incorrect number of returned observables: {len(output)} vs {get_n_observables()}")
        sys.exit(0)

    return output

def observable_names():
    return [
        "water_pos",
        "ligand_pos",
        "prot_lig_contacts",
        "prot_lip_contacts",
        "prot_wat_contacts",
        "lig_lip_contacts",
        "lig_wat_contacts",
        "lig_saltbridge_hbonds",
        "boxdims",
        "membrane_thickness"
        ]