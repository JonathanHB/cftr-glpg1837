#Jonathan Borowsky
#Grabe lab
#3/2/26
##########################################################

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import distances



def tmd_query():
    segment_resis = [[77, 149], [192, 245], [298, 362], [988, 1034], [857, 889], [900, 942], [1094, 1154]]
    print("color blue, " + " or ".join([f"resi {sr[0]}-{sr[1]}" for sr in segment_resis]))

    segment_resis_all = [i for sr in segment_resis for i in range(sr[0], sr[1]+1)]
    query = " or ".join([f"resid {sr}" for sr in segment_resis_all])

    #indices = frame.top.select(f"protein and ({query})")

    # print_indices = frame.top.select(f"protein and name CA and ({query})")
    # print("+".join([str(i+1) for i in print_indices]))
    return f"protein and ({query})"


def contacts_bin(group_1, group_2, cutoff=5.0):
    dists = distances.distance_array(group_1, group_2, box=u.dimensions)
    contacts = np.where(dists < cutoff, 1, 0)
    contact_bin = np.where(np.sum(contacts, axis=1) > 0, 1, 0)

    return contact_bin


def select_heavy_atoms(u):
    """Select heavy atoms based on a given selection string."""

    # ============================
    # SELECTION PARAMETERS
    # ============================

    # List of protein residue numbers (edit as needed)
    #protein_resids = [10, 25, 48]   # <-- USER FILLS THIS
    protein_atoms = f"{tmd_query()} and not name H*"
    #protein_atoms = "resid 207 and resname GLN and not name H*"

    # Ligand residue name
    ligand_resname = "LJP"          # <-- USER FILLS THIS

    # Lipid phosphate atom name (common choices: "P", "PO4", etc.)
    phosphate_selection = "(resname PA or resname OL) and name C12"  #"resname PC and name P31"

    #membrane_core_buffer = 0


    # ============================
    # PROTEIN HEAVY ATOMS
    # ============================

    #resid_string = " ".join(str(r) for r in protein_resids)

    protein_sel = u.select_atoms(protein_atoms)

    # ============================
    # LIGAND HEAVY ATOMS
    # ============================

    ligand_sel = u.select_atoms(
        f"resname {ligand_resname} and not name H*"
    )

    # ============================
    # MEMBRANE PHOSPHATES
    # ============================

    phosphates = u.select_atoms(phosphate_selection)

    if len(phosphates) == 0:
        raise ValueError("No phosphate atoms found. Check lipid resnames and phosphate atom name.")

    z_positions = phosphates.positions[:, 2]
    z_mean = np.mean(z_positions)

    # Separate into upper and lower leaflets
    upper_leaflet = phosphates[z_positions > z_mean]
    lower_leaflet = phosphates[z_positions <= z_mean]

    if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
        raise ValueError("Leaflet separation failed.")

    z_upper_avg = np.mean(upper_leaflet.positions[:, 2])
    z_lower_avg = np.mean(lower_leaflet.positions[:, 2])

    z_min = min(z_upper_avg, z_lower_avg)
    z_max = max(z_upper_avg, z_lower_avg)

    # ============================
    # WATERS INSIDE MEMBRANE
    # ============================

    waters = u.select_atoms("resname TP3 and name O")

    water_z = waters.positions[:, 2]

    inside_mask = (water_z > z_min) & (water_z < z_max)
    waters_inside = waters[inside_mask]

    #print(np.where(inside_mask))
    #water_in_inds = np.where(inside_mask)[0]
    print("color blue, index " + "+".join([str(i+1) for i in waters_inside.indices]))
    # ============================
    # SUMMARY
    # ============================

    print(f"Protein heavy atoms: {len(protein_sel)}")
    print(f"Ligand heavy atoms: {len(ligand_sel)}")
    print(f"Waters inside membrane (heavy atoms): {len(waters_inside)}")

    # ============================
    # CALCULATE CONTACTS
    # ============================

    #water contacts are averaged across waters so we get a 1d contact array
    protein_water_contacts = contacts_bin(protein_sel, waters_inside, cutoff=5.0)
    ligand_water_contacts = contacts_bin(ligand_sel, waters_inside, cutoff=5.0)

    #protein-ligand contacts are not averaged so we get a 2d contact array
    protein_ligand_dists = distances.distance_array(protein_sel,
                             ligand_sel,
                             box=u.dimensions)

    protein_ligand_contacts = np.where(protein_ligand_dists < 5, 1, 0)

    return waters_inside.positions, protein_water_contacts, ligand_water_contacts, protein_ligand_contacts, protein_sel, ligand_sel




def main(ref_frame_path, gro_file, xtc_file):
    # ============================
    # LOAD TRAJECTORY
    # ============================

    #ref_frame_path_x01 = '/media/X01Raid01/Data_Backup/home/csheen/cftr-project/wstp_cftr_1_degrabo/bstates/input/min.gro'
    #ref_frame_path = '/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/topology/input.gro'

    ref_frame = mda.Universe(ref_frame_path)

    #gro_file = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/topology/input.gro"
    #xtc_file = "/home/jonathan/Documents/grabelab/cftr/independent-partial-dissociation/nonlip_glpg_1/001913-000187-trj-pbcmol-centered-tmd-rot-s10.xtc"

    # Frame index to analyze
    frame_index = 0

    u = mda.Universe(gro_file, xtc_file)

    align.AlignTraj(u, ref_frame, select=f"{tmd_query()} and name CA")

    wat_pos, pwc, lwc, plc, prot_sel, lig_sel = select_heavy_atoms(u)

    return wat_pos, pwc, lwc, plc
    
