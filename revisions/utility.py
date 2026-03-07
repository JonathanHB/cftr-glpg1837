import numpy as np
import os
import MDAnalysis as mda

#Chatgpt prompt:
"""write a python function that takes as its arguments a path to a pdb file, an mdanalysis atomgroup, and a list of numbers equal in length to the atom group, 
and saves a copy of the pdb file with the b factors of the atoms in the atom group set equal to the numbers in the list. 
Match atoms using residue names and numbers since the indices do not match between the atomgroup and pdb file"""


def write_bfactors_by_residue_match(pdb_path, atomgroup, values, output_pdb=None):
    """
    Write a copy of a PDB file with B-factors set for atoms matching an AtomGroup.
    Matching is done using (resname, resid, atom name), NOT indices.

    Parameters
    ----------
    pdb_path : str
        Path to input PDB file.
    atomgroup : MDAnalysis AtomGroup
        AtomGroup derived from (possibly different) structure.
    values : array-like
        Numbers equal in length to atomgroup.
    output_pdb : str
        Path to output PDB file.
    """

    if output_pdb is None:
        base = os.path.splitext(pdb_path)[0]
        output_pdb = base + "_bfactored.pdb"

    values = np.asarray(values)

    if len(atomgroup) != len(values):
        raise ValueError("Length of values must equal AtomGroup length.")

    # Load fresh universe from PDB to preserve original ordering
    u = mda.Universe(pdb_path)

    # Ensure tempfactors exist
    if not hasattr(u.atoms, "tempfactors"):
        u.add_TopologyAttr("tempfactors")

    # Start from existing B-factors (or zeros)
    try:
        b = u.atoms.tempfactors.copy()
    except AttributeError:
        b = np.zeros(len(u.atoms))

    # Build lookup dictionary from PDB universe
    pdb_lookup = {}
    for atom in u.atoms:
        key = (atom.resname, atom.resid, atom.name)
        pdb_lookup[key] = atom.index

    # Assign B-factors by matching keys
    for atom, value in zip(atomgroup, values):
        key = (atom.resname, atom.resid, atom.name)

        if key not in pdb_lookup:
            raise ValueError(
                f"Atom not found in PDB: resname={atom.resname}, "
                f"resid={atom.resid}, name={atom.name}"
            )

        pdb_index = pdb_lookup[key]
        b[pdb_index] = value

    # Reassign full array (required by MDAnalysis)
    u.atoms.tempfactors = b

    # Write output
    with mda.Writer(output_pdb, multiframe=False) as W:
        W.write(u.atoms)
