#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script calculates McGowan Volume.
"""

from scripts import header as h

global atom_mcgowan_volumes
atom_mcgowan_volumes = {
    "H": 8.71, "He": 6.75,

    "Li": 22.23, "Be": 20.27, "B": 18.31, "C": 16.35, "N": 14.39, "O": 12.43, "F": 10.47,
    "Ne": 8.51,

    "Na": 32.71, "Mg": 30.75, "Al": 28.79, "Si": 26.83, "P": 24.87, "S": 22.91, "Cl": 20.95,
    "Ar": 18.99,

    "K": 51.89, "Ca": 50.28, "Sc": 48.68, "Ti": 47.07, "V": 45.47, "Cr": 43.86, "Mn": 42.26,
    "Fe": 40.65, "Co": 39.05, "Ni": 37.44, "Cu": 35.84, "Zn": 34.23, "Ga": 32.63, "Ge": 31.02,
    "As": 29.42, "Se": 27.81, "Br": 26.21, "Kr": 24.60,

    "Rb": 60.22, "Sr": 58.61, "Y": 57.01, "Zr": 55.40, "Nb": 53.80, "Mo": 52.19, "Tc": 50.59,
    "Ru": 48.98, "Rh": 47.38, "Pd": 45.77, "Ag": 44.17, "Cd": 42.56, "In": 40.96, "Sn": 39.35,
    "Sb": 37.75, "Te": 36.14, "I": 34.54, "Xe": 32.93,

    "Cs": 77.25, "Ba": 76.00, "La": 74.75, "Hf": 55.97, "Ta": 54.71, "W": 53.46, "Re": 52.21,
    "Os": 50.96, "Ir": 49.71, "Pt": 48.45, "Au": 47.20, "Hg": 45.95, "Tl": 44.70, "Pb": 43.45,
    "Bi": 42.19, "Po": 40.94, "At": 39.69, "Rn": 38.44,

}
def get_mcgowan_volume(m):
    # m = h.MolFromSmiles(smiles)
    m = h.AddHs(m)

    num_atoms = h.rdMolDescriptors.CalcNumAtoms(m)

    atoms = []
    bonds = []

    atoms_sum = 0
    bonds_sum = 0
    mcgowan_volume = 0

    for a, atom in enumerate(m.GetAtoms()):
        atoms.append(atom)
        # print(atom.GetSymbol(), end =", ")
        atoms_sum += atom_mcgowan_volumes[atom.GetSymbol()]

    # print(f"Atoms sum = {atoms_sum}")

    for b, bond in enumerate(m.GetBonds()):
        bonds.append(bond)
        bonds_sum += 6.56
        # print(bond.GetBondType(), end=", ")

    mcgowan_volume = atoms_sum - bonds_sum

    # print(f"Bonds sum = {bonds_sum}")
    # print(f"McGowanVolume = {mcgowan_volume}")

    return mcgowan_volume