#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains functions for evaluating cheminformatic values using RDKits suite.
"""

from scripts import header as h

# "Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "# Atoms", "# Heteroatoms", "# Amide Bonds",
# "# Rings", "# Aromatic Rings", "# Rotatable Bonds", "# HBA", "# HBD", "Largest Ring Size"


def get_molecule_from_smiles(smiles_string):
    """Returns RDKit molecule object from a SMILES string using the MolFromSmiles function."""

    molecule = h.MolFromSmiles(smiles_string)

    molecule = h.AddHs(molecule)

    return molecule


def get_exact_mass(molecule):
    """Returns exact mass using RDKit's ExactMolWt function"""
    exact_mass = h.Descriptors.ExactMolWt(molecule)

    return exact_mass

def get_a_log_p_86(smiles):
    # h.os.environ['CLASSPATH'] = 'cdk-2.7.1.jar'
    # h.os.environ['JPYPE_JVM'] = h.jpype.getDefaultJVMPath()

    molecule = h.readstring(smiles)

    a_log_p = 0
    return a_log_p


def get_a_log_p_89(smiles):

    a_log_p = 0
    return a_log_p

def get_a_log_p_99(molecule):
    """Returns Crippen's ALogP using RDKit's MolLogP function."""

    a_log_p = h.Crippen.MolLogP(molecule)

    # molecule = h.indigo.loadMolecule(smiles)
    # a_log_p = molecule.logP()

    return a_log_p

def get_mcgowan_volume(molecule):
    """Returns McGowan Volume."""

    mcgowan_volume = h.mcgowan_volume.get_mcgowan_volume(molecule)
    return mcgowan_volume


def get_s_log_p(molecule):
    """Returns SLogP using RDKit's SlogP_VSA_ function."""

    s_log_p = h.rdMolDescriptors.SlogP_VSA_(molecule)

    return s_log_p


def get_sasa(molecule):
    """Returns SASA using RDKit's CalcSASA function"""
    h.AllChem.EmbedMolecule(molecule)
    radii = h.rdFreeSASA.classifyAtoms(molecule)
    sasa = h.rdFreeSASA.CalcSASA(molecule, radii)

    return sasa


def get_tpsa(molecule):
    """Returns TPSA using RDKit's CalcTPSA function"""

    tpsa = h.rdMolDescriptors.CalcTPSA(molecule)

    return tpsa


def get_num_atoms(molecule):
    """Returns number of atoms using RDKit's CalcNumAtoms function"""

    num_atoms = h.rdMolDescriptors.CalcNumAtoms(molecule)

    return num_atoms


def get_num_heteroatoms(molecule):
    """Returns number of atoms using RDKit's CalcNumHeteroatoms function"""

    num_heteroatoms = h.rdMolDescriptors.CalcNumHeteroatoms(molecule)

    return num_heteroatoms


def get_num_amide_bonds(molecule):
    """Returns number of amide bonds using RDKit's CalcNumAmideBonds function"""

    num_amide_bonds = h.rdMolDescriptors.CalcNumAmideBonds(molecule)

    return num_amide_bonds


def get_num_rings(molecule):
    """Returns number of rings using RDKit's CalcNumRings function"""

    num_rings = h.rdMolDescriptors.CalcNumRings(molecule)

    return num_rings


def get_num_aromatic_rings(molecule):
    """Returns number of aromatic rings using RDKit's CalcNumAromaticRings function"""

    num_aromatic_rings = h.rdMolDescriptors.CalcNumAromaticRings(molecule)

    return num_aromatic_rings


def get_num_rotatable_bonds(molecule):
    """Returns number of rotatable bonds using RDKit's CalcNumRotatableBonds function"""

    num_rotatable_bonds = h.rdMolDescriptors.CalcNumRotatableBonds(molecule)

    return num_rotatable_bonds


def get_num_hba(molecule):
    """Returns number of hydrogen bond acceptors using RDKit's CalcNumHBA function"""

    num_hba = h.rdMolDescriptors.CalcNumHBA(molecule)

    return num_hba


def get_num_hbd(molecule):
    """Returns number of hydrogen bond donors using RDKit's CalcNumHBD function"""

    num_hbd = h.rdMolDescriptors.CalcNumHBD(molecule)

    return num_hbd


def get_largest_ring_size(molecule):
    """Returns number of rotatable bonds using RDKit's CalcNumRotatableBonds function"""

    ri = molecule.GetRingInfo()
    largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)

    return largest_ring_size


def get_chemometrics(molecule):
    """Returns the exact mass, tpsa, and ALogP which is
    calculated using RDKit without any force field
    optimization"""

    mol_h = h.AddHs(molecule)

    exact_mass = get_exact_mass(mol_h)
    tpsa = get_tpsa(mol_h)
    a_log_p = get_a_log_p(mol_h)

    return exact_mass, tpsa, a_log_p


def get_chemometrics_with_uff_optimization(molecule):
    """Returns the exact mass, tpsa, and ALogP which is
    calculated using RDKit with UFF force field
    optimization"""

    molecule_h = h.AddHs(molecule)

    h.AllChem.EmbedMolecule(molecule_h)

    h.AllChem.UFFOptimizeMolecule(molecule_h)

    exact_mass = get_exact_mass(molecule_h)
    tpsa = get_tpsa(molecule_h)
    a_log_p = get_a_log_p(molecule_h)

    return exact_mass, tpsa, a_log_p


def get_chemometrics_with_mmff_optimization(molecule):
    """Returns the exact mass, tpsa, and ALogP which is
    calculated using RDKit with MMFF force field
    optimization"""

    molecule_h = h.AddHs(molecule)

    h.AllChem.EmbedMolecule(molecule_h)

    h.AllChem.MMFFOptimizeMolecule(molecule_h)

    exact_mass = get_exact_mass(molecule_h)
    tpsa = get_tpsa(molecule_h)
    a_log_p = get_a_log_p(molecule_h)

    return exact_mass, tpsa, a_log_p
