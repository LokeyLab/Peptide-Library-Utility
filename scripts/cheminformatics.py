#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains functions for evaluating cheminformatic values using RDKits suite.
"""

from scripts import header as h


def get_molecule_from_smiles(smilesString):
    """Returns RDKit molecule object from a SMILES string using the MolFromSmiles function."""
    
    molecule = h.MolFromSmiles(smilesString)

    return molecule


def get_exact_mass(molecule):
    """Returns exact mass using RDKit's ExactMolWt function"""
    
    exact_mass = h.Descriptors.ExactMolWt(molecule)

    return exact_mass


def get_tpsa(molecule):
    """Returns TPSA using RDKit's CalcTPSA function"""
    
    tpsa = h.rdMolDescriptors.CalcTPSA(molecule)

    return tpsa


def Geta_log_p(molecule):
    """Returns Crippen's ALogP using RDKit's MolLogP function."""
    
    a_log_p = h.Crippen.MolLogP(molecule)

    return a_log_p


def get_chemometrics(molecule):
    """Returns the exact mass, tpsa, and ALogP which is calculated using RDKit without any force field optimization"""

    mol_h = h.AddHs(molecule)

    exact_mass = get_exact_mass(mol_h)
    tpsa = get_tpsa(mol_h)
    a_log_p = Geta_log_p(mol_h)

    return exact_mass, tpsa, a_log_p


def get_chemometrics_with_uff_optimization(molecule):
    """Returns the exact mass, tpsa, and ALogP which is calculated using RDKit with UFF force field optimization"""

    moleculeH = h.AddHs(molecule)

    h.AllChem.EmbedMolecule(moleculeH)

    h.AllChem.UFFOptimizeMolecule(moleculeH)

    exact_mass = get_exact_mass(moleculeH)
    tpsa = get_tpsa(moleculeH)
    a_log_p = Geta_log_p(moleculeH)

    return exact_mass, tpsa, a_log_p


def get_chemometrics_with_mmff_optimization(molecule):
    """Returns the exact mass, tpsa, and ALogP which is calculated using RDKit with MMFF force field optimization"""

    moleculeH = h.AddHs(molecule)

    h.AllChem.EmbedMolecule(moleculeH)

    h.AllChem.MMFFOptimizeMolecule(moleculeH)

    exact_mass = get_exact_mass(moleculeH)
    tpsa = get_tpsa(moleculeH)
    a_log_p = Geta_log_p(moleculeH)

    return exact_mass, tpsa, a_log_p
