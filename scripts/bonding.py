#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains all the functions that are used in
the bonding process of subunits to others as well as the
functions that are used for cyclization.
"""


def bond_head_to_tail(subunits_to_bond):
    """This functions bonds the subunits in a list of
    subunits in a head to tail format to create a
    peptide."""

    temp_peptide = subunits_to_bond[len(subunits_to_bond) - 1].smiles_string

    for subunit in range(len(subunits_to_bond) - 1, 0, -1):

        if subunit >= (1):
            temp_peptide = temp_peptide[0: - 1] + subunits_to_bond[subunit - 1].smiles_string

    return temp_peptide


def cyclize_head_to_tail(peptide):
    """This function cyclizes a peptide in a head to
    tail format."""

    temp_peptide = []

    temp_peptide = str(peptide[0]) + "%99" + str(peptide[1:-1]) + "%99"

    return temp_peptide


def cyclize_triazole(peptide, subunit_library):
    """This function cyclizes a peptide in a trizaole format."""

    temp_peptide = []

    if peptide.endswith(subunit_library["Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["Prg"].smiles_string)]) + \
                       "N[C@@H](CC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("Prg")

    elif peptide.endswith(subunit_library["ß-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["ß-Prg"].smiles_string)]) + \
                       "N[C@H](CC%98=C-N(CC(=O)%99)-N=N%98)CC(=O)O"
        print("ß-Prg")

    elif peptide.endswith(subunit_library["H-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["H-Prg"].smiles_string)]) + \
                       "N[C@H](CCC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("H-Prg")

    elif peptide.endswith(subunit_library["D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["D-Prg"].smiles_string)]) + \
                       "N[C@H](CC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("D-Prg")
    elif peptide.endswith(subunit_library["ß-D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["ß-D-Prg"].smiles_string)]) + \
                       "N[C@@H](CC%98=C-N(CC(=O)%99)-N=N%98)CC(=O)O"
        print("ß-D-Prg")

    elif peptide.endswith(subunit_library["H-D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["H-D-Prg"].smiles_string)]) + \
                       "N[C@@H](CCC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("H-D-Prg")

    elif peptide.endswith(subunit_library["(NMe)-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-Prg"].smiles_string)]) + \
                       "N(C)[C@H](CC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("(NMe)-Prg")

    elif peptide.endswith(subunit_library["(NMe)-ß-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-ß-Prg"].smiles_string)]) + \
                       "N(C)[C@H](CC%98=C-N(CC(=O)%99)-N=N%98)CC(=O)O"
        print("(NMe)-ß-Prg")

    elif peptide.endswith(subunit_library["(NMe)-H-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-H-Prg"].smiles_string)]) + \
                       "N(C)[C@H](CCC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("(NMe)-H-Prg")

    elif peptide.endswith(subunit_library["(NMe)-D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-D-Prg"].smiles_string)]) + \
                       "N(C)[C@H](CC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("(NMe)-D-Prg")

    elif peptide.endswith(subunit_library["(NMe)-ß-D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-ß-D-Prg"].smiles_string)]) + \
                       "N(C)[C@@H](CC%98=C-N(CC(=O)%99)-N=N%98)CC(=O)O"
        print("(NMe)-ß-D-Prg")

    elif peptide.endswith(subunit_library["(NMe)-H-D-Prg"].smiles_string):
        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-len(subunit_library["(NMe)-H-D-Prg"].smiles_string)]) + \
                       "N(C)[C@@H](CCC%98=C-N(CC(=O)%99)-N=N%98)C(=O)O"
        print("(NMe)-H-D-Prg")

    return temp_peptide


def cyclize_depsipeptide(peptide, subunit_library):
    """This function cyclizes a peptide in a macrolactone
    format and required that the last subunit be a
    threonine."""

    temp_peptide = []

    if peptide.startswith(subunit_library["Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["ß-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["H-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["D-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["ß-D-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["H-D-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-ß-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-H-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-D-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-ß-D-Thr"].smiles_string[:-1]) \
            or peptide.startswith(subunit_library["(NMe)-H-D-Thr"].smiles_string[:-1]):

        oh_pos = peptide.index("H](O)") + 4

        temp_peptide = str(peptide[0:oh_pos]) + "%99" + str(peptide[oh_pos:-1]) + "%99"

    return temp_peptide


def cyclize_thioether(peptide, subunit_library):
    """This function cyclizes a peptide in a thioether
    format and required that the last subunit be a
    threonine."""

    temp_peptide = []

    if peptide.endswith(subunit_library["Cys"].smiles_string) \
            or peptide.endswith(subunit_library["ß-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["H-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["D-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["ß-D-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["H-D-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-ß-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-H-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-D-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-ß-D-Cys"].smiles_string) \
            or peptide.endswith(subunit_library["(NMe)-H-D-Cys"].smiles_string):

        sulfur_pos = peptide.rfind("S)") + 1

        temp_peptide = peptide[0] + "%99" + peptide[1:sulfur_pos] + "CC(=O)%99" + peptide[
                                                                                 sulfur_pos:]
    print(temp_peptide)
    return temp_peptide


def cyclize_side_chain_to_side_chain(peptide):
    """This function cyclizes a peptide in a side chain to side chain format."""

    temp_peptide = []

    return temp_peptide


def cyclize_head_to_side_chain(peptide):
    """This function cyclizes a peptide in a head to side chain format."""

    temp_peptide = []

    return temp_peptide


def cyclize_tail_to_side_chain(peptide):
    """This function cyclizes a peptide in a tail to side chain format."""

    temp_peptide = []

    return temp_peptide
