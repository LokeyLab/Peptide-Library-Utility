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

    temp_peptide = subunits_to_bond[0].smiles_string

    for subunit in range(0, len(subunits_to_bond)):

        if subunit != (len(subunits_to_bond) - 1) and subunit != len(subunits_to_bond):
            temp_peptide = temp_peptide[0: - 1] + subunits_to_bond[subunit + 1].smiles_string

    return temp_peptide


def cyclize_head_to_tail(peptide):
    """This function cyclizes a peptide in a head to
    tail format."""

    temp_peptide = []

    temp_peptide = str(peptide[0]) + "%99" + str(peptide[1:-1]) + "%99"

    return temp_peptide

# FIXME


def cyclize_huisgen(peptide):
    """This function cyclizes a peptide in a Huisgen
    ("Click") format."""

    temp_peptide = []

    return temp_peptide


def cyclize_thioether(peptide, subunit_library):
    """This function cyclizes a peptide in a thioether
    format and required that the last subunit be a
    threonine."""

    temp_peptide = []

    if peptide.endswith(subunit_library["Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-24]) + \
                       "O[C@H](C)[C@@H](N)C(=O)%99"

    elif peptide.endswith(subunit_library["D-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-23]) + \
                       "O[C@H](C)[C@H](N)C(=O)%99"

    elif peptide.endswith(subunit_library["ß-H-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-25]) + \
                       "O[C@H](C)[C@@H](N)CC(=O)%99"

    elif peptide.endswith(subunit_library["ß-H-D-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-24]) + \
                       "O[C@H](C)[C@H](N)CC(=O)%99"

    elif peptide.endswith(subunit_library["(NMe)-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-27]) + \
                       "O[C@H](C)[C@@H](N(C))C(=O)%99"

    elif peptide.endswith(subunit_library["(NMe)-D-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-26]) + \
                       "O[C@H](C)[C@H](N(C))C(=O)%99"

    elif peptide.endswith(subunit_library["(NMe)-ß-H-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-28]) + \
                       "O[C@H](C)[C@@H](N(C))CC(=O)%99"

    elif peptide.endswith(subunit_library["(NMe)-ß-H-D-Thr"].smiles_string):

        temp_peptide = str(peptide[0]) + "%99" + \
                       str(peptide[1:-27]) + \
                       "O[C@H](C)[C@H](N(C))CC(=O)%99"

    else:

        print("\nError. Peptide may not end in a threonine.")

    return temp_peptide


def cyclize_side_chain_to_side_chain(peptide):
    """This function cyclizes a peptide in a Huisgen ("Click") format."""

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
