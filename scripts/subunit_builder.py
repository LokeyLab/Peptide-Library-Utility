#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains the classes that are used in the various stages of creating peptide libraries. Throughout the
process, data is held in objects in order to easily be able to the different associate attributes of the large
quantity of subunits that might be used in creating the peptide libraries. This script also contains the various
functions that are used in the creation of the various subunits.
"""

from scripts import header as h


class AminoAcid:
    """Creates an amino acid object that holds a name, multiple letter acronym, the SMILES string for the amino group,
    side chain, stereocenter, carboxyl group, and finally, the complete smiles string."""

    def __init__(self, name, multiple_letter, amino_group, side_chain, stereocenter, carboxyl_group):
        """Initializes the attributes of the amino acid object."""

        self.name = name
        self.multiple_letter = multiple_letter
        self.amino_group = amino_group
        self.side_chain = side_chain
        self.stereocenter = stereocenter
        self.carboxyl_group = carboxyl_group
        self.smiles_string = (amino_group + stereocenter + side_chain + "C" + carboxyl_group).rstrip()


class Peptoid:
    """Creates a peptoid object that holds a name, multiple letter acronym, and the complete smiles string."""

    def __init__(self, name, multiple_letter, smiles_string):
        """Initializes the attributes of the peptoid object."""

        self.name = name
        self.multiple_letter = multiple_letter
        self.smiles_string = smiles_string.rstrip()


class Misc:
    """Creates a miscellaneous object that holds a name, multiple letter acronym, and the complete smiles string."""

    def __init__(self, name, multiple_letter, smiles_string):
        """Initializes the attributes of the miscellaneous object."""

        self.name = name
        self.multiple_letter = multiple_letter
        self.smiles_string = smiles_string.rstrip()


def generate_base_amino_acid(parentheses, name, multiple_letter, amino_group, side_chain, carboxyl_group):
    """Generates a base amino acid object for the canonical amino acids without stereochemistry from which the
    stereoisomers, and homo-, beta-, homo-beta, N-methyl, etc., variants can be derived from."""

    # Handles cases where there are either parentheses needed or not.
    if parentheses == "y":

        base_amino_acid = h.classes.AminoAcid(name, multiple_letter, amino_group, "(" + side_chain + ")", "C",
                                              carboxyl_group)

    elif parentheses == "n":

        base_amino_acid = h.classes.AminoAcid(name, multiple_letter, amino_group, side_chain, "C",
                                              carboxyl_group)

    else:

        print("Parentheses not specified in spread sheet!")

    print("\n" + base_amino_acid.name.capitalize() + "s")
    print("----------------------------------------------------------")

    return base_amino_acid


def keep_base_variant(base_amino_acid):
    """Keeps th base variant in the case where stereochemistry is not needed."""

    base_variant = base_amino_acid

    print(base_variant.name + "     " + base_variant.multiple_letter + "     " + base_variant.smiles_string)

    return base_variant


def generate_l_variants(base_amino_acid):
    """Generates the L- variants."""

    l_variant = h.classes.AminoAcid(base_amino_acid.name, base_amino_acid.multiple_letter,
                                    base_amino_acid.amino_group, base_amino_acid.side_chain, "[C@@H]",
                                    base_amino_acid.carboxyl_group)

    print(l_variant.name + "     " + l_variant.multiple_letter + "     " + l_variant.smiles_string)

    return l_variant


def generate_d_variants(base_amino_acid):
    """Generates the D- variants."""

    d_variant = AminoAcid("D-" + base_amino_acid.name, "D-" + base_amino_acid.multiple_letter,
                          base_amino_acid.amino_group, base_amino_acid.side_chain, "[C@H]",
                          base_amino_acid.carboxyl_group)

    print(d_variant.name.capitalize() + "     " + d_variant.multiple_letter + "     " + d_variant.smiles_string)

    return d_variant


# FIXME


""" 
def GenerateBetaVariants(base_amino_acid):
    
    side_chain = base_amino_acid.side_chain
    betaside_chain = base_amino_acid.side_chain[0] + base_amino_acid.side_chain[2:]

    betaVariant = AminoAcid("ß-" + base_amino_acid.name, "ß-" + base_amino_acid.multiple_letter, base_amino_acid.amino_group,
                            betaside_chain, base_amino_acid.stereocenter, "C" + base_amino_acid.carboxyl_group)

    print(betaVariant.name + "      " + betaVariant.smiles_string)

    return betaVariant


def GenerateHomoVariants(base_amino_acid):
    side_chain = base_amino_acid.side_chain
    homoside_chain = base_amino_acid.side_chain[0] + "C" + base_amino_acid.side_chain[1:]

    homoVariant = AminoAcid("homo-" + base_amino_acid.name, "H-" + base_amino_acid.multiple_letter, base_amino_acid.amino_group,
                         homoside_chain, base_amino_acid.stereocenter, base_amino_acid.carboxyl_group)

    print(homoVariant.name + "      " + homoVariant.smiles_string)

    return homoVariant

"""


def generate_beta_homo_variants(base_amino_acid):
    """Generates the beta- homo- variants."""

    if base_amino_acid.name.startswith("D-"):

        beta_homo_variant = AminoAcid("ß-homo-" + base_amino_acid.name, "ß-H-" + base_amino_acid.multiple_letter,
                                      base_amino_acid.amino_group,
                                      base_amino_acid.side_chain, "[C@H]",
                                      "C" + base_amino_acid.carboxyl_group)

    else:

        beta_homo_variant = AminoAcid("ß-homo-" + base_amino_acid.name, "ß-H-" + base_amino_acid.multiple_letter,
                                      base_amino_acid.amino_group,
                                      base_amino_acid.side_chain, "[C@@H]",
                                      "C" + base_amino_acid.carboxyl_group)

    print(beta_homo_variant.name + "      " + beta_homo_variant.multiple_letter + "     " +
          beta_homo_variant.smiles_string)

    return beta_homo_variant


def generate_nme_variants(base_amino_acid):
    """Generates the N-methylated variants."""

    nme_variant = AminoAcid("NMe-" + base_amino_acid.name, "(NMe)-" + base_amino_acid.multiple_letter,
                            base_amino_acid.amino_group + "(C)", base_amino_acid.side_chain,
                            base_amino_acid.stereocenter, base_amino_acid.carboxyl_group)

    print(nme_variant.name + "       " + nme_variant.multiple_letter + "     " + nme_variant.smiles_string)

    return nme_variant


def generate_amino_acids(df):
    """Generates all the amino acids according to the specifications in the canonical and noncanonical
    amino acids .csv."""

    amino_acid_list = []

    for i in range(0, len(df)):

        temp_amino_acid_list = []

        data = df.loc[i]

        name = data.loc["Name"]
        multiple_letter = data.loc["Multiple Letter"]
        amino_group = data.loc["Amino Group"]
        side_chain = data.loc["Side Chain"]
        carboxyl_group = data.loc["Carboxyl Group"]
        parentheses = data.loc["Side Chain Parentheses? y/n"]
        l_and_d = data.loc["L- and D-? y/n"]
        beta_homo = data.loc["ß-Homo-? y/n"]
        nme = data.loc["NMe? y/n"]

        # Special case for the prolines. ß AND HOMO TO BE ADDED.
        if multiple_letter == "Pro":

            print("\nProlines")

            print("----------------------------------------------------------")

            amino_acid = AminoAcid(name, multiple_letter, amino_group, side_chain, "", carboxyl_group)

            temp_amino_acid_list += [amino_acid]

            print(amino_acid.name + "       " + amino_acid.multiple_letter + "     " + amino_acid.smiles_string)

            amino_acid_list += temp_amino_acid_list

            temp_amino_acid_list.clear()

        elif multiple_letter == "ß-H-Pro":

            amino_acid = AminoAcid(name, multiple_letter, amino_group, side_chain, "", "C" + carboxyl_group)

            temp_amino_acid_list += [amino_acid]

            print(amino_acid.name + "       " + amino_acid.multiple_letter + "     " + amino_acid.smiles_string)

            amino_acid_list += temp_amino_acid_list

            temp_amino_acid_list.clear()

        elif multiple_letter == "D-Pro":

            amino_acid = AminoAcid(name, multiple_letter, amino_group, side_chain, "", carboxyl_group)

            temp_amino_acid_list += [amino_acid]

            print(amino_acid.name + "       " + amino_acid.multiple_letter + "     " + amino_acid.smiles_string)

            amino_acid_list += temp_amino_acid_list

            temp_amino_acid_list.clear()

        elif multiple_letter == "ß-H-D-Pro":

            amino_acid = AminoAcid(name, multiple_letter, amino_group, side_chain, "", "C" + carboxyl_group)

            temp_amino_acid_list += [amino_acid]

            print(amino_acid.name + "       " + amino_acid.multiple_letter + "     " + amino_acid.smiles_string)

            amino_acid_list += temp_amino_acid_list

            temp_amino_acid_list.clear()

        # Case for the rest of the canonical amino acids
        else:

            base_amino_acid = generate_base_amino_acid(parentheses, name, multiple_letter, amino_group,
                                                       side_chain, carboxyl_group)

            if l_and_d == "y":

                temp_amino_acid_list.append(generate_l_variants(base_amino_acid))
                temp_amino_acid_list.append(generate_d_variants(base_amino_acid))

            else:

                temp_amino_acid_list.append(keep_base_variant(base_amino_acid))

            if beta_homo == "y":

                beta_homoamino_acid_list = []

                for amino_acid in temp_amino_acid_list:
                    beta_homoamino_acid_list.append(generate_beta_homo_variants(amino_acid))

                temp_amino_acid_list += beta_homoamino_acid_list

            if nme == "y":

                nmeamino_acid_list = []

                for amino_acid in temp_amino_acid_list:
                    nmeamino_acid_list.append(generate_nme_variants(amino_acid))

                temp_amino_acid_list += nmeamino_acid_list

            amino_acid_list += temp_amino_acid_list

    print("----------------------------------------------------------")

    return amino_acid_list


def generate_peptoids(df):
    """Generates the peptoids according to the specifications in the amines.csv. This function automatically binds
    the amines to bromoacetic acid in order to create a peptoid. This function and the proper creation of the peptoids
    in general rely heavily on the proper formatting of the SMILES string in the amines.csv so that they can be bound 
    to the bromoacetic acid and work properly in the peptide bonding process."""

    peptoid_list = []

    for i in range(0, len(df)):

        data = df.loc[i]

        name = str(data.loc["Name"])
        multiple_letter = str(data.loc["Multiple Letter"])
        smiles_string = str(data.loc["SMILES String"])
        bromoacetic_acid = "BrCC(=O)O"

        if h.pd.isnull(data.loc["Multiple Letter"]):

            peptoid = Peptoid(name, "       ", "N(" + smiles_string[1:] + ")" + bromoacetic_acid[2:])

        else:

            peptoid = Peptoid(name, multiple_letter, "N(" + smiles_string[1:] + ")" + bromoacetic_acid[2:])

        print("Bromoacetic Acid + " + peptoid.name + "     " + peptoid.multiple_letter + "     " +
              peptoid.smiles_string)

        peptoid_list += [peptoid]

    print("----------------------------------------------------------")

    return peptoid_list


def generate_miscellaneous(df):
    """This function generates on object that holds any miscellaneous subunit that might be added to the misc.csv.
    This function heavily depends on the user inputting the SMILES string properly so that it behaves the way that 
    they want it to. If the outcome the user is expecting is not happening, they should check their SMILES string
    input."""

    miscellaneous_list = []

    for i in range(0, len(df)):

        data = df.loc[i]

        name = str(data.loc["Name"])
        multiple_letter = str(data.loc["Multiple Letter"])
        smiles_string = str(data.loc["SMILES String"])

        if h.pd.isnull(data.loc["Multiple Letter"]):

            misc = Misc(name, "       ", smiles_string)

        else:

            misc = Misc(name, multiple_letter, smiles_string)

        print(misc.name + "     " + misc.multiple_letter + "     " + misc.smiles_string)

        miscellaneous_list += [misc]

    print("----------------------------------------------------------")

    return miscellaneous_list


def generate_subunit_library():
    """This is the main function that calls all the other functions that create the various subunits from the
    parameters defined in their respective .csv files."""

    subunits = {}

    canonical_amino_acid_dataframe = h.utilities.csv_to_dataframe(
        "subunits/Canonical Amino Acids-Canonical Amino Acids.csv")
    noncanonical_amino_acid_dataframe = h.utilities.csv_to_dataframe(
        "subunits/Non-Canonical Amino Acids-Non-Canonical Amino Acids.csv")
    amine_dataframe = h.utilities.csv_to_dataframe(
        "subunits/Amines-Amines.csv")
    miscellaneous_dataframe = h.utilities.csv_to_dataframe(
        "subunits/Misc-Misc.csv")

    print("\nCanonical Amino Acids")
    print("----------------------------------------------------------")

    canonical_amino_acid_object_list = generate_amino_acids(canonical_amino_acid_dataframe)

    for i in range(0, len(canonical_amino_acid_object_list)):
        subunits.update({canonical_amino_acid_object_list[i].multiple_letter: canonical_amino_acid_object_list[i]})

    print("\nNoncanonical Amino Acids")
    print("----------------------------------------------------------")

    noncanonical_amino_acid_object_list = generate_amino_acids(noncanonical_amino_acid_dataframe)

    for i in range(0, len(noncanonical_amino_acid_object_list)):
        subunits.update(
            {noncanonical_amino_acid_object_list[i].multiple_letter: noncanonical_amino_acid_object_list[i]})

    print("\nPeptoids")
    print("----------------------------------------------------------")

    peptoid_object_list = generate_peptoids(amine_dataframe)

    for i in range(0, len(peptoid_object_list)):
        subunits.update({peptoid_object_list[i].multiple_letter: peptoid_object_list[i]})

    print("\nMisc")
    print("----------------------------------------------------------")

    miscellaneous_object_list = generate_miscellaneous(miscellaneous_dataframe)

    for i in range(0, len(miscellaneous_object_list)):
        subunits.update({miscellaneous_object_list[i].multiple_letter: miscellaneous_object_list[i]})

    print("\nNumber of Valid subunits Generated = " + str(len(subunits)))

    return subunits
