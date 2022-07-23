"""
This script contains the classes that are used in the
various stages of populating the subunit library and creating
the peptides.
"""


class AminoAcid:
    """The amino acid class."""

    try:

        def __init__(self, name, multiple_letter, amino_group, side_chain,
                     stereocenter, alpha_carbon, carboxyl_group):

            self.name = name
            self.multiple_letter = multiple_letter
            self.amino_group = amino_group
            self.side_chain = side_chain
            self.stereocenter = stereocenter
            self.alpha_carbon = alpha_carbon
            self.carboxyl_group = carboxyl_group
            self.smiles_string = (amino_group + stereocenter +
                                  side_chain + alpha_carbon + carboxyl_group).rstrip()

    except:

        print("There was an error initializing an amino acid subunit.")





class Peptoid:
    """The peptoid class."""

    try:

        def __init__(self, name, multiple_letter, smiles_string):

            self.name = name
            self.multiple_letter = multiple_letter
            self.smiles_string = smiles_string.rstrip()

    except:

        print("There was an error initializing a peptoid subunit.")


class Miscellaneous:
    """The miscellaneous class. This class is used for
    subunits that do not fit into the other categories and do not
    need to follow special rules for their creation."""

    try:

        def __init__(self, name, multiple_letter, smiles_string):

            self.name = name
            self.multiple_letter = multiple_letter
            self.smiles_string = smiles_string.rstrip()

    except:

        print("There was an error initializing a miscellaneous subunit.")


class Peptide:
    """The peptide class."""

    try:

        def __init__(self, name, exact_mass, tpsa, a_log_p, smiles_string):

            self.name = name
            self.exact_mass = exact_mass
            self.tpsa = tpsa
            self.a_log_p = a_log_p
            self.smiles_string = smiles_string

    except:

        print("There was an error initializing a peptide.")
