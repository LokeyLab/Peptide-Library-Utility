"""
This script contains the classes that are used in the
various stages of populating the subunit library and creating
the peptides.
"""


class AminoAcid:
    """The amino acid class."""

    try:

        def __init__(self, name, multiple_letter, amino_group, side_chain,
                     stereocenter, c_terminus, carboxyl_group):

            self.name = name
            self.multiple_letter = multiple_letter
            self.amino_group = amino_group
            self.side_chain = side_chain
            self.stereocenter = stereocenter
            self.c_terminus = c_terminus
            self.carboxyl_group = carboxyl_group
            self.smiles_string = (amino_group + stereocenter +
                                  side_chain + c_terminus + carboxyl_group).rstrip()

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
    # "Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "McGowan Volume", "# Atoms", "# Heteroatoms",
    # "# Amide Bonds",
    #                   "# Rings", "# Aromatic Rings", "# Rotatable Bonds", "# HBA", "# HBD"

    def __init__(self, name, smiles_string, exact_mass=None,
                 a_log_p=None, s_log_p=None, sasa=None, tpsa=None, mcgowan_volume=None,
                 num_atoms=None, num_heteroatoms=None, num_amide_bonds=None, num_rings=None,
                 num_aromatic_rings=None, num_rotatable_bonds=None, num_hba=None,
                 num_hbd=None, largest_ring_size=None):

        self.name = name
        self.smiles_string = smiles_string
        self.exact_mass = exact_mass
        self.a_log_p = a_log_p
        self.s_log_p = s_log_p
        self.sasa = sasa
        self.tpsa = tpsa
        self.mcgowan_volume = mcgowan_volume
        self.num_atoms = num_atoms
        self.num_heteroatoms = num_heteroatoms
        self.num_amide_bonds = num_amide_bonds
        self.num_rings = num_rings
        self.num_aromatic_rings = num_aromatic_rings
        self.num_rotatable_bonds = num_rotatable_bonds
        self.num_hba = num_hba
        self.num_hbd = num_hbd
        self.largest_ring_size = largest_ring_size
