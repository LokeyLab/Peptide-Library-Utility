#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script calculates AlogP using the Ghose-Crippen 1986 method.
"""
import itertools

import matplotlib.pyplot as plt
import pandas as pd

from scripts import header as h



global fragment_values
fragment_values = {
                    # C in:
                    "CH3R": -0.6327,                # 1
                    "CH4": -0.6327,                 # |
                    "CH2R2": -0.3998,               # 2
                    "CHR3": -0.2793,                # 3
                    "CR4": 0.2202,                  # 4
                    "CH3X": -1.1461,                # 5
                    "CH2RX": -0.9481,               # 6
                    "CH2X2": 0.2394,                # 7
                    "CHR2X": -0.9463,               # 8
                    "CHRX2": 0.5822,                # 9
                    "CHX3": 0.7245,                 # 10
                    ":CR3X": -1.0777,               # 11
                    "CR2X2": 1.1220,                # 12
                    ":CRX3": 0.6278,                # 13
                    "CX4": 1.2558,                  # 14
                    "=CH2": -0.2633,                # 15
                    "=CHR": -0.0460,                # 16
                    "=CR2": 0.3496,                 # 17
                    "=CHX": -0.3053,                # 18
                    "=CRX": -0.4451,                # 19
                    "=CX2": -0.1915,                # 20
                    "#CH": 0.1785,                  # 21
                    "#CR": 0.1541,                  # 22
                    "R=C=R":0.1541,                 # |
                    "#CX": None,                    # 23
                    "R--CH--R": -0.0548,            # 24
                    "R--CR--R": 0.3345,             # 25
                    "R--CX--R": -0.1153,            # 26
                    "R--CH--X": 0.0219,             # 27
                    "R--CR--X": 0.2093,             # 28
                    "R--CX--X": -0.1378,            # 29
                    "X--CH--X": -0.2686,            # 30
                    "X--CR--X":0.7376,              # 31
                    "X--CX--X": 0.0339,             # 32
                    "R--CH...X": 0.0230,            # 33
                    "R--CR...X": 0.2455,            # 34
                    "R--CX...X": -0.1883,           # 35
                    "Al-CH=X": 0.7853,              # 36
                    "Ar-CH=X": 0.1682,              # 37
                    "Al-C(=X)-Al": -0.4349,         # 38
                    "Ar-C(=X)-R": -0.2392,          # 39
                    "R-C(=X)-X": -0.1703,           # 40
                    "R-C#X": -0.1703,               # |
                    "X=C=X": -0.1703,               # |
                    "X-C(=X)-X": 0.0340,            # 41
                    "X--CH...X": -0.7231,           # 42
                    "X--CR...X": 0.2256,            # 43
                    "X--CX...X": -0.2692,           # 44
                                                    # 45 ununsed

                    # H attached to:
                    "H-C^0_sp3": 0.4307,              # 46
                    "H-C^1_sp3": 0.3722,              # 47
                    "H-C^0_sp2": 0.3722,              # |
                    "H-C^2_sp3": 0.0065,              # 48
                    "H-C^1_sp2": 0.0065,              # |
                    "H-C^0_sp" : 0.0065,               # |
                    "H-C^3_sp3": 0.4307,              # 49
                    "H-C^2_sp2": 0.4307,              # |
                    "H-C^3_sp2": 0.4307,              # |
                    "H-C^3_sp" : 0.4307,               # |
                    "H-Hereroatom": -0.3703,          # 50
                    "H-⍺-C": 0.2421,                  # 51
                                                    # 52 unused
                                                    # 53 unused
                                                    # 54 unused
                                                    # 55 unused

                    # O in:
                    "Alcohol": -0.0517,             # 56
                    "Phenol": 0.5212,               # 57
                    "Enol": 0.5212,                 # |
                    "Carboxyl": 0.5212,             # |
                    "=O": -0.1729,                  # 58
                    "Al-O-Al": 0.0407,              # 59
                    "Al-O-Ar": 0.3410,              # 60
                    "Ar_2O": 0.3410,                # |
                    "R...O...R": 0.3410,            # |
                    "R-O-C=X": 0.3410,              # |
                    "--O": 1.8020,                  # 61
                                                    # 62 unused
                                                    # 63 unused
                                                    # 64 unused
                                                    # 65 unused

                    # N in:
                    "Al-NH_2": 0.2658,              # 66
                    "Al_2NH": 0.2817,               # 67
                    "Al_3N": 0.3990,                # 68
                    "Ar-NH_2": 0.4442,              # 69
                    "X-NH_2": 0.4442,               # |
                    "Ar-NH-Al": 1.0841,             # 70
                    "Ar-NAl_2": 0.6632,             # 71
                    "RCO-N<": 0.1414,               # 72
                    ">N-X=X": 0.1414,               # |
                    "Ar_2NH": 0.3493,               # 73
                    "Ar_3N": 0.3493,                # |
                    "Ar_2N-Al": 0.3493,             # |
                    "R...N...R": 0.3493,            # |
                    "R#N": -0.1201,                 # 74
                    "R--N--R": 0.1757,              # 75
                    "R--N--X": 0.1757,              # |
                    "Ar-NO_2": -3.1516,             # 76
                    "R--N(--R)--O": -3.1516,        # |
                    "Ar-NO_2": -3.1516,             # |
                    "Al-NO_2": -3.3332,             # 77
                    "Ar-N=X": 0.1709,               # 78
                    "X-N=X": 0.1709,                # |
                                                    # 79 unused
                                                    # 80 unused

                    # F attached to:
                    "F-C^1_sp3": 0.4649,            # 81
                    "F-C^2_sp3": -0.1701,           # 82
                    "F-C^3_sp3": 0.1172,            # 83
                    "F-C^1_sp2": 0.6035,            # 84
                    "F-C^2_sp2": 0.4752,            # 85
                    "F-C^3_sp2": 0.4752,            # |
                    "F-C^4_sp2": 0.4752,            # |
                    "F-C^1_sp": 0.4752,             # |
                    "F-C^4_sp": 0.4752,             # |
                    "F-X": 0.4752,                  # |

                    # Cl attached to:
                    "Cl-C^1_sp3": 1.0723,            # 86
                    "Cl-C^2_sp3": 0.3027,           # 87
                    "Cl-C^3_sp3": 0.4108,            # 88
                    "Cl-C^1_sp2": 1.0278,            # 88
                    "Cl-C^2_sp2": 0.6972,            # 90
                    "Cl-C^3_sp2": 0.6972,            # |
                    "Cl-C^4_sp2": 0.6972,            # |
                    "Cl-C^1_sp": 0.6972,             # |
                    "Cl-C^4_sp": 0.6972,             # |
                    "Cl-X": 0.6972,                  # |

                    # Br attached to:
                    "Br-C^1_sp3": 1.0966,            # 91
                    "Br-C^2_sp3": 0.4292,           # 92
                    "Br-C^3_sp3": None,            # 93
                    "Br-C^1_sp2": 1.3224,            # 94
                    "Br-C^2_sp2": 0.9987,            # 95
                    "Br-C^3_sp2": 0.9987,            # |
                    "Br-C^4_sp2": 0.9987,            # |
                    "Br-C^1_sp": 0.9987,             # |
                    "Br-C^4_sp": 0.9987,             # |
                    "Br-X": 0.9987,                  # |

                    # I attached to:
                    "I-C^1_sp3": 1.4334,            # 96
                    "I-C^2_sp3": None,           # 97
                    "I-C^3_sp3": None,            # 98
                    "I-C^1_sp2": 1.8282,            # 99
                    "I-C^2_sp2": 1.0735,            # 100
                    "I-C^3_sp2": 1.0735,            # |
                    "I-C^4_sp2": 1.0735,            # |
                    "I-C^1_sp": 1.0735,             # |
                    "I-C^4_sp": 1.0735,             # |
                    "I-X": 1.0735,                  # |
                                                    # 101 unused
                                                    # 102 unused
                                                    # 103 unused
                                                    # 104 unused
                                                    # 105 unused

                    # S in:
                    "R-SH": 1.0152,                 # 106
                    "R_2S": 1.0339,                 # 107
                    "RS-SR": 1.0339,                 # |
                    "R=S": 0.0727,                 # 108
                    "R-SO-R": -0.3332,                 # 109
                    "R-O2-R": -0.1005,                 # 110
}

global smiles_fragments
smiles_fragments = {
                    # C in:
                    "C(-R)(-H)(-H)(-H)": fragment_values["CH3R"],           # 1
                    "C(-H)(-H)(-H)(-H)": fragment_values["CH4"],            # |
                    "C(-R)(-R)(-H)(-H)": fragment_values["CH2R2"],          # 2
                    "C(-R)(-R)(-R)(-H)": fragment_values["CHR3"],           # 3
                    "C(-R)(-R)(-R)(-R)": fragment_values["CR4"],            # 4
                    "C(-X)(-H)(-H)(-H)": fragment_values["CH3X"],           # 5
                    "C(-X)(-R)(-H)(-H)": fragment_values["CH3X"],           # 6
                    "C(-X)(-X)(-H)(-H)": fragment_values["CH2X2"],          # 7
                    "C(-X)(-R)(-R)(-H)": fragment_values["CHR2X"],          # 8
                    "C(-X)(-X)(-R)(-H)": fragment_values["CHRX2"],          # 9
                    "C(-X)(-X)(-X)(-H)": fragment_values["CHX3"],           # 10
                    "C(-X)(-R)(-R)(-R)": fragment_values[":CR3X"],           # 11
                    "C(-X)(-X)(-R)(-R)": fragment_values["CR2X2"],          # 12
                    "C(-X)(-X)(-X)(-R)": fragment_values[":CRX3"],           # 13
                    "C(-X)(-X)(-X)(-X)": fragment_values["CX4"],            # 14

                    "C(=R)(-H)(-H)": fragment_values["=CH2"],               # 15
                    "C(=X)(-H)(-H)": fragment_values["=CH2"],               # |
                    "C(=R)(-R)(-H)": fragment_values["=CHR"],               # 16
                    "C(=X)(-R)(-H)": fragment_values["=CHR"],               # |
                    "C(=R)(-R)(-R)": fragment_values["=CHR"],               # 17
                    "C(=X)(-R)(-R)": fragment_values["=CHR"],               # |
                    "C(=R)(-X)(-H)": fragment_values["=CHX"],               # 18
                    "C(=X)(-X)(-H)": fragment_values["=CHX"],               # |
                    "C(=R)(-X)(-R)": fragment_values["=CRX"],               # 19
                    # "C(=X)(-X)(-R)": fragment_values["=CRX"],             # |
                    "C(=R)(-X)(-X)": fragment_values["=CX2"],               # 20
                    # "C(=X)(-X)(-X)": fragment_values["=CX2"],             # |

                    "C(#X)(-H)": fragment_values["#CH"],                    # 21
                    "C(#R)(-H)": fragment_values["#CH"],                    # |
                    "C(#R)(-R)": fragment_values["#CR"],                    # 22
                    # "C(#X)(-R)": fragment_values["#CR"],  # |
                    "C(=R)(=R)": fragment_values["R=C=R"],                  # |
                    "C(#X)(-X)": fragment_values["#CX"],                    # 23
                    "C(#R)(-X)": fragment_values["#CX"],                    # |

                    "C(*R)(*R)(-H)": fragment_values["R--CH--R"],           # 24
                    "C(*R)(*R)(-R)": fragment_values["R--CR--R"],           # 25
                    "C(-X)(*R)(*R)": fragment_values["R--CX--R"],           # 26
                    "C(*X)(*R)(*R)": fragment_values["R--CX--R"],           # |
                    "C(*X)(*R)(-H)": fragment_values["R--CH--X"],           # 27
                    "C(*X)(*R)(-R)": fragment_values["R--CR--X"],           # 28
                    "C(*X)(-X)(*R)": fragment_values["R--CX--X"],           # 29
                    "C(*X)(*X)(-H)": fragment_values["X--CH--X"],           # 30
                    "C(*X)(*X)(-R)": fragment_values["X--CR--X"],           # 31
                    "C(*X)(*X)(*R)": fragment_values["X--CR--X"],           # |

                    "C(*X)(*X)(-X)": fragment_values["X--CX--X"],           # 32

                    # "C(*X)(*R)(-H)": fragment_values["R--CH...X"],          # 33
                    # "C(*X)(*R)(-R)": fragment_values["R--CR...X"],          # 34
                    # "C(*X)(-X)(*R)": fragment_values["R--CX...X"],          # 35

                    "C(=X)(-Al)(-H)": fragment_values["Al-CH=X"],          # 36
                    "C(=X)(-Ar)(-H)": fragment_values["Ar-CH=X"],          # 37
                    "C(=X)(-Al)(-Al)": fragment_values["Al-C(=X)-Al"],         # 38
                    "C(=X)(-Ar)(-R)": fragment_values["Ar-C(=X)-R"],        # 39
                    "C(=X)(-X)(-R)": fragment_values["R-C(=X)-X"],          # 40

                    "C(#X)(-R)": fragment_values["R-C#X"],                  # |
                    "C(=X)(=X)": fragment_values["X=C=X"],                  # |
                    "C(=X)(-X)(-X)": fragment_values["X-C(=X)-X"],          # 41

                    # "C(*X)(*X)(-H)": fragment_values["X--CH...X"],          # 42
                    # "C(*X)(*X)(-R)": fragment_values["X--CR...X"],          # 43
                    # "C(*X)(*X)(-X)": fragment_values["X--CX...X"],          # 44

                    # H attached to:
                    "H(-C^0_sp3)": fragment_values["H-C^0_sp3"],          # 46
                    "H(-C^1_sp3)": fragment_values["H-C^1_sp3"],          # 47
                    "H(-C^0_sp2)": fragment_values["H-C^0_sp2"],          # |
                    "H(-C^2_sp3)": fragment_values["H-C^2_sp3"],          # 48
                    "H(-C^1_sp2)": fragment_values["H-C^1_sp2"],          # |
                    "H(-C^0_sp)" : fragment_values["H-C^0_sp"],           # |
                    "H(-C^3_sp3)": fragment_values["H-C^3_sp3"],          # 49
                    "H(-C^2_sp2)": fragment_values["H-C^2_sp2"],          # |
                    "H(-C^3_sp2)": fragment_values["H-C^3_sp2"],          # |
                    "H(-C^3_sp)" : fragment_values["H-C^3_sp"],           # |

                    "H(-X)": fragment_values["H-Hereroatom"],             # 50
                    # "H(-⍺C)": fragment_values["H-⍺-C"],                     # 51 FIXME

                    # O in:
                    "O(-R)(-H)": fragment_values["Alcohol"],                    # 56
                    "O(-H)(-R)": fragment_values["Alcohol"],                    # |
                    "O(-H)(-Al)": fragment_values["Alcohol"],                    # |

                    # "O(-H)": fragment_values["Phenol"],                     # 57 FIXME
                    # "O(-H)": fragment_values["Enol"],                       # | FIXME
                   #  "O(-H)": fragment_values["Carboxyl"],                   # | FIXME

                    "O(=X)": fragment_values["=O"],                         # 58
                    "O(=R)": fragment_values["=O"],                         # |
                    "O(=Al)": fragment_values["=O"],                        # |

                    "O(-Al)(-Al)": fragment_values["Al-O-Al"],              # 59
                    "O(-Ar)(-Al)": fragment_values["Al-O-Ar"],              # 60
                    "O(-Ar)(-Ar)": fragment_values["Ar_2O"],                # |
                    "O(*R)(*R)": fragment_values["R...O...R"],              # |
                    "O(-R)(-C=X)": fragment_values["R-O-C=X"],              # |
                    "O(*X)": fragment_values["--O"],                         # 61
                    "O(*R)": fragment_values["--O"],                        # |

                    # N in:
                    "N(-Al)(-H)(-H)": fragment_values["Al-NH_2"],           # 66
                    "N(-Al)(-Al)(-H)": fragment_values["Al_2NH"],          # 67
                    "N(-Al)(-Al)(-Al)": fragment_values["Al_3N"],           # 68
                    "N(-Ar)(-H)(-H)": fragment_values["Ar-NH_2"],           # 69
                    "N(-X)(-H)(-H)": fragment_values["X-NH_2"],            # |
                    "N(-Ar)(-Al)(-H)": fragment_values["Ar-NH-Al"],         # 70
                    "N(-Ar)(-Al)(-Al)": fragment_values["Ar-NAl_2"],         # 71
                    # "N(-Ar)(-Al)(-Al)": fragment_values["RCO-N<"],          # 72 FIXME
                    # "N(-Ar)(-Al)(-Al)": fragment_values["RCO-N<"],          # | FIXME
                    "N(-Ar)(-Ar)(-H)": fragment_values["Ar_2NH"],           # 73
                    "N(*Ar)(*Ar)(-H)": fragment_values["Ar_2NH"],           # |
                    "N(-Ar)(-Ar)(-Ar)": fragment_values["Ar_3N"],           # |
                    "N(*Ar)(*Ar)(*Ar)": fragment_values["Ar_3N"],                # |
                    "N(-Ar)(-Ar)(-Al)": fragment_values["Ar_2N-Al"],         # |
                    "N(*Ar)(*Ar)": fragment_values["R...N...R"],                # |
                    "N(#R)": fragment_values["R#N"],                        # 74
                    "N(*R)(*R)": fragment_values["R--N--R"],                 # 75
                    "N(*X)(*X)": fragment_values["R--N--X"],                 # |
                    "N(-Ar)(-O)(-O)": fragment_values["Ar-NO_2"],            # 76 FIXME
                    "N(*R)(*R)(*O)": fragment_values["Ar-NO_2"],            # |
                    "N(-Al)(-O)(-O)": fragment_values["Al-NO_2"],           # 77
                    "N(-Ar)(=X)": fragment_values["Ar-N=X"],               # 78
                    "N(-Al)(=X)": fragment_values["Ar-N=X"],               # |
}

def get_fragment_contribution(fragment):

    contribution = 0

    atom = fragment[0]

    neighbors = []
    neighbor = ""

    for i in range(1, len(fragment)):

        if not fragment[i] == ")":
            neighbor += fragment[i]

        else:
            neighbor += fragment[i]
            neighbors.append(neighbor)
            neighbor = ""

    # print(neighbors)

    neighbor_variants = []

    for p in h.itertools.permutations(neighbors, len(neighbors)):

        neighbor_variants.append(atom + ''.join(p))

    # print(neighbor_variants)

    contribution_list = []

    for v in neighbor_variants:
        try:
            contribution = smiles_fragments[v]

        except KeyError:
            pass

    return contribution


def get_a_log_p_86(smiles):
    # smiles = "N41[C@H](CCC1)C(=O)N(CCC)CC(=O)N[C@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N2[C@@H](
    # CCC2)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O"
    # target_value = 2.54080000000001

    # smiles = "C1COCCN1N=O"
    # target_value = - 0.635

    # smiles = "C1=NC2=NC=NC(=C2N1)N"
    # target_value = -0.383

    m = h.MolFromSmiles(smiles)
    m = h.AddHs(m)

    num_atoms = h.rdMolDescriptors.CalcNumAtoms(m)

    atoms = []

    fragments = []

    for a, atom in enumerate(m.GetAtoms()):
        fragment = ""

        atoms.append(atom)

        # print(atom.GetSymbol(), end ="")

        fragment += atom.GetSymbol()

        for n, neighbor in enumerate(atom.GetNeighbors()):
            bond_string = ""
            bond = str(m.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType())

            if bond == "SINGLE":
                bond_string = "-"
            if bond == "DOUBLE":
                bond_string = "="
            if bond == "TRIPLE":
                bond_string = "#"
            if bond == "QUADRUPLE":
                bond_string = ":"
            if bond == "AROMATIC":
                bond_string = "*"

            # print("(" + bond_string + neighbor.GetSymbol() + ")", end ="")

            neighbor_symbol = neighbor.GetSymbol()

            if atom.GetSymbol() == "C":
                if neighbor_symbol == "C":
                    neighbor_symbol = "R"

                elif neighbor_symbol == "H":
                    neighbor_symbol = "H"

                else:
                    neighbor_symbol = "X"

                fragment += "(" + bond_string + neighbor_symbol + ")"

            if atom.GetSymbol() == "N":
                if neighbor_symbol == "C":
                    if neighbor.GetIsAromatic():
                        neighbor_symbol = "Ar"
                    else:
                        neighbor_symbol = "Al"

                elif neighbor_symbol == "H":
                    neighbor_symbol = "H"

                elif neighbor_symbol == "O":
                    neighbor_symbol = "X"

                elif neighbor_symbol == "N":
                    if neighbor.GetIsAromatic():
                        neighbor_symbol = "Ar"
                    else:
                        neighbor_symbol = "Al"

                fragment += "(" + bond_string + neighbor_symbol + ")"

            if atom.GetSymbol() == "O":
                if neighbor_symbol == "C":
                    if neighbor.GetIsAromatic():
                        neighbor_symbol = "Ar"
                    else:
                        neighbor_symbol = "Al"

                elif neighbor_symbol == "H":
                    neighbor_symbol = "H"

                else:
                    neighbor_symbol = "X"

                fragment += "(" + bond_string + neighbor_symbol + ")"

            if atom.GetSymbol() == "H" or \
               atom.GetSymbol() == "F" or \
               atom.GetSymbol() == "Cl" or \
               atom.GetSymbol() == "Br" or \
               atom.GetSymbol() == "I":

                if neighbor_symbol == "C":
                    hybridiation_state = neighbor.GetHybridization()
                    oxidation_number = 0

                    carbon_neighbors = neighbor.GetNeighbors()

                    # print()

                    # print("Neighbors: ")

                    for c_n, carbon_neighbor in enumerate(carbon_neighbors):

                        # print(carbon_neighbor.GetSymbol())

                        if carbon_neighbor.GetSymbol() == "H":
                            oxidation_number -= 1

                        if carbon_neighbor.GetSymbol() == "C":
                            oxidation_number += 0

                        if carbon_neighbor.GetSymbol() in ["O", "N", "F",
                                                                "S", "Cl",
                                                                     "Br",
                                                                     "I"]:
                            oxidation_number += 1

                    oxidation_number = abs(oxidation_number)

                    # print(oxidation_number)

                    # print(hybridiation_state)

                    fragment += f"({bond_string}{neighbor_symbol}^{oxidation_number}_" \
                                f"{str(hybridiation_state).lower()})"

                else:
                    neighbor_symbol = "X"

                    fragment += "(" + bond_string + neighbor_symbol + ")"


        fragments.append(fragment)

    a_log_p = 0

    contributions = []
    contribution_over_time = []

    for i, frag in enumerate(fragments):
        # print(frag)
        contribution = get_fragment_contribution(frag)
        contributions.append(contribution)

        """
        if contribution == 0:
            print(f"{frag} contribution = !!!! ERROR !!!!")

        else:
            print(f"{frag} contribution = {contribution}")
        # print()"""
        a_log_p += contribution

        contribution_over_time.append(a_log_p)

    """
    h.plt.plot(contributions)
    h.plt.plot(contribution_over_time)
    h.plt.show()
    """

    """print()
    print(f"Calculated AlogP = {a_log_p}"
          f"\nTarget AlogP = {target_value}"
          f"\nDifference= {target_value - a_log_p}")
    print()"""

    return a_log_p

def a_log_p_test():

    sample_smiles = pd.read_csv("/Users/Adam/Desktop/Programming/github-repos/Peptide-Library-Utility/scripts/a_log_p/naylor2018.csv")["SMILE"]

    sample_a_log_p = pd.read_csv("/Users/Adam/Desktop/Programming/github-repos/Peptide-Library-Utility/scripts/a_log_p/naylor2018.csv")["ALogP_86"]

    # print(sample_smiles, sample_a_log_p)
    test_data = []

    for i in h.tqdm.tqdm(range(len(sample_smiles))):
        test_data.append(-get_a_log_p_86(sample_smiles[i]))

    # print(test_data)

    slope, intercept, r_value, p_value, std_err = h.stats.linregress(sample_a_log_p, test_data)
    entropy = h.stats.entropy(sample_a_log_p, test_data)

    diff = test_data - sample_a_log_p
    mean_diff = h.np.mean(diff)

    fig, (ax1, ax2) = h.plt.subplots(2,1)

    ax1.plot(sample_a_log_p)
    ax1.plot(test_data)

    ax1.text(0.95, 0.95, s=f"r = {round(r_value, 4)}",
             transform=ax1.transAxes,
             verticalalignment='top', horizontalalignment="right")

    ax1.set_ylabel("AlogP")

    ax1.legend(["Original", "Test"], loc=2)

    ax2.plot(diff)

    # ax2.set_ylim([-1,1])

    ax2.text(0.95, 0.95, s=f"Mean Diff. = {round(mean_diff, 4)}", transform=ax2.transAxes,
             verticalalignment='top', horizontalalignment="right")

    ax2.set_ylabel("Difference")

    h.plt.tight_layout()

    h.plt.show()

# a_log_p_test()
