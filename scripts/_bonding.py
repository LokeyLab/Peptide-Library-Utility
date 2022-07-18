def bond_head_to_tail(subunitsToBond):
    tempPeptide = subunitsToBond[0].smiles_string

    for a in range(0, len(subunitsToBond)):

        if a != (len(subunitsToBond) - 1) and a != len(subunitsToBond):
            tempPeptide = tempPeptide[0: - 1] + subunitsToBond[a + 1].smiles_string

    return tempPeptide


def cyclize_head_to_tail(peptide):
    tempPeptide = []

    for i in range(0, len(peptide)):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-1]) + "%99"

    return tempPeptide


def cyclize_huisgen(peptide):
    tempPeptide = []

    return tempPeptide


def cyclize_thioether(peptide, subunitLibrary):

    tempPeptide = []

    if peptide.endswith(subunitLibrary["Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-24]) + "O[C@H](C)[C@@H](N)C(=O)%99"

    elif peptide.endswith(subunitLibrary["D-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-23]) + "O[C@H](C)[C@H](N)C(=O)%99"

    elif peptide.endswith(subunitLibrary["ß-H-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-25]) + "O[C@H](C)[C@@H](N)CC(=O)%99"

    elif peptide.endswith(subunitLibrary["ß-H-D-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-24]) + "O[C@H](C)[C@H](N)CC(=O)%99"

    elif peptide.endswith(subunitLibrary["(NMe)-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-27]) + "O[C@H](C)[C@@H](N(C))C(=O)%99"

    elif peptide.endswith(subunitLibrary["(NMe)-D-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-26]) + "O[C@H](C)[C@H](N(C))C(=O)%99"

    elif peptide.endswith(subunitLibrary["(NMe)-ß-H-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-28]) + "O[C@H](C)[C@@H](N(C))CC(=O)%99"

    elif peptide.endswith(subunitLibrary["(NMe)-ß-H-D-Thr"].smiles_string):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-27]) + "O[C@H](C)[C@H](N(C))CC(=O)%99"

    else:
        print("\nError. Peptide may not end in a threonine.") # N[C@@H](Cc%98ccc(F)cc%98)C(=O)O

    return tempPeptide

def CyclizeSideChainToSideChain(peptides):
    tempPeptide = []

    return peptides

def CyclizeHeadToSizeChain(peptides):
    tempPeptide = []

    return peptides
