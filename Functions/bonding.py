def BondHeadToTail(subunitsToBond):
    tempPeptide = subunitsToBond[0].smilesString

    for a in range(0, len(subunitsToBond)):

        if a != (len(subunitsToBond) - 1) and a != len(subunitsToBond):
            tempPeptide = tempPeptide[0: - 1] + subunitsToBond[a + 1].smilesString

    return tempPeptide


def CyclizeHeadToTail(peptide):
    tempPeptide = []

    for i in range(0, len(peptide)):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-1]) + "%99"

    return tempPeptide


def CyclizeHuisgen(peptide):
    tempPeptide = []

    return tempPeptide


def CyclizeThioether(peptide, subunitLibrary):

    tempPeptide = []

    if peptide.endswith(subunitLibrary["Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-24]) + "O[C@H](C)[C@@H](N)C(=O)O%99"
    elif peptide.endswith(subunitLibrary["D-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-23]) + "O[C@H](C)[C@H](N)C(=O)O%99"
    elif peptide.endswith(subunitLibrary["ß-H-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-25]) + "O[C@H](C)[C@@H](N)CC(=O)O%99"
    elif peptide.endswith(subunitLibrary["ß-H-D-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-24]) + "O[C@H](C)[C@H](N)CC(=O)O%99"
    elif peptide.endswith(subunitLibrary["(NMe)-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-27]) + "O[C@H](C)[C@@H](N(C))C(=O)O%99"
    elif peptide.endswith(subunitLibrary["(NMe)-D-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-26]) + "O[C@H](C)[C@H](N(C))C(=O)O%99"
    elif peptide.endswith(subunitLibrary["(NMe)-ß-H-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-28]) + "O[C@H](C)[C@@H](N(C))CC(=O)O%99"
    elif peptide.endswith(subunitLibrary["(NMe)-ß-H-D-Thr"].smilesString):
        tempPeptide = str(peptide[0]) + "%99" + str(peptide[1:-27]) + "O[C@H](C)[C@H](N(C))CC(=O)O%99"

    return tempPeptide

def CyclizeSideChainToSideChain(peptides):
    tempPeptide = []

    return peptides

def CyclizeHeadToSizeChain(peptides):
    tempPeptide = []

    return peptides
