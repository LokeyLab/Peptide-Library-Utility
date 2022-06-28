# External


def BondHeadToTail(subunitsToBond):
    tempPeptide = subunitsToBond[0].smilesString

    for a in range(0, len(subunitsToBond)):

        if a != (len(subunitsToBond) - 1) | a != len(subunitsToBond):
            tempPeptide = tempPeptide[0: - 1] + subunitsToBond[a + 1].smilesString

    return tempPeptide


def CyclizeHeadToTail(peptide):
    tempPeptide = []

    for i in range(0, len(peptide)):
        tempPeptide = str(peptide[0]) + "9" + str(peptide[1:-1]) + "9"

    return tempPeptide


def CyclizeHuisgen(peptides):
    tempPeptide = []

    return peptides


def CyclizeThioether(peptides):
    tempPeptide = []

    return peptides

def CyclizeSideChainToSideChain(peptides):
    tempPeptide = []

    return peptides

def CyclizeHeadToSizeChain(peptides):
    tempPeptide = []

    return peptides