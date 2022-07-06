# Internal
from Functions import bonding as b, chemoinformatics as c, utilities as u

# External


class Peptide:
    def __init__(self, name, exactMass, TPSA, ALogP, predictedPapp, smilesString):
        self.name = name
        self.exactMass = exactMass
        self.TPSA = TPSA
        self.ALogP = ALogP
        self.predictedPapp = predictedPapp
        self.smilesString = smilesString

def GetSubunits(desiredSubunits, possibleSubunits):
    peptideSubunits = []

    for i in range(0, len(desiredSubunits)):
        subunit = desiredSubunits[i]

        peptideSubunits.append(possibleSubunits[subunit])

    return peptideSubunits

def GetNameOfCyclicPeptide(subunits):
    name = ""
    for subunit in subunits:
        name += subunit.name + " "

    return name

def GenerateCyclicPeptide(desiredSubunits, possibleSubunits):
    subunits = GetSubunits(desiredSubunits, possibleSubunits)
    name = GetNameOfCyclicPeptide(subunits)
    peptide = b.BondHeadToTail(subunits)
    smilesString = b.CyclizeHeadToTail(peptide)
    molecule = c.GetMoleculeFromSmiles(smilesString)
    exactMass = u.GetExactMass(molecule)
    logP = u.GetALogP(molecule)
    cyclicPeptideObj = Peptide(name, exactMass, logP, smilesString)

    return cyclicPeptideObj