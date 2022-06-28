# Internal
import bonding as b
import subunitBuilder as s
import userInterface as ui
import utilities as u

# External
from rdkit.Chem import MolFromSmiles


class Peptide:
    def __init__(self, name, exactMass, logP, smilesString):
        self.name = name
        self.exactMass = exactMass
        self.logP = logP
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
    molecule = MolFromSmiles(smilesString)
    exactMass = u.GetExactMass(molecule)
    logP = u.GetALogP(molecule)
    cyclicPeptideObj = Peptide(name, exactMass, logP, smilesString)

    return cyclicPeptideObj