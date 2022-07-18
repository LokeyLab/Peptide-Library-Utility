# Internal
from scripts import bonding as bnd, cheminformatics as cmi, utilities as utl, classes as cls


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
    peptide = bnd.bond_head_to_tail(subunits)
    smilesString = bnd.cyclize_head_to_tail(peptide)
    molecule = cmi.get_molecule_from_smiles(smilesString)
    exactMass = utl.GetExactMass(molecule)
    logP = utl.GetALogP(molecule)
    cyclicPeptideObj = cls.Peptide(name, exactMass, logP, smilesString)

    return cyclicPeptideObj