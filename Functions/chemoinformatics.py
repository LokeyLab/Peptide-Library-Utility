# External
from rdkit.Chem import Descriptors, AllChem, AddHs
from rdkit.Chem import MolFromSmiles

def GetMoleculeFromSmiles(smilesString):
    molecule = MolFromSmiles(smilesString)

    return molecule

def GetExactMass(molecule):
    exactMass = Descriptors.ExactMolWt(molecule)

    return exactMass

def GetTPSA(molecule):
    TPSA = Descriptors.TPSA(molecule)

    return TPSA

def GetALogP(molecule):
    aLogP = Descriptors.MolLogP(molecule)

    return aLogP

def GetPredictedPapp(molecule):
    predictedPapp = "N/A"

    return predictedPapp

def GetChemometrics(molecule):

    moleculeH = AddHs(molecule)

    exactMass = GetExactMass(moleculeH)
    TPSA = GetTPSA(moleculeH)
    aLogP = GetALogP(moleculeH)
    predictedPapp = GetPredictedPapp(moleculeH)

    return exactMass, TPSA, aLogP, predictedPapp

def GetChemometricsWithUFFOptimazation(molecule):
    moleculeH = AddHs(molecule)

    AllChem.EmbedMolecule(moleculeH)

    AllChem.UFFOptimizeMolecule(moleculeH)

    exactMass = GetExactMass(moleculeH)
    TPSA = GetTPSA(moleculeH)
    aLogP = GetALogP(moleculeH)
    predictedPapp = GetPredictedPapp(moleculeH)

    return exactMass, TPSA, aLogP, predictedPapp

def GetChemometricsWithMMFFOptimazation(molecule):
    moleculeH = AddHs(molecule)

    AllChem.EmbedMolecule(moleculeH)

    AllChem.MMFFOptimizeMolecule(moleculeH)

    exactMass = GetExactMass(moleculeH)
    TPSA = GetTPSA(moleculeH)
    aLogP = GetALogP(moleculeH)
    predictedPapp = GetPredictedPapp(moleculeH)

    return exactMass, TPSA, aLogP, predictedPapp