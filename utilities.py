# External
from rdkit.Chem import Descriptors
import csv

def PrintPeptides(cyclicPeptides):
    print("Number of peptides = " + str(len(cyclicPeptides)) + "\n")

    print("Name -> Exact Mass -> LogP -> SMILES String")
    for cyclicPeptide in cyclicPeptides:
        print(cyclicPeptide.name + "   " + str(cyclicPeptide.exactMass) + "   " + str(cyclicPeptide.logP) + "   " + cyclicPeptide.smilesString)

    print("")

def PeptidesToCSV(cyclicPeptides):
    f = open('cyclicPeptides.csv', 'w')
    writer = csv.writer(f)


    writer.writerow(['Name', 'Exact Mass', 'LogP', 'SMILES String'])

    for cyclicPeptide in cyclicPeptides:
        writer.writerow([cyclicPeptide.name, str(cyclicPeptide.exactMass), str(cyclicPeptide.logP), cyclicPeptide.smilesString])

    f.close()

def GetALogP(molecule):
    aLogP = Descriptors.MolLogP(molecule)

    return aLogP

def GetExactMass(molecule):
    exactMass = Descriptors.ExactMolWt(molecule)

    return exactMass