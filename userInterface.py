# Internal
import subunitBuilder as s
import peptideBuilder as pb
import utilities as u

# External
import csv
from rdkit.Chem import MolFromSmiles

def CyclicPeptides():

    cyclicPeptides = []

    print("\nInput Type:")
    print("1) CSV File (Place the file in this folder)")
    print("2) Terminal Input")
    print("")

    choice = input("Enter a number: ")

    print("")

    if choice == "1":

        subunits = []

        f = input("Enter File Name: ") + ".csv"

        print("")

        reader = csv.reader(open(f))

        print("File Contents: ")

        for row in reader:
            print(row)
            subunits.append(row)

        possibleSubunits = s.GenerateSubunits()

        for subunitList in subunits:
            cyclicPeptides.append(pb.GenerateCyclicPeptide(subunitList, possibleSubunits))

        return cyclicPeptides

    elif choice == "2":

        cyclicPeptideCount = int(input("How many cyclic peptides?: "))

        for i in range(0, cyclicPeptideCount):
            subunits = []

            subunitCount = int(input("How many subunits?: "))
            print("")

            for subunit in range(0, subunitCount):

                subunits.append(input("Enter subunit (Ex: 'Leu' or 'D-B-H-Ala'): "))

            possibleSubunits = s.GenerateSubunits()

            cyclicPeptides.append(pb.GenerateCyclicPeptide(subunits[0], possibleSubunits))

        return cyclicPeptides

    else:
        MainUILoop()

def CombinatorialCyclicPeptides():

    return

def Introduction():
    print("1) Generate Cyclic Peptide(s)?")
    print("2) Generate Cyclic Peptide Library Combinatorially?")
    print("3) Calculate AlogP of a Cyclic Peptide SMILE String?")
    print("4) Calculate AlogP of a Cyclic Peptide SMILE String Library?")
    print("5) Check a Subunit?")
    print("6) Quit...\n")

    choice = input("Enter a number: ")

    return choice

def MainUILoop():

    choice = Introduction()

    if choice == "1":
        output = CyclicPeptides()

    elif choice == "2":
        output = CombinatorialCyclicPeptides()

    elif choice == "3":
        True

    elif choice == "4":
        True

    elif choice == "5":
        subunit = input("Enter subunit (Ex: 'Leu' or 'D-B-H-Ala'): ")

        possibleSubunits = s.GenerateSubunits()

        subunit = possibleSubunits[subunit]

        print(subunit.smilesString)

        MainUILoop()

    elif choice == "6":
        exit()

    else:
        print("Invalid input")

        MainUILoop()

    return output
