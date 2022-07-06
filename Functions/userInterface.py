# Internal
from Functions import bonding as b, chemoinformatics as ch, peptideBuilder as pb, subunitBuilder2 as sb2, \
    utilities as u, combinatronics as c

# External
import numpy as np
import pandas as pd


def CyclicPeptideInTerminal(subunitLibrary):
    cyclicPeptides = []

    while True:
        subunitCount = input("\nHow many subunits?: ")
        if subunitCount.isnumeric():
            subunitCount = int(subunitCount)
            break

        else:
            print("\nInvalid input, try again.")
            continue

    subunits = []

    print("\nExample input: (NMe)-ß-H-D-Leu or Leu")

    for i in range(0, subunitCount):
        while True:
            try:
                subunit = subunitLibrary[input("\nEnter subunit: ")]
                break

            except:
                print("\nInvalid input, try again.")

                continue

        print("\nSubunit " + str(i + 1) + ": " + subunit.name)

        subunits.append(subunit)

    print("")

    name = ""

    for i in range(0, len(subunits)):
        name += subunits[i].multipleLetter + " "

    peptide = b.BondHeadToTail(subunits)

    print("Cyclization Types:")
    print("1) Head to Tail")
    print("2) Huisgen FIX ME")
    print("3) Thioether (Must end in Threonine)")

    choice = input("\nEnter a number: ")
    while True:
        if choice == "1":
            print("\nCyclizing Head to Tail...")

            cyclicPeptide = b.CyclizeHeadToTail(peptide)

            print("\nCyclization Complete.")

            break

        elif choice == "2":
            print("\nCyclizing Huisgen...")

            cyclicPeptide = b.CyclizeHuisgen(peptide)

            print("\nCyclization Complete.")

            break

        elif choice == "3":
            print("\nCyclizing Thioether...")

            cyclicPeptide = b.CyclizeThioether(peptide, subunitLibrary)

            print("\nCyclization Complete.")

            while True:
                choice = input("\nCap exposed Amino Group? y/n: ")

                if choice == "y":
                    tempPeptide = str(cyclicPeptide[0:cyclicPeptide.rfind("N") + 1]) + "(" + \
                                  subunitLibrary[input("\nEnter subunit: ")].smilesString + ")" + \
                                  str(cyclicPeptide[cyclicPeptide.rfind("N") + 1:])

                    cyclicPeptide = tempPeptide
                    break

                elif choice == "n":
                    break

                else:
                    print("\nInvalid input, try again.")

        else:
            print("\nInvalid input, try again.")

            continue

    molecule = ch.GetMoleculeFromSmiles(cyclicPeptide)

    exactMass, TPSA, aLogP, predictedPapp = ch.GetChemometrics(molecule)

    cyclicPeptideObject = pb.Peptide(name, exactMass, TPSA, aLogP, predictedPapp, cyclicPeptide)

    df = u.CyclicPeptidesToDataframe([cyclicPeptideObject])

    u.PrintDataframe(df)

    u.DataframeToCSV(df)

    return cyclicPeptide


def CyclicPeptideLibraryFromCSV(subunitLibrary):
    while True:
        try:
            df = u.ReadFileToDataframe("Input/" + input("\nEnter file name: ") + ".csv")

            break

        except:
            print("\nInvalid input, try again.")

            continue

    u.PrintDataframe(df)

    potSizes = []
    pots = []
    for i in range(0, df.shape[1]):
        tempPot = []
        count = 0

        for j in range(0, df.shape[0]):
            if pd.notnull(df.iloc[j, i]):
                tempPot.append(df.iloc[j, i])
                count += 1

        potSizes.append(count)
        pots.append(tempPot)

    cyclicPeptideCount = potSizes[0]

    for k in range(1, len(potSizes)):
        cyclicPeptideCount *= potSizes[k]

    print("\nCyclic Peptides To Be Generated: " + str(cyclicPeptideCount))

    cartesianProduct = c.CartesianProduct(pots)

    while True:
        print("\nCyclization Types:")
        print("1) Head to Tail")
        print("2) Huisgen FIX ME")
        print("3) Thioether (Must end in Threonine)")

        choice = input("\nEnter a number: ")

        if choice == "1":
            print("\nCyclizing Head to Tail...")

            cyclizationType = "Head to Tail"

            break

        elif choice == "2":
            print("\nCyclizing Huisgen...")

            cyclizationType = "Huisgen"

            break

        elif choice == "3":
            print("\nCyclizing Thioether...")

            cyclizationType = "Thioether"

            break

        else:
            print("\nInvalid input, try again.")

            continue

    capBool = False

    while True:
        choice = input("\nCap exposed Amino Group? y/n: ")

        if choice == "y":
            cap = input("\nEnter subunit: ")

            capBool = True
            break

        elif choice == "n":
            capBool = False
            break
        else:
            print("\nInvalid input, try again.")

    cyclicPeptides = []

    for l in cartesianProduct:
        subunits = []
        for m in l:
            subunits.append(subunitLibrary[m])
        peptide = b.BondHeadToTail(subunits)
        if cyclizationType == "Head to Tail":
            cyclicPeptide = b.CyclizeHeadToTail(peptide)
        if cyclizationType == "Huisgen":
            cyclicPeptide = b.CyclizeHuisgen(peptide)
        if cyclizationType == "Thioether":
            cyclicPeptide = b.CyclizeThioether(peptide, subunitLibrary)

            if capBool:
                subunits.append(subunitLibrary[cap])
                tempPeptide = str(cyclicPeptide[0:cyclicPeptide.rfind("N") + 1]) + "(" + \
                              subunitLibrary[cap].smilesString + ")" + \
                              str(cyclicPeptide[cyclicPeptide.rfind("N") + 1:])

                cyclicPeptide = tempPeptide

        name = ""

        for i in range(0, len(subunits)):
            name += subunits[i].multipleLetter + " "

        molecule = ch.GetMoleculeFromSmiles(cyclicPeptide)

        exactMass, TPSA, aLogP, predictedPapp = ch.GetChemometrics(molecule)

        cyclicPeptideObject = pb.Peptide(name, exactMass, TPSA, aLogP, predictedPapp, cyclicPeptide)

        cyclicPeptides.append(cyclicPeptideObject)

    print("\nCyclizations Complete.")

    df = u.CyclicPeptidesToDataframe(cyclicPeptides)

    u.PrintDataframe(df)

    u.DataframeToCSV(df)

    u.PlotDataframe(df)

    return


def Introduction():
    print("0) Populate Subunit Library")
    print("1) Create Cyclic Peptide in Terminal FIX ME (add Huisgen)")
    print("2) Create Cyclic Peptide Library from .CSV FIX ME (add Huisgen)")
    print("3) Calculate Chemometrics of a SMILES string")
    print("4) Check a Subunit")
    print("5) Quit\n")

    choice = input("Enter number: ")

    return choice


def MainUILoop(subunitLibrary):
    print("\n-------------------- Main Menu --------------------")
    choice = Introduction()

    if choice == "0":
        print("\nPopulating Subunit Library...\n")

        subunitLibrary = sb2.GenerateSubunitLibrary()

        print("\nSubunit Library Populated.")

    elif choice == "1":
        print("\nCreating Cyclic Peptide...\n")

        cyclicPeptide = CyclicPeptideInTerminal(subunitLibrary)

        print("\nCyclic Peptide Created.")

    elif choice == "2":
        print("\nCreating Cyclic Peptide Library...")

        cyclicPeptides = CyclicPeptideLibraryFromCSV(subunitLibrary)

        print("\nCyclic Peptide Library Created.")

    elif choice == "3":
        while True:
            try:
                cyclicPeptide = input("\nEnter SMILES String: ")

                molecule = ch.GetMoleculeFromSmiles(cyclicPeptide)

                exactMass, TPSA, aLogP, predictedPapp = ch.GetChemometrics(molecule)

                cyclicPeptideObject = pb.Peptide("", exactMass, TPSA, aLogP, predictedPapp, cyclicPeptide)

                df = u.CyclicPeptidesToDataframe([cyclicPeptideObject])

                u.PrintDataframe(df)

                break

            except:
                print("\nInvalid input, try again.")

    elif choice == "4":
        subunit = subunitLibrary[input("Enter subunit: ")]

        print("")

        molecule = ch.GetMoleculeFromSmiles(subunit.smilesString)

        exactMass, TPSA, aLogP, predictedPapp = ch.GetChemometrics(molecule)

        columns = ["Name", "Exact Mass", "TPSA", "ALogP", "Predicted PappE-6", "SMILES String"]

        df = pd.DataFrame(columns=columns)

        df.loc[len(df.index)] = {"Name": subunit.multipleLetter,
                                 "Exact Mass": exactMass,
                                 "TPSA": TPSA,
                                 "ALogP": aLogP,
                                 "Predicted PappE-6": 0,
                                 "SMILES String": subunit.smilesString}

        u.PrintDataframe(df)

    elif choice == "5":
        print("\nGoodbye!")

        exit()

    MainUILoop(subunitLibrary)
