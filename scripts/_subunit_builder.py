# Internal
from scripts import _utilities as u

# External
import pandas as pd


class AminoAcid:
    def __init__(self, name, multipleLetter, aminoGroup, sideChain, stereocenter, carboxylGroup):
        self.name = name
        self.multipleLetter = multipleLetter
        self.aminoGroup = aminoGroup
        self.sideChain = sideChain
        self.stereocenter = stereocenter
        self.carboxylGroup = carboxylGroup
        self.smilesString = (aminoGroup + stereocenter + sideChain + "C" + carboxylGroup).rstrip()

class Peptoid:
    def __init__(self, name, multipleLetter, smilesString):
        self.name = name
        self.multipleLetter = multipleLetter
        self.smilesString = (smilesString).rstrip()

class Misc:
    def __init__(self, name, multipleLetter, smilesString):
        self.name = name
        self.multipleLetter = multipleLetter
        self.smilesString = (smilesString).rstrip()


def GenerateBaseAminoAcid(parentheses, name, multipleLetter, aminoGroup, sideChain, carboxylGroup):
    if parentheses == "y":
        baseAminoAcid = AminoAcid(name, multipleLetter, aminoGroup, "(" + sideChain + ")", "C",
                               carboxylGroup)

    elif parentheses == "n":
        baseAminoAcid = AminoAcid(name, multipleLetter, aminoGroup, sideChain, "C",
                               carboxylGroup)

    else:
        print("Parentheses not specified in spread sheet!")

    print("\n" + baseAminoAcid.name.capitalize()  + "s")
    print("----------------------------------------------------------")

    return baseAminoAcid


def KeepBaseVariant(baseAminoAcid):
    baseVariant = baseAminoAcid

    print(baseVariant.name + "     " + baseVariant.multiple_letter + "     " + baseVariant.smiles_string)

    return baseVariant


def GenerateLVariants(baseAminoAcid):
    lVariant = AminoAcid(baseAminoAcid.name, baseAminoAcid.multiple_letter, baseAminoAcid.aminoGroup,
                         baseAminoAcid.sideChain, "[C@@H]", baseAminoAcid.carboxylGroup)

    print(lVariant.name + "     " + lVariant.multipleLetter + "     " + lVariant.smilesString)

    return lVariant


def GenerateDVariants(baseAminoAcid):
    dVariant = AminoAcid("D-" + baseAminoAcid.name, "D-" + baseAminoAcid.multiple_letter, baseAminoAcid.aminoGroup,
                         baseAminoAcid.sideChain, "[C@H]", baseAminoAcid.carboxylGroup)

    print(dVariant.name.capitalize() + "     " + dVariant.multipleLetter + "     " + dVariant.smilesString)

    return dVariant


"""
def GenerateBetaVariants(baseAminoAcid):
    sideChain = baseAminoAcid.sideChain
    betaSideChain = baseAminoAcid.sideChain[0] + baseAminoAcid.sideChain[2:]

    betaVariant = AminoAcid("ß-" + baseAminoAcid.name, "ß-" + baseAminoAcid.multipleLetter, baseAminoAcid.aminoGroup,
                            betaSideChain, baseAminoAcid.stereocenter, "C" + baseAminoAcid.carboxylGroup)

    print(betaVariant.name + "      " + betaVariant.smilesString)

    return betaVariant


def GenerateHomoVariants(baseAminoAcid):
    sideChain = baseAminoAcid.sideChain
    homoSideChain = baseAminoAcid.sideChain[0] + "C" + baseAminoAcid.sideChain[1:]

    homoVariant = AminoAcid("homo-" + baseAminoAcid.name, "H-" + baseAminoAcid.multipleLetter, baseAminoAcid.aminoGroup,
                         homoSideChain, baseAminoAcid.stereocenter, baseAminoAcid.carboxylGroup)

    print(homoVariant.name + "      " + homoVariant.smilesString)

    return homoVariant

"""


def GenerateBetaHomoVariants(baseAminoAcid):
    if baseAminoAcid.name.startswith("D-"):
        betaHomoVariant = AminoAcid("ß-homo-" + baseAminoAcid.name, "ß-H-" + baseAminoAcid.multiple_letter,
                                    baseAminoAcid.aminoGroup,
                                    baseAminoAcid.sideChain, "[C@H]",
                                    "C" + baseAminoAcid.carboxylGroup)

    else:
        betaHomoVariant = AminoAcid("ß-homo-" + baseAminoAcid.name, "ß-H-" + baseAminoAcid.multiple_letter,
                                    baseAminoAcid.aminoGroup,
                                    baseAminoAcid.sideChain, "[C@@H]",
                                    "C" + baseAminoAcid.carboxylGroup)

    print(betaHomoVariant.name + "      " + betaHomoVariant.multipleLetter + "     " + betaHomoVariant.smilesString)

    return betaHomoVariant


def GenerateNMeVariants(baseAminoAcid):
    nMeVariant = AminoAcid("NMe-" + baseAminoAcid.name, "(NMe)-" + baseAminoAcid.multiple_letter, baseAminoAcid.aminoGroup + "(C)",
                           baseAminoAcid.sideChain, baseAminoAcid.stereocenter, baseAminoAcid.carboxylGroup)

    print(nMeVariant.name + "       " + nMeVariant.multipleLetter + "     " + nMeVariant.smilesString)

    return nMeVariant


def GenerateAminoAcids(df):

    aAList = []

    for i in range(0, len(df)):
        tempAAList = []

        data = df.loc[i]

        name = data.loc["Name"]
        multipleLetter = data.loc["Multiple Letter"]
        aminoGroup = data.loc["Amino Group"]
        sideChain = data.loc["Side Chain"]
        carboxylGroup = data.loc["Carboxyl Group"]

        parentheses = data.loc["Side Chain Parentheses? y/n"]
        lAndD = data.loc["L- and D-? y/n"]
        betaHomo = data.loc["ß-Homo-? y/n"]
        nMe = data.loc["NMe? y/n"]

        if multipleLetter == "Pro":
            print("\nProlines")
            print("----------------------------------------------------------")

            aminoAcid = AminoAcid(name, multipleLetter, aminoGroup, sideChain, "",  carboxylGroup)
            tempAAList += [aminoAcid]
            print(aminoAcid.name + "       " + aminoAcid.multipleLetter + "     " + aminoAcid.smilesString)

            aAList += tempAAList

            tempAAList = []

        elif multipleLetter == "ß-H-Pro":
            aminoAcid = AminoAcid(name, multipleLetter, aminoGroup, sideChain, "",
                                  "C" + carboxylGroup)
            tempAAList += [aminoAcid]
            print(aminoAcid.name + "       " + aminoAcid.multipleLetter + "     " + aminoAcid.smilesString)

            aAList += tempAAList

            tempAAList = []

        elif multipleLetter == "D-Pro":
            aminoAcid = AminoAcid(name, multipleLetter, aminoGroup, sideChain, "", carboxylGroup)
            tempAAList += [aminoAcid]
            print(aminoAcid.name + "       " + aminoAcid.multipleLetter + "     " + aminoAcid.smilesString)

            aAList += tempAAList

            tempAAList = []

        elif multipleLetter == "ß-H-D-Pro":
            aminoAcid = AminoAcid(name, multipleLetter, aminoGroup, sideChain, "",
                                  "C" + carboxylGroup)
            tempAAList += [aminoAcid]
            print(aminoAcid.name + "       " + aminoAcid.multipleLetter + "     " + aminoAcid.smilesString)

            aAList += tempAAList

            tempAAList = []

        else:
            baseAminoAcid = GenerateBaseAminoAcid(parentheses, name, multipleLetter, aminoGroup,
                                                  sideChain, carboxylGroup)

            if lAndD == "y":
                tempAAList.append(GenerateLVariants(baseAminoAcid))
                tempAAList.append(GenerateDVariants(baseAminoAcid))

            else:
                tempAAList.append(KeepBaseVariant(baseAminoAcid))

            if betaHomo == "y":
                betaHomoAAList = []

                for aminoAcid in tempAAList:
                    betaHomoAAList.append(GenerateBetaHomoVariants(aminoAcid))

                tempAAList += betaHomoAAList

            if nMe == "y":
                nMeAAList = []

                for aminoAcid in tempAAList:
                    nMeAAList.append(GenerateNMeVariants(aminoAcid))

                tempAAList += nMeAAList

            aAList += tempAAList

    print("----------------------------------------------------------")

    return aAList


def GeneratePeptoids(df):
    pList = []

    for i in range(0, len(df)):
        data = df.loc[i]

        name = str(data.loc["Name"])
        multipleLetter = str(data.loc["Multiple Letter"])
        smilesString = str(data.loc["SMILES String"])
        bromoaceticAcid = "BrCC(=O)O"

        if pd.isnull(data.loc["Multiple Letter"]):
            peptoid = Peptoid(name, "       ",  "N(" + smilesString[1:] + ")" + bromoaceticAcid[2:])

        else:
            peptoid = Peptoid(name, multipleLetter, "N(" + smilesString[1:] + ")" + bromoaceticAcid[2:])

        print("Bromoacetic Acid + " + peptoid.name + "     " + peptoid.multipleLetter + "     " + peptoid.smilesString)

        pList += [peptoid]

    print("----------------------------------------------------------")

    return pList


def GenerateMiscs(df):
    mList = []

    for i in range(0, len(df)):
        data = df.loc[i]

        name = str(data.loc["Name"])
        multipleLetter = str(data.loc["Multiple Letter"])
        smilesString = str(data.loc["SMILES String"])

        if pd.isnull(data.loc["Multiple Letter"]):
            misc = Misc(name, "       ", smilesString)

        else:
            misc = Misc(name, multipleLetter, smilesString)

        print(misc.name + "     " + misc.multipleLetter + "     " + misc.smilesString)

        mList += [misc]

    print("----------------------------------------------------------")
    return mList


def GenerateSubunitLibrary():

    subunits = {}

    canonicalAminoAcidDataframe = u.csv_to_dataframe("Subunits/Canonical Amino Acids-Canonical Amino Acids.csv")
    noncanonicalAminoAcidDataframe = u.csv_to_dataframe("Subunits/Non-Canonical Amino Acids-Non-Canonical Amino Acids.csv")
    amineDataframe = u.csv_to_dataframe("Subunits/Amines-Amines.csv")
    miscDataframe = u.csv_to_dataframe("Subunits/Misc-Misc.csv")

    print("\nCanonical Amino Acids")
    print("----------------------------------------------------------")

    canonicalAminoAcidObjectList = GenerateAminoAcids(canonicalAminoAcidDataframe)

    for i in range(0,len(canonicalAminoAcidObjectList)):
        subunits.update({canonicalAminoAcidObjectList[i].multiple_letter : canonicalAminoAcidObjectList[i]})

    print("\nNoncanonical Amino Acids")
    print("----------------------------------------------------------")

    noncanonicalAminoAcidObjectList = GenerateAminoAcids(noncanonicalAminoAcidDataframe)

    for i in range(0, len(noncanonicalAminoAcidObjectList)):
        subunits.update({noncanonicalAminoAcidObjectList[i].multiple_letter : noncanonicalAminoAcidObjectList[i]})

    print("\nPeptoids")
    print("----------------------------------------------------------")

    peptoidObjectList = GeneratePeptoids(amineDataframe)

    for i in range(0, len(peptoidObjectList)):
        subunits.update({peptoidObjectList[i].multiple_letter : peptoidObjectList[i]})

    print("\nMisc")
    print("----------------------------------------------------------")

    miscObjectList = GenerateMiscs(miscDataframe)

    for i in range(0, len(miscObjectList)):
        subunits.update({miscObjectList[i].multiple_letter : miscObjectList[i]})

    print("\nNumber of Valid Subunits Generated = " + str(len(subunits)))

    return subunits
