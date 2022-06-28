# External
import time

# Generic Structures
aAmineGroup = "N"
carboxylGroup = "(=O)O"
aCarbon = "C"
bCarbon = "C"
aNull = ""
aL = "[C@@H]"
aD = "[C@H]"
bL = "[C@H]"
bD = "[C@@H]"

# Canonical Amino Acids

# Side Chains
alaSC = "(C)"
argSC = "(CCCN=C(N)N)"
asnSC = "(CC(=O)N)"
aspSC = "(CC(=O)O)"
cysSC = "(CS)"
hisSC = "(CC1=CN=C-N1)"
ileSC = "([C@@H](C)CC)"
leuSC = "(CC(C)C)"
lysSC = "(CCCCN)"
glnSC = "(CCC(=O)N)"
gluSC = "(C(CC(=O)O))"
glySC = "C"
metSC = "(CSCC)"
pheSC = "(Cc2ccccc2)"
proSC = "3CCC[C@H]3"
dproSC = "3CCC[C@@H]3"
thrSC = "([C@H](O)C)"
trpSC = "(Cc4=c[NH]c5c4cccc5)"
tyrSC = "(Cc1ccc(O)cc1)"
serSC = "(C(O))"
valSC = "(C(C)C)"

# Nc1ccc(cc1)C(=O)O

# Noncanonical Amino Acids

# Side Chains
chaSC = ""

class AlphaAminoAcid:

    def __init__(self, name, aAmineGroup, aL, sideChain, aCarbon, carboxylGroup):
        self.name = name
        self.aAmineGroup = aAmineGroup
        self.aCarbon = aCarbon
        self.sideChain = sideChain
        self.carboxylGroup = carboxylGroup
        self.smilesString = aAmineGroup + aL + sideChain + aCarbon + carboxylGroup

"""
class BetaAminoAcid(AlphaAminoAcid):
	def __init__():

		self.smilesString = aAmineGroup + aCarbon + sideChain + carboxylGroup


class HomoAminoAcid(AminoAcid):

	def __init__():

		self.smilesString = aAmineGroup + aCarbon + sideChain + carboxylGroup

class BetaHomoAminoAcid(HomoAminoAcid, BetaAminoAcid):

	def __init__():

		self.smilesString = aAmineGroup + aCarbon + sideChain + carboxylGroup
"""


def GenerateCanonicalAminoAcids():
    # Alanines
    lAla = AlphaAminoAcid("Ala", aAmineGroup, aL, alaSC, aCarbon, carboxylGroup)
    dAla = AlphaAminoAcid("D-Ala", aAmineGroup, aD, alaSC, aCarbon, carboxylGroup)

    # Aspartic Acids
    lAsp = AlphaAminoAcid("Asp", aAmineGroup, aL, aspSC, aCarbon, carboxylGroup)
    dAsp = AlphaAminoAcid("D-Asp", aAmineGroup, aD, aspSC, aCarbon, carboxylGroup)

    # Cysteines
    lCys = AlphaAminoAcid("Cys", aAmineGroup, aL, cysSC, aCarbon, carboxylGroup)
    dCys = AlphaAminoAcid("D-Cys", aAmineGroup, aD, cysSC, aCarbon, carboxylGroup)

    # Isoleucines
    lIle = AlphaAminoAcid("Ile", aAmineGroup, aL, ileSC, aCarbon, carboxylGroup)
    dIle = AlphaAminoAcid("D-Ile", aAmineGroup, aD, ileSC, aCarbon, carboxylGroup)

    # Leucines
    lLeu = AlphaAminoAcid("Leu", aAmineGroup, aL, leuSC, aCarbon, carboxylGroup)
    dLeu = AlphaAminoAcid("D-Leu", aAmineGroup, aD, leuSC, aCarbon, carboxylGroup)

    # Glutamines
    lGlu = AlphaAminoAcid("Glu", aAmineGroup, aL, gluSC, aCarbon, carboxylGroup)
    dGlu = AlphaAminoAcid("D-Glu", aAmineGroup, aD, gluSC, aCarbon, carboxylGroup)

    # Glycines
    lGly = AlphaAminoAcid("Gly", aAmineGroup, "", glySC, aCarbon, carboxylGroup)
    dGly = AlphaAminoAcid("D-Gly", aAmineGroup, "", glySC, aCarbon, carboxylGroup)

    # Methionines
    lMet = AlphaAminoAcid("Met", aAmineGroup, aL, metSC, aCarbon, carboxylGroup)
    dMet = AlphaAminoAcid("D-Met", aAmineGroup, aD, metSC, aCarbon, carboxylGroup)

    # Phenylalanines
    lPhe = AlphaAminoAcid("Phe", aAmineGroup, aL, pheSC, aCarbon, carboxylGroup)
    dPhe = AlphaAminoAcid("D-Phe", aAmineGroup, aD, pheSC, aCarbon, carboxylGroup)

    # Prolines

    # FIX ME
    lPro = AlphaAminoAcid("Pro", aAmineGroup, aNull, proSC, aCarbon, carboxylGroup)
    dPro = AlphaAminoAcid("D-Pro", aAmineGroup, aNull, dproSC, aCarbon, carboxylGroup)

    # Threonines
    lThr = AlphaAminoAcid("Thr", aAmineGroup, aL, thrSC, aCarbon, carboxylGroup)
    dThr = AlphaAminoAcid("D-Thr", aAmineGroup, aD, thrSC, aCarbon, carboxylGroup)

    # Tryptophans
    lTrp = AlphaAminoAcid("Trp", aAmineGroup, aL, trpSC, aCarbon, carboxylGroup)
    dTrp = AlphaAminoAcid("D-Trp", aAmineGroup, aD, trpSC, aCarbon, carboxylGroup)

    # Tyrosines
    lTyr = AlphaAminoAcid("Tyr", aAmineGroup, aL, tyrSC, aCarbon, carboxylGroup)
    dTyr = AlphaAminoAcid("D-Tyr", aAmineGroup, aD, tyrSC, aCarbon, carboxylGroup)

    # Serines
    lSer = AlphaAminoAcid("Ser", aAmineGroup, aL, serSC, aCarbon, carboxylGroup)
    dSer = AlphaAminoAcid("D-Ser", aAmineGroup, aD, serSC, aCarbon, carboxylGroup)

    # Valines
    lVal = AlphaAminoAcid("Val", aAmineGroup, aL, valSC, aCarbon, carboxylGroup)
    dVal = AlphaAminoAcid("D-Val", aAmineGroup, aD, valSC, aCarbon, carboxylGroup)

    # print(dAla.smilesString)

    return {lAla.name : lAla, dAla.name : dAla,
            lAsp.name : lAsp, dAsp.name : dAsp,
            lCys.name : lCys, dCys.name : dCys,
            lIle.name : lIle, dIle.name : dIle,
            lLeu.name : lLeu, dLeu.name : dLeu,
            lGlu.name : lGlu, dGlu.name : dGlu,
            lGly.name : lGly, dGly.name : dGly,
            lMet.name : lMet, dMet.name : dMet,
            lPhe.name : lPhe, dPhe.name : dPhe,
            lPro.name : lPro, dPro.name : dPro,
            lThr.name : lThr, dThr.name : dThr,
            lTrp.name : lTrp, dTrp.name : dTrp,
            lTyr.name : lTyr, dTyr.name : dTyr,
            lSer.name : lSer, dSer.name : dSer,
            lVal.name : lVal, dVal.name : dVal
            }


def GenerateNoncanonicalAminoAcids():
    # Cyclohexylalanines
    lCha = AlphaAminoAcid("Cha", aAmineGroup, aL, chaSC, aCarbon, carboxylGroup)
    dCha = AlphaAminoAcid("D-Cha", aAmineGroup, aD, chaSC, aCarbon, carboxylGroup)

    return {lCha.name : lCha, dCha.name : dCha
            }

def GenerateSubunits():
    print("\nGenerating subunit library...")

    start = time.time()

    subunits = {}

    subunits.update(GenerateCanonicalAminoAcids())

    subunits.update(GenerateNoncanonicalAminoAcids())

    print("Subunits in library: " + str(len(subunits)))

    end = time.time()

    print("Time elapsed: " + str((end - start) / 60) + " minutes\n")

    return subunits
