# internal
import bonding as b

# External
import itertools
import time


def CartesianProduct(aminos, aminosPerPeptide):
    listToIterate = []

    for i in range(0, len(aminos)):
        listToIterate.append(i)

    cartesianProduct = []

    for j in itertools.product(listToIterate, repeat=aminosPerPeptide):
        cartesianProduct.append(j)

    return cartesianProduct


def AminoAcidCartesianProduct(aminos, aminosPerPeptide):
    peptides = []

    print("Finding Cartesian Product...")

    start = time.time()

    cartesianProduct = CartesianProduct(aminos, aminosPerPeptide)

    print("Combinations found: " + str(len(cartesianProduct)))

    end = time.time()

    print("Time elapsed: " + str((end - start) / 60) + " minutes\n")

    print("Bonding Residues...")

    start = time.time()

    for i in range(0, len(cartesianProduct)):

        peptide = []

        for j in range(0, aminosPerPeptide):

            if peptide not in peptides:
                peptide.append(aminos[cartesianProduct[i][j]].molecule)

        peptides.append(b.BondHeadToTail(peptide))

    end = time.time()

    print("Time elapsed: " + str((end - start) / 60) + " minutes\n")

    return peptides
