import py3Dmol as py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
import rdkit.DistanceGeometry as DG
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.ForceField import rdForceField

import matplotlib.pyplot as plt

from tqdm import tqdm

import statistics


smiles = ["N4(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O",
"N4(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O",
"N4(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N1[C@H](CCC1)C(=O)N(Cc2ccccc2)CC(=O)N[C@@H](Cc3ccccc3)C4=O",
"N5(C)[C@H](C)C(=O)N[C@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N5(C)[C@H](C)C(=O)N[C@@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N4(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O",
"N4(C)[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O",
"N5(C)[C@H](C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N5(C)[C@H](C)C(=O)N(Cc1ccccc1)CC(=O)N[C@@H](CC2=CC=CC=C2)CC(=O)N3[C@H](CCC3)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc4ccccc4)C5=O",
"N5(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N4(C)[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O",
"N5(C)[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N5(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)CC(=O)N2[C@H](CCC2)C(=O)N(Cc3ccccc3)CC(=O)N[C@@H](Cc4ccccc4)C5=O",
"N5(C)[C@H](C)C(=O)N(Cc1ccccc1)CC(=O)N[C@@H](CC2=CC=CC=C2)CC(=O)N3[C@H](CCC3)C(=O)N(C)[C@H](CC(C)C)C(=O)N[C@@H](Cc4ccccc4)C5=O",
"N6(C)[C@H](C)C(=O)N(Cc1ccccc1)CC(=O)N(Cc2ccccc2)CC(=O)N3[C@H](CCC3)C(=O)N(Cc4ccccc4)CC(=O)N[C@@H](Cc5ccccc5)C6=O"
]

mols = []

energiesInWater = []
energiesInChloroform = []

mol = Chem.MolFromSmiles(smiles[0])

params = AllChem.ETKDGv3()

# Water
for i in tqdm(range(0, 72)):
    mol = Chem.AddHs(mol)

    Chem.SanitizeMol(mol)

    AllChem.EmbedMultipleConfs(mol, 1, params)

    AllChem.MMFFOptimizeMolecule(mol)

    mp = AllChem.MMFFGetMoleculeProperties(mol)

    mp.SetMMFFDielectricConstant(dielConst=80.1)

    ffm = AllChem.MMFFGetMoleculeForceField(mol, mp)

    energy = ffm.CalcEnergy()

    energiesInWater.append(energy)

    mol = Chem.RemoveHs(mol)

    mols.append(mol)

#Chloroform
for i in tqdm(range(0, 72)):
    mol = Chem.AddHs(mol)

    Chem.SanitizeMol(mol)

    AllChem.EmbedMultipleConfs(mol, 1, params)

    AllChem.MMFFOptimizeMolecule(mol)

    mp = AllChem.MMFFGetMoleculeProperties(mol)

    mp.SetMMFFDielectricConstant(dielConst=4.81)

    ffm = AllChem.MMFFGetMoleculeForceField(mol, mp)

    energy = ffm.CalcEnergy()

    energiesInChloroform.append(energy)

    mol = Chem.RemoveHs(mol)

    mols.append(mol)

img = MolsToGridImage(mols[:], molsPerRow=16, subImgSize=(100, 100))
img.save(fp="Conformers.png")
img.show()

def myFunc(e):
    return e['energy']

# mols.sort(key=myFunc)

# energiesInWater.sort()

# print(mols)

plt.violinplot([energiesInWater, energiesInChloroform], showextrema=True, showmeans=True)
plt.title("Energy of Conformers")
plt.ylabel("Energy (kcal/mol)")
plt.xlabel("Water vs. Chloroform")
# plt.axhline(y=statistics.median(energys), color='r', linestyle='--')

plt.savefig("Energy of Conformers.png")
plt.show()


