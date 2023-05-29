from rdkit import Chem
from rdkit.Chem import Draw
import networkx as nx
import matplotlib.pyplot as plt


CARBOXYLIC_ACID = "(=O)O"
ALPHA_CARBON = "C"
CHIRAL_CARBON = "[CH]"
CW_CHIRAL_CARBON = "[C@@H]"
CCW_CHIRAL_CARBON = "[C@H]"

AMINE = "N"

side_chains = {
"ala": "C",
"gly": "([H])",
"bip": "(Cc1ccc(c2ccccc2)cc1)",
"leu": "(CC(C)(C))",
"phe": "(Cc1ccccc1)",
"val": "(C(C)(C)"
}

amino_acid = AMINE + CW_CHIRAL_CARBON + side_chains["ala"] + ALPHA_CARBON + CARBOXYLIC_ACID

"""def get_graph(mol):
  Chem.Kekulize(mol)
  atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
  am = Chem.GetAdjacencyMatrix(mol,useBO=True)
  for i,atom in enumerate(atoms):
    am[i,i] = atom
  G = nx.from_numpy_matrix(am)
  return G

print(amino_acid)

mol = Chem.MolFromSmiles(amino_acid)

G = get_graph(mol)

nx.draw_networkx(G)  # pyplot draw()
plt.show()"""

mol = Chem.MolFromSmiles(amino_acid)
mol_img = Draw.MolToImage(mol)
mol_img.show()
