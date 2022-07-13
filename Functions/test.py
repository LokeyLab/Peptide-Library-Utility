from rdkit.Chem import AllChem
from rdkit import Chem
import networkx as nx
import matplotlib.pyplot as plt

# define the smiles string and covert it into a molecule sturcture ------------
smiles = 'N4[C@@H](CC1=CC=CC=C1)CC(=O)N[C@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N2[C@H](CCC2)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O'
mol = Chem.MolFromSmiles(smiles)
molH = AllChem.AddHs(mol)
AllChem.EmbedMolecule(molH)
AllChem.MMFFOptimizeMolecule(molH)
distMat3D = AllChem.Get3DDistanceMatrix(molH)
adjMat = AllChem.GetAdjacencyMatrix(molH)

print(Chem.MolToXYZBlock(molH))
print(adjMat)


# define the function for coverting rdkit object to networkx object -----------
def mol_to_nx(mol):
    G = nx.Graph()

    nodes = []
    edges = []
    adjacencies = []

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   atom_symbol=atom.GetSymbol())
        nodes.append(atom.GetAtomicNum())

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
        edges.append(bond.GetBondType())
        adjacencies.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

    print("Atoms = " + str(nodes))
    print("Edges = " + str(edges))
    print("Adjacencies = " + str(adjacencies))

    return G


# conver rdkit object to networkx object --------------------------------------
caffeine_nx = mol_to_nx(molH)

caffeine_atom = nx.get_node_attributes(caffeine_nx, 'atom_symbol')

color_map = {'H': '#DDDDDD',
             'C': 'gray',
             'O': 'red',
             'N': 'blue'}

caffeine_colors = []
for idx in caffeine_nx.nodes():
    if (caffeine_nx.nodes[idx]['atom_symbol'] in color_map):
        caffeine_colors.append(color_map[caffeine_nx.nodes[idx]['atom_symbol']])
    else:
        caffeine_colors.append('gray')

nx.draw(caffeine_nx,
        labels=caffeine_atom,
        with_labels=True,
        node_color=caffeine_colors,
        node_size=200)

# plt.show()
# print out the adjacency matrix ----------------------------------------------
adjMatrix = nx.to_numpy_matrix(caffeine_nx)
plt.matshow(adjMatrix)
plt.show()
print(adjMatrix)

