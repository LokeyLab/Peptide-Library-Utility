import numpy as np
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt


class molMatricies:
    def __init__(self, adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix, distanceMatrix):

        self.adjacencyMatrix = adjacencyMatrix
        self.atomicNumberMatrix = atomicNumberMatrix
        # self.formalChargeMatrix = formalChargeMatrix
        self.aromaticityMatrix = aromaticityMatrix
        self.chiralityMatrix = chiralityMatrix
        self.distanceMatrix = distanceMatrix

def molToNx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomicNum = atom.GetAtomicNum(),
                   isAromatic = atom.GetIsAromatic(),
                   chirality = atom.GetChiralTag(),
                   isomer = atom.GetIsotope(),
                   atomSymbol = atom.GetSymbol(),
                   atomCharge = atom.GetFormalCharge())

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())

    return G

def molToAtomicNumberMatrix(molNx):
    atomsList = []

    for i in range(0, len(molNx)):
        atomList = []
        atomList.append(np.zeros(i))
        atomList.append([float(molNx.nodes.data('atomicNum')[i])])
        atomList.append(np.zeros(len(molNx) - len(atomList[0]) - 1))
        flatList = [item for sublist in atomList for item in sublist]
        atomList = flatList
        atomsList.append(atomList)

    atomicNumberMatrix = np.matrix(atomsList)

    return atomicNumberMatrix

def molToFormalChargeMatrix(molNx):
    atomsList = []

    for i in range(0, len(molNx)):
        atomList = []
        atomList.append(np.zeros(i))
        atomList.append([float(molNx.nodes.data('atomCharge')[i])])
        atomList.append(np.zeros(len(molNx) - len(atomList[0]) - 1))
        flatList = [item for sublist in atomList for item in sublist]
        atomList = flatList
        atomsList.append(atomList)

    formalChargeMatrix = np.matrix(atomsList)

    return formalChargeMatrix

def molToAromaticityMatrix(molNx):
    atomsList = []

    for i in range(0, len(molNx)):
        atomList = []
        atomList.append(np.zeros(i))
        atomList.append([float(molNx.nodes.data('isAromatic')[i])])
        atomList.append(np.zeros(len(molNx) - len(atomList[0]) - 1))
        flatList = [item for sublist in atomList for item in sublist]
        atomList = flatList
        atomsList.append(atomList)

    aromaticityMatrix = np.matrix(atomsList)

    return aromaticityMatrix

def molToChiralityMatrix(molNx):
    atomsList = []

    for i in range(0, len(molNx)):
        atomList = []
        atomList.append(np.zeros(i))
        atomList.append([float(molNx.nodes.data('chirality')[i])])
        atomList.append(np.zeros(len(molNx) - len(atomList[0]) - 1))
        flatList = [item for sublist in atomList for item in sublist]
        atomList = flatList
        atomsList.append(atomList)

    chiralityMatrix = np.matrix(atomsList)

    return chiralityMatrix

def smilesToMatricies(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    molH = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(molH)
    AllChem.MMFFOptimizeMolecule(molH)

    molNx = molToNx(molH)

    atomicNumberMatrix = molToAtomicNumberMatrix(molNx)
    # formalChargeMatrix = molToFormalChargeMatrix(molNx)
    aromaticityMatrix = molToAromaticityMatrix(molNx)
    chiralityMatrix = molToChiralityMatrix(molNx)
    distanceMatrix = np.matrix(AllChem.Get3DDistanceMatrix(molH))
    adjacencyMatrix = nx.to_numpy_matrix(molNx)

    matricies = molMatricies(adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix, distanceMatrix)

    # return matricies
    return molNx, adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix, distanceMatrix

def plotCyclicPeptideMatricies(molNx, adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix,
                               distanceMatrix):
    atom = nx.get_node_attributes(molNx, 'atomSymbol')

    color_map = {'H': '#DDDDDD', 'C': 'gray', 'N': 'blue', 'O': 'red', 'Cl': 'green', 'S': 'yellow', 'P': 'orange'}

    mol_colors = []

    for idx in molNx.nodes():
        if (molNx.nodes[idx]['atomSymbol'] in color_map):
            mol_colors.append(color_map[molNx.nodes[idx]['atomSymbol']])
        else:
            mol_colors.append('#A62929')

    fig, axs = plt.subplots(2,3)

    axs[0,0].set_title("Graph")
    nx.draw(molNx,
            labels=atom,
            with_labels=True,
            font_size = 4,
            node_color=mol_colors,
            node_size=25, ax = axs[0,0])

    axs[0,1].set_title("Adjacency")
    axs[0,1].matshow(adjacencyMatrix)

    axs[0,2].set_title("Atomic Number")
    axs[0,2].matshow(atomicNumberMatrix)

    # axs[1,0].set_title("Formal Charge")
    # axs[1,0].matshow(formalChargeMatrix)

    axs[1,0].set_title("Aromaticity")
    axs[1,0].matshow(aromaticityMatrix)

    axs[1,1].set_title("Chirality")
    axs[1,1].matshow(chiralityMatrix)

    axs[1,2].set_title("3D Distance")
    axs[1,2].matshow(distanceMatrix)

    plt.tight_layout()
    plt.savefig("Data Structure.png")

    plt.show()

"""
smiles = "N4(C)[C@H](C)C(=O)N(C)[C@H](CC(C)C)C(=O)N(Cc1ccccc1)CC(=O)N2[C@H](CCC2)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](Cc3ccccc3)C4=O"

molNx, adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix, distanceMatrix = smilesToMatricies(smiles)

plotCyclicPeptideMatricies(molNx, adjacencyMatrix, atomicNumberMatrix, aromaticityMatrix, chiralityMatrix,
                               distanceMatrix)
"""