import matplotlib.pyplot as plt
from matplotlib import font_manager
import networkx as nx


class Atom:
    """Atom class"""
    def __init__(self, node_index=-1, atom_type="Z"):
        """Initializes atom object"""
        self.node_index = node_index
        self.atom_type = atom_type



class Electron:
    """Electron class"""
    def __init__(self, node_index=-1, spin=0):
        """Initializes electron object"""
        self.node_index = node_index
        self.spin = spin


G = nx.Graph()

atoms = []
labels = []
G.add_node(0)
G.add_node(1)
G.add_node(2)
G.add_node(3)

G.add_edge(0,1)
G.add_edge(1,2)
G.add_edge(1,3)
G.add_edge(3,1)


print("Number of Nodes = %i" % G.number_of_nodes())
print("Number of Edges = %i" % G.number_of_edges())

nx.draw_networkx(G, node_color=['black', 'black', 'black', 'red'], labels={0:'C', 1: "C", 2: "C",
                                                                           3:"O"},
                 font_color='white', font_weight="bold")

plt.tight_layout()

plt.show()
