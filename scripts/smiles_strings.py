"""
This script is not used in the program, but rather, just
holds the SMILES strings of the canonical amino acids which
had to be rearranged from their original form that was
taken from Pub Chem. The Pub Chem version is commented next
to each amino acid.

Standardized format:
α-amine group | carbon atom C2 | ( side chain ) | α-carboxyl group

Taken from:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149970/

Template:
N[C@@H]()C(=O)O
"""


ALA = "N[C@@H](C)C(=O)O"  # C[C@@H](C(=O)O)N
ARG = "N[C@@H](CCCN=C(N)N)C(=O)O"  # C(C[C@@H](C(=O)O)N)CN=C(N)N
ASN = "N[C@@H](CC(=O)N)C(=O)N"  # C([C@@H](C(=O)O)N)C(=O)N
ASP = "N[C@@H](CC(=O)O)C(=O)O"  # C([C@@H](C(=O)O)N)C(=O)O
CYS = "N[C@@H](CS)C(=O)O"  # C([C@@H](C(=O)O)N)S
HIS = "N[C@@H](CC1=C(NC=N1))C(=O)O"  # C1=C(NC=N1)C[C@@H](C(=O)O)N
ILE = "N[C@@H]([C@@H](C)CC)C(=O)O"  # CC[C@H](C)[C@@H](C(=O)O)N
LEU = "N[C@@H](CC(C)C)C(=O)O"  # CC(C)C[C@@H](C(=O)O)N
LYS = "N[C@@H](CCCCN)C(=O)O"  # C(CCN)C[C@@H](C(=O)O)N
GLN = "N[C@@H](CCC(=O)N)C(=O)O"  # C(CC(=O)N)[C@@H](C(=O)O)N
GLU = "N[C@@H](C(CC(=O)O))C(=O)O"  # C(CC(=O)O)[C@@H](C(=O)O)N
GLY = "NCC(=O)O"  # C(C(=O)O)N
MET = "N[C@@H](CSCC)C(=O)O"  # CSCC[C@@H](C(=O)O)N
PHE = "N[C@@H](Cc1ccccc1)C(=O)O"  # C1=CC=C(C=C1)C[C@@H](C(=O)O)N
PRO = "N1CCC[C@H]1C(=O)O"  # C1C[C@H](NC1)C(=O)O
THR = "N[C@@H]([C@H](O)C)C(=O)O"  # C[C@H]([C@@H](C(=O)O)N)O
TRP = "N[C@@H](Cc1=c[nH]c2c1cccc2)C(=O)O"  # C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N
TYR = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"  # C1=CC(=CC=C1C[C@@H](C(=O)O)N)O
SER = "N[C@@H](C(O))C(=O)O"  # C([C@@H](C(=O)O)N)O
VAL = "N[C@@H](C(C)C)C(=O)O"  # CC(C)[C@@H](C(=O)O)N
