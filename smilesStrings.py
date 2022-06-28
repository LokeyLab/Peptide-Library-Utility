# AMINO ACIDS

# Standardized format: α-amine group | carbon atom C2 | ( side chain ) | α-carboxyl group
# Template: N[C@@H]()C(=O)O

ala = "N[C@@H](C)C(=O)O"  # C[C@@H](C(=O)O)N
arg = "N[C@@H](CCCN=C(N)N)C(=O)O" # C(C[C@@H](C(=O)O)N)CN=C(N)N
asn = "N[C@@H](CC(=O)N)C(=O)N"  # C([C@@H](C(=O)O)N)C(=O)N
asp = "N[C@@H](CC(=O)O)C(=O)O" # C([C@@H](C(=O)O)N)C(=O)O
cys = "N[C@@H](CS)C(=O)O"  # C([C@@H](C(=O)O)N)S
his = "N[C@@H](CC1=C(NC=N1))C(=O)O" # C1=C(NC=N1)C[C@@H](C(=O)O)N
ile = "N[C@@H]([C@@H](C)CC)C(=O)O"  # CC[C@H](C)[C@@H](C(=O)O)N
leu = "N[C@@H](CC(C)C)C(=O)O"  # CC(C)C[C@@H](C(=O)O)N # N[C@@H](CC(C)C)C(=O)O # CC(C)[C@@H](N)CCC(O)=O
lys = "N[C@@H](CCCCN)C(=O)O" # C(CCN)C[C@@H](C(=O)O)N
gln = "N[C@@H](CCC(=O)N)C(=O)O"  # C(CC(=O)N)[C@@H](C(=O)O)N
glu = "N[C@@H](C(CC(=O)O))C(=O)O" # C(CC(=O)O)[C@@H](C(=O)O)N
gly = "NCC(=O)O"  # C(C(=O)O)N
met = "N[C@@H](CSCC)C(=O)O"  # CSCC[C@@H](C(=O)O)N
phe = "N[C@@H](Cc1ccccc1)C(=O)O"  # C1=CC=C(C=C1)C[C@@H](C(=O)O)N
pro = "N1CCC[C@H]1C(=O)O"  # C1C[C@H](NC1)C(=O)O
thr = "N[C@@H]([C@H](O)C)C(=O)O"  # C[C@H]([C@@H](C(=O)O)N)O
trp = "N[C@@H](Cc1=c[nH]c2c1cccc2)C(=O)O"  # C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N
tyr = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"  # C1=CC(=CC=C1C[C@@H](C(=O)O)N)O
ser = "N[C@@H](C(O))C(=O)O"  # C([C@@H](C(=O)O)N)O
val = "N[C@@H](C(C)C)C(=O)O"  # CC(C)[C@@H](C(=O)O)N

# NC(CC)C(=O)O
