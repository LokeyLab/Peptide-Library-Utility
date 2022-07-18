"""
This script contains various miscellaneous utility functions.
"""

from scripts import header as h


def csv_to_dataframe(f):
    """Reads a CSV file to a Pandas dataframe."""
    df = h.pd.read_csv(f)

    return df


def peptides_to_dataframe(cyclicPeptides):
    """Writes the parameters from a list of cyclic peptide objects to a Pandas dataframe."""

    columns = ["Name", "Exact Mass", "TPSA", "ALogP", "SMILES String"]

    df = h.pd.DataFrame(columns=columns)

    for i in range(0,len(cyclicPeptides)):
        df.loc[len(df.index)] = { "Name" : cyclicPeptides[i].name,
                                  "Exact Mass" : cyclicPeptides[i].exactMass,
                                  "TPSA" : cyclicPeptides[i].TPSA,
                                  "ALogP" : cyclicPeptides[i].ALogP,
                                  "SMILES String" : cyclicPeptides[i].smiles_string}

    return df


def print_dataframe(df):
    """Prints a Pandas dataframe."""

    print("\n" + df.to_string())


def dataframe_to_csv(df):
    """Writes a Pandas dataframe to a CSV file in the output folder"""

    df.to_csv("output/output.csv")


def plot_exact_mass_tpsa_alogp(df):
    """Displays a 3D scatter plot showing Exact Mass, TPSA, and ALogP"""

    fig = h.plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')

    x = df["Exact Mass"]
    y = df["TPSA"]
    z = df["ALogP"]
    p = (x + y + z) # df["Predicted PappE-6"]

    ax.scatter(x, y, z, c = p, cmap="viridis", edgecolors='black', s=20)

    ax.plot(x, z, "r+", zdir='y', zs=y.max())
    ax.plot(y, z, "g+", zdir='x', zs=x.min())
    ax.plot(x, y, "b+", zdir='z', zs=z.min())

    ax.set_xlabel('Exact Mass')
    ax.set_ylabel('TPSA')
    ax.set_zlabel('ALogP')

    h.plt.show()
