#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains various miscellaneous utility functions.
"""

from scripts import header as h


def csv_to_dataframe(f):
    """Reads a CSV file to a Pandas dataframe."""

    df = h.pd.read_csv(f)

    return df


def peptides_to_dataframe(peptides):
    """Writes the parameters from a list of cyclic peptide objects to a Pandas dataframe."""

    df = h.pd.DataFrame(columns=["Name"] + h.CHEMINFORMATICS_OPTIONS + ["SMILES String"])

    for i in range(0, len(peptides)):

        # "Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "McGowan Volume",
        # "# Atoms", "# Heteroatoms", "# Amide Bonds",
        # "# Rings", "# Aromatic Rings", "# Rotatable Bonds",
        # "# Stereocenters", "# HBA", "# HBD", "Largest Ring Size"

        row_dict = {}
        
        row_dict.update({"Name": peptides[i].name})
        row_dict.update({"SMILES String": peptides[i].smiles_string})
        if "Exact Mass" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"Exact Mass": peptides[i].exact_mass})
        if "ALogP" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"ALogP": peptides[i].a_log_p})
        if "SLogP" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"SLogP": peptides[i].s_log_p})
        if "SASA" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"SASA": peptides[i].sasa})
        if "TPSA" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"TPSA": peptides[i].tpsa})
        if "McGowan Volume" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"McGowan Volume": peptides[i].mcgowan_volume})
        if "# Atoms" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Atoms": peptides[i].num_atoms})
        if "# Heteroatoms" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Heteroatoms": peptides[i].num_heteroatoms})
        if "# Amide Bonds" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Amide Bonds": peptides[i].num_amide_bonds})
        if "# Rings" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Rings": peptides[i].num_rings})
        if "# Aromatic Rings" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Aromatic Rings": peptides[i].num_aromatic_rings})
        if "# Rotatable Bonds" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Rotatable Bonds": peptides[i].num_rotatable_bonds})
        if "# Stereocenters" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# Stereocenters": peptides[i].num_stereocenters})
        if "# HBA" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# HBA": peptides[i].num_hba})
        if "# HBD" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"# HBD": peptides[i].num_hbd})
        if "Largest Ring Size" in h.CHEMINFORMATICS_OPTIONS:
            row_dict.update({"Largest Ring Size": peptides[i].largest_ring_size})

        df.loc[len(df.index)] = row_dict

    return df

def subunits_to_dataframe():
    """Writes the subunit library to a Pandas dataframe."""

    df = h.pd.DataFrame(
        columns=["Name", "Multiple Letter", "SMILES String"])

    for i in range(len(h.SUBUNIT_LIBRARY)):
        key = list(h.SUBUNIT_LIBRARY)[i]
        row_dict = {}

        row_dict.update({"Name": h.SUBUNIT_LIBRARY[key].name})
        row_dict.update({"Multiple Letter": h.SUBUNIT_LIBRARY[key].multiple_letter})
        row_dict.update({"SMILES String": h.SUBUNIT_LIBRARY[key].smiles_string})

        df.loc[len(df.index)] = row_dict

    return df

def molecules_to_dataframe(molecules_list):

    columns = ["Molecules"]
    df = h.pd.DataFrame(columns=columns)

    for i in range(0,len(molecules_list)):
        df.loc[len(df.index)] = {"Molecules": molecules_list[i]}

    return df


def print_dataframe(df):
    """Prints a Pandas dataframe."""

    print("\n" + df.to_string())


# ".CSV", ".JSON", ".SQL", ".TEX", ".TXT", ".XLSX",
# ".XML", ".SDF (Only Mol Data)", ".SMI (Only SMILES)"


def dataframe_to_csv(df):
    """Writes a Pandas dataframe to a .CSV file in the output folder"""

    file_name = "output_data.csv"
    path = h.os.path.join(h.OUTPUT_DIR, file_name)

    df.to_csv(path)

def subunits_to_csv(df):
    """Writes a Pandas dataframe to a .CSV file in the output folder"""

    file_name = "subunit_library.csv"
    path = h.os.path.join(h.OUTPUT_DIR, file_name)

    df.to_csv(path)


def dataframe_to_txt(df):
    """Writes a Pandas dataframe to a .TXT file in the output folder"""

    file_name = "output_data.txt"
    path = h.os.path.join(h.OUTPUT_DIR, file_name)

    df.to_csv(path)


def dataframe_to_sdf(df):
    """Writes a Pandas dataframe to an .SDF file in the output folder"""

    file_name = "output_molecules.sdf"
    path = h.os.path.join(h.OUTPUT_DIR, file_name)

    h.PandasTools.WriteSDF(df, path, molColName="Molecules", properties=list(df.columns))


def dataframe_to_smi(df):
    """Writes a Pandas dataframe to a .SMI file in the output folder"""

    file_name = "output_smiles.smi"
    path = h.os.path.join(h.OUTPUT_DIR, file_name)

    df["SMILES String"].to_csv(path)


def plot_exact_mass_tpsa_a_log_p(x, y, z, c, cmap):
    """Displays a 3D scatter plot showing Exact Mass, TPSA, and a_log_p"""

    fig = h.plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')

    ax.scatter(x, y, z, c = c, cmap=cmap, edgecolors='black', s=20)

    ax.plot(x, z, "r+", zdir='y', zs=y.max())
    ax.plot(y, z, "g+", zdir='x', zs=x.min())
    ax.plot(x, y, "b+", zdir='z', zs=z.min())

    ax.set_xlabel('Exact Mass')
    ax.set_ylabel('TPSA')
    ax.set_zlabel('ALogP')

    h.plt.show()
