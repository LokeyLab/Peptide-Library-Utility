#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains all the functions used for the GUI.
"""

from scripts import header as h


def single_peptide_create_button(single_peptide_window,
                                 console_list, subunit_library,
                                 subunit_list, single_peptide_progress_bar):
    """
    try:
        single_peptide_progress_bar.start(10)
        console_list.insert("end", "Peptide Creation Succeeded.")
        console_list.see("end")
        single_peptide_progress_bar.stop()

    except KeyError:
        console_list.insert("end", "Peptide Creation Failed")
        console_list.see("end")


    finally:
        single_peptide_window.destroy()
        console_list.insert("end", "Single Peptide Creator Exited")
        console_list.see("end")
    """


def subunit_scan_key(event):
    val = event.widget.get()

    # print(val)

    if val == '':

        data = h.SUBUNIT_STRING_DATA
    else:

        data = []
        for item in h.SUBUNIT_STRING_DATA:
            if val.lower() in item.lower():
                data.append(item)

    subunit_update(data)


def subunit_update(data):
    h.SUBUNIT_LIST.clear()

    for item in data:
        h.SUBUNIT_LIST.append((item))
        h.SUBUNIT_COUNT = f"n = {len(data)}"


def subunit_library_update(console_list, subunit_list):
    h.SUBUNIT_LIBRARY = h.subunit_builder.generate_subunit_library()
    h.SUBUNIT_COUNT = len(h.SUBUNIT_LIBRARY)

    console_list.insert("end", "Subunit Library Populated.")
    console_list.insert("end", "")
    console_list.see("end")

    subunit_list.delete(0, "end")

    h.SUBUNIT_STRING_DATA.clear()

    for i in range(0, len(h.SUBUNIT_LIBRARY)):
        key = list(h.SUBUNIT_LIBRARY)[i]

        h.SUBUNIT_STRING_DATA.append(h.SUBUNIT_LIBRARY[key].name
                                     + "        " +
                                     h.SUBUNIT_LIBRARY[key].multiple_letter)

        subunit_list.insert("end", h.SUBUNIT_STRING_DATA[i])

    return len(h.SUBUNIT_LIBRARY)



def peptide_library_creation(console_list, output_list, subunit_list,
                             peptide_library_progress_bar, cyclization_type,
                             amine_cap_subunit, c_terminus_cap_subunit,
                             peptide_library_window, subunit_count):

    if not h.SUBUNIT_LIBRARY_POPULATED:
        subunit_library_update(console_list, subunit_list)
        subunit_count.set("n = " + str(h.SUBUNIT_COUNT))


    subunit_library = h.SUBUNIT_LIBRARY

    input_file = h.os.path.abspath(h.INPUT_FILE)

    if h.os.path.exists(input_file):
        df = h.utilities.csv_to_dataframe(input_file)
        console_list.insert("end", f"Input File: {input_file}")
        console_list.insert("end", "")
        console_list.see("end")

    else:
        h.messagebox.showinfo("Error Message", "Input file does not exist.")

    # h.utilities.print_dataframe(df)
    df_list = df.to_string().split("\n")

    for i in range(0, df.count().max() + 1):
        # row_list = df_list[i]
        print(df_list[i])
        console_list.insert("end", df_list[i])
        console_list.see("end")

    console_list.insert("end", "")

    df = df.iloc[:, ::-1]

    pot_sizes = []

    pots = []

    for i in range(0, df.shape[1]):

        temp_pot = []

        count = 0

        for j in range(0, df.shape[0]):

            if h.pd.notnull(df.iloc[j, i]):
                temp_pot.append(df.iloc[j, i])

                count += 1

        pot_sizes.append(count)

        pots.append(temp_pot)

    peptide_count = pot_sizes[0]

    for k in range(1, len(pot_sizes)):
        peptide_count *= pot_sizes[k]

    console_list.insert("end", "Peptide Count = " + str(peptide_count))
    console_list.see("end")

    cartesian_product = h.combinatronics.cartesian_product(pots)

    peptides = []

    molecules_list = []

    for i in cartesian_product:

        peptide_library_window.update_idletasks()

        subunits = []

        for j in i:

            try:
                subunit = subunit_library[j]

                subunits.append(subunit)

            except KeyError:

                h.messagebox.showinfo("Error Message", f"There was an issue with subunit: {j}.")

                break

        peptide = h.bonding.bond_head_to_tail(subunits)

        if cyclization_type == "Head to Tail":

            peptide = h.bonding.cyclize_head_to_tail(peptide)

        elif cyclization_type == "Depsipeptide":

            peptide = h.bonding.cyclize_depsipeptide(peptide, subunit_library)

            if not amine_cap_subunit == '':
                subunits.append(h.SUBUNIT_LIBRARY[amine_cap_subunit])

                tempPeptide = str(peptide[0]) + "(" + \
                              h.SUBUNIT_LIBRARY[amine_cap_subunit].smiles_string + ")" + \
                              str(peptide[1:])

                peptide = tempPeptide

        elif cyclization_type == "Thioether":

            peptide = h.bonding.cyclize_thioether(peptide, subunit_library)

            if not c_terminus_cap_subunit == '':
                subunits.append(h.SUBUNIT_LIBRARY[c_terminus_cap_subunit])

                tempPeptide = str(peptide[0:-1]) + \
                              h.SUBUNIT_LIBRARY[c_terminus_cap_subunit].smiles_string

                peptide = tempPeptide

        elif cyclization_type == "Triazole":

            peptide = h.bonding.cyclize_triazole(peptide, subunit_library)

            if not c_terminus_cap_subunit == '':
                subunits.append(h.SUBUNIT_LIBRARY[c_terminus_cap_subunit])

                tempPeptide = str(peptide[0:-1]) + \
                    h.SUBUNIT_LIBRARY[c_terminus_cap_subunit].smiles_string

                peptide = tempPeptide

        elif cyclization_type == "Linear":

            if not amine_cap_subunit == '':
                subunits.append(h.SUBUNIT_LIBRARY[amine_cap_subunit])

                tempPeptide = str(peptide[0]) + "(" + \
                              h.SUBUNIT_LIBRARY[amine_cap_subunit].smiles_string + ")" + \
                              str(peptide[1:])

                peptide = tempPeptide

            if not c_terminus_cap_subunit == '':
                subunits.append(h.SUBUNIT_LIBRARY[c_terminus_cap_subunit])

                tempPeptide = str(peptide[0:-1]) + \
                              h.SUBUNIT_LIBRARY[c_terminus_cap_subunit].smiles_string

                peptide = tempPeptide

            else:
                peptide = peptide

        print(peptide)

        name = ""

        subunits.reverse()

        for j in range(0, len(subunits)):

            name += subunits[j].multiple_letter + " "

        molecule = h.cheminformatics.get_molecule_from_smiles(peptide)

        # "Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "McGowan Volume",
        # "# Atoms", "# Heteroatoms", "# Amide Bonds",
        # "# Rings", "# Aromatic Rings", "# Rotatable Bonds",
        # "# HBA", "# HBD", "Largest Ring Size"

        if not h.CHEMINFORMATICS_OPTIONS == '':
            if "Exact Mass" in h.CHEMINFORMATICS_OPTIONS:
                exact_mass = h.cheminformatics.get_exact_mass(molecule)
            else:
                exact_mass = None

            if "ALogP" in h.CHEMINFORMATICS_OPTIONS:
                a_log_p = h.cheminformatics.get_a_log_p_99(molecule)
            else:
                a_log_p = None

            if "SLogP" in h.CHEMINFORMATICS_OPTIONS:
                s_log_p = h.cheminformatics.get_s_log_p(molecule)
            else:
                s_log_p = None

            if "SASA" in h.CHEMINFORMATICS_OPTIONS:
                sasa = h.cheminformatics.get_sasa(molecule)
            else:
                sasa = None

            if "TPSA" in h.CHEMINFORMATICS_OPTIONS:
                tpsa = h.cheminformatics.get_tpsa(molecule)
            else:
                tpsa = None

            if "McGowan Volume" in h.CHEMINFORMATICS_OPTIONS:
                mcgowan_volume = h.cheminformatics.get_mcgowan_volume(molecule)
            else:
                mcgowan_volume = None

            if "# Atoms" in h.CHEMINFORMATICS_OPTIONS:
                num_atoms = h.cheminformatics.get_num_atoms(molecule)
            else:
                num_atoms = None

            if "# Heteroatoms" in h.CHEMINFORMATICS_OPTIONS:
                num_heteroatoms = h.cheminformatics.get_num_heteroatoms(molecule)
            else:
                num_heteroatoms = None

            if "# Amide Bonds" in h.CHEMINFORMATICS_OPTIONS:
                num_amide_bonds = h.cheminformatics.get_num_amide_bonds(molecule)
            else:
                num_amide_bonds = None

            if "# Rings" in h.CHEMINFORMATICS_OPTIONS:
                num_rings = h.cheminformatics.get_num_rings(molecule)
            else:
                num_rings = None

            if "# Aromatic Rings" in h.CHEMINFORMATICS_OPTIONS:
                num_aromatic_rings = h.cheminformatics.get_num_aromatic_rings(molecule)
            else:
                num_aromatic_rings = None

            if "# Rotatable Bonds" in h.CHEMINFORMATICS_OPTIONS:
                num_rotatable_bonds = h.cheminformatics.get_num_rotatable_bonds(molecule)
            else:
                num_rotatable_bonds = None

            if "# HBA" in h.CHEMINFORMATICS_OPTIONS:
                num_hba = h.cheminformatics.get_num_hba(molecule)
            else:
                num_hba = None

            if "# HBD" in h.CHEMINFORMATICS_OPTIONS:
                num_hbd = h.cheminformatics.get_num_hbd(molecule)
            else:
                num_hbd = None

            if "Largest Ring Size" in h.CHEMINFORMATICS_OPTIONS:
                largest_ring_size = h.cheminformatics.get_largest_ring_size(molecule)
            else:
                largest_ring_size = None

        peptide_object = h.classes.Peptide(name, peptide, exact_mass, a_log_p, s_log_p, sasa,
                                           tpsa, mcgowan_volume,
                                           num_atoms, num_heteroatoms, num_amide_bonds, num_rings,
                                           num_aromatic_rings, num_rotatable_bonds, num_hba,
                                           num_hbd, largest_ring_size)

        peptides.append(peptide_object)

        molecules_list.append(molecule)

    df = h.utilities.peptides_to_dataframe(peptides)

    df.dropna(how='all', axis=1, inplace=True)

    df2 = h.pd.concat([h.utilities.molecules_to_dataframe(molecules_list), df], axis=1)

    df_list = df.to_string().split("\n")

    date_time = h.datetime.datetime.now().strftime("%m%d%Y%H%M")

    h.OUTPUT_DIR = h.os.path.join(h.OUTPUT_DIR, date_time)

    if not len(peptides) == peptide_count:
        h.messagebox.showinfo("Error Message", f"Number of peptides ({len(peptides)}) output does "
                                               f"not equal the Cartesian product "
                                               f"({peptide_count}). Only valid molecules were "
                                               "output. There may have been an issue with one or "
                                               "more subunits.")

    temp_num = 2

    if h.os.path.exists(h.OUTPUT_DIR):
        temp_path = h.OUTPUT_DIR + "_" + str(temp_num)
        while h.os.path.exists(temp_path):
            temp_num += 1
            temp_path = h.OUTPUT_DIR + "_" + str(temp_num)
        h.OUTPUT_DIR = temp_path
        h.os.mkdir(h.OUTPUT_DIR)

    else:
        h.os.mkdir(h.OUTPUT_DIR)

    for i in range(0, len(df)):
        output_list.insert("end", df_list[i])
        output_list.see("end")

    console_list.insert("end", "")

    if not h.OUTPUT_DIR == '':
        if ".CSV" in h.OUTPUT_TYPES:
            h.utilities.dataframe_to_csv(df)

        if ".TXT" in h.OUTPUT_TYPES:
            h.utilities.dataframe_to_txt(df)

        if ".SDF (Only Mol Data)" in h.OUTPUT_TYPES:
            h.utilities.dataframe_to_sdf(df2)

        if ".SMI (Only SMILES)" in h.OUTPUT_TYPES:
            h.utilities.dataframe_to_smi(df)

    # h.utilities.plot_exact_mass_tpsa_a_log_p(df)


def plot_generator(x_axis, y_axis, z_axis, console_list, progress_bar):

    input_file = h.os.path.abspath(h.INPUT_FILE)

    if h.os.path.exists(input_file):
        df = h.utilities.csv_to_dataframe(input_file)
        console_list.insert("end", f"Input File: {input_file}")
        console_list.insert("end", "")
        console_list.see("end")

    else:
        h.messagebox.showinfo("Error Message", "There was a problem with your input file.")
        progress_bar.stop()

    try:
        x = df[x_axis]

    except:
        h.messagebox.showinfo("Error Message", "There was a problem with your X axis input.")
        progress_bar.stop()

    try:
        y = df[y_axis]

    except:
        h.messagebox.showinfo("Error Message", "There was a problem with your Y axis input.")
        progress_bar.stop()

    try:
        z = df[z_axis]

    except:
        h.messagebox.showinfo("Error Message", "There was a problem with your Z axis input.")
        progress_bar.stop()


    fig = h.plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, c=(x+y+z), cmap="viridis", edgecolors='black', s=20)
    ax.plot(x, z, "r+", zdir='y', zs=y.max())
    ax.plot(y, z, "g+", zdir='x', zs=x.min())
    ax.plot(x, y, "b+", zdir='z', zs=z.min())

    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    ax.set_zlabel(z_axis)

    date_time = h.datetime.datetime.now().strftime("%m%d%Y%H%M")

    output_dir = h.os.path.abspath(h.OUTPUT_DIR)

    if h.os.path.exists(output_dir):
        h.OUTPUT_DIR = h.os.path.join(h.OUTPUT_DIR, date_time + "_figures")

        temp_num = 2

        if h.os.path.exists(h.OUTPUT_DIR):
            temp_path = h.OUTPUT_DIR + "_" + str(temp_num)
            while h.os.path.exists(temp_path):
                temp_num += 1
                temp_path = h.OUTPUT_DIR + "_" + str(temp_num)
            h.OUTPUT_DIR = temp_path
            h.os.mkdir(h.OUTPUT_DIR)

        else:
            h.os.mkdir(h.OUTPUT_DIR)

        if not len(h.OUTPUT_IMG_TYPES) == 0:
            for output_type in h.OUTPUT_IMG_TYPES:
                if not output_type == '':
                    h.plt.savefig(h.os.path.join(h.OUTPUT_DIR, "3d_plot" + output_type.lower()),
                                  format=(output_type.lower()[1:]), transparent=True)

    progress_bar.stop()

    h.plt.show()


def quote_generator():
    quote = h.random.choice(h.QUOTE_LIST)

    return quote
