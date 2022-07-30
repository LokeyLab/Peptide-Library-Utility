#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains all the user interface related
functions and most of the function calls that are
related to the creation of peptide libraries.
"""

import header as h


def cyclic_peptide_in_terminal(subunit_library):
    """This function guides the user through
    creating a peptide using the text-based user
    interface."""

    subunit_count = 0

    while True:

        try:

            subunit_count = int(input("\nHow many subunits?: "))

            if subunit_count > 0:

                break

            else:

                print("\nInput was less than 1. Try again.")

                continue

        except ValueError:

            print("\nInput was not a number. Try again.")

            continue

    subunits = []

    print("\nExample input: (NMe)-ß-H-D-Leu or Leu")

    for i in range(0, subunit_count):

        while True:

            try:

                subunit = subunit_library[input("\nEnter subunit: ")]

                break

            except KeyError:

                print("\nInvalid subunit. Check the subunit list and try again.")

                continue

        print("\nSubunit " + str(i + 1) + ": " + subunit.name)

        subunits.append(subunit)

    print("")

    name = ""

    for i in range(0, len(subunits)):
        name += subunits[i].multiple_letter + " "

    peptide = h.bonding.bond_head_to_tail(subunits)

    print("Cyclization Types:")
    print("1) Head to Tail")
    print("2) Huisgen FIX ME")
    print("3) Macrolactone (Must end in Threonine)")
    print("4) Keep Linear")

    while True:

        try:

            choice = int(input("\nEnter a number: "))

            if choice < 1 or choice > 4:

                print("\nInput was not a valid choice. Try again.")

                continue

            else:

                break

        except ValueError:

            print("\nInput was not a number. Try again.")

            continue

    if choice == 1:

        print("\nCyclizing Head to Tail...")

        cyclic_peptide = h.bonding.cyclize_head_to_tail(peptide)

        print("\nCyclization Complete.")

    elif choice == 2:

        print("\nCyclizing Huisgen...")

        cyclic_peptide = h.bonding.cyclize_huisgen(peptide)

        print("\nCyclization Complete.")

    elif choice == 3:

        print("\nCyclizing Macrolactone...")

        cyclic_peptide = h.bonding.cyclize_Macrolactone(peptide, subunit_library)

        print("\nCyclization Complete.")

        while True:

            try:

                choice = input("\nCap exposed Amino Group? y/n: ")

                if not choice == "y" or choice == "n":

                    print("\nInput was not a valid choice. Try again.")

                    continue

                else:

                    break

            except ValueError:

                print("\nInput was not either \"y\" or \"n\". Try again.")

                continue

        if choice == "y":

            while True:

                try:

                    subunit = subunit_library[input("\nEnter subunit: ")]

                    break

                except KeyError:

                    print("\nInvalid subunit. Check the subunit list and try again.")

                    continue

            temp_peptide = str(cyclic_peptide[0:cyclic_peptide.rfind("N") + 1]) + "(" + \
                               subunit.smiles_string + ")" + \
                               str(cyclic_peptide[cyclic_peptide.rfind("N") + 1:])

            cyclic_peptide = temp_peptide

        elif choice == "n":

            pass

    elif choice == 4:

        cyclic_peptide = peptide

    else:

        print("\nSomething went wrong.")

    while True:

        try:

            molecule = h.cheminformatics.get_molecule_from_smiles(
                cyclic_peptide)

            exact_mass, tpsa, a_log_p = h.cheminformatics.get_chemometrics(
                molecule)

            cyclic_peptide_object = h.classes.Peptide(name, exact_mass, tpsa,
                                                      a_log_p, cyclic_peptide)

            df = h.utilities.peptides_to_dataframe([cyclic_peptide_object])

            h.utilities.print_dataframe(df)

            h.utilities.dataframe_to_csv(df)

            break

        except TypeError:

            print("\n\nInput was not a valid SMILES string. Try again.")

            continue

    return cyclic_peptide


def peptide_library_from_csv(subunit_library):
    """Creates a peptide library from a user provided .csv."""

    while True:

        try:

            file = h.os.path.normcase("input/" + input("\nEnter file name: ") + ".csv")

            print("\n" + str(file))

            df = h.utilities.csv_to_dataframe(file)

            break

        except FileNotFoundError:

            print("\nFile not found. Try again.")

            continue

    h.utilities.print_dataframe(df)

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

    cyclic_peptide_count = pot_sizes[0]

    for k in range(1, len(pot_sizes)):
        cyclic_peptide_count *= pot_sizes[k]

    print("\nCyclic Peptides To Be Generated: " + str(cyclic_peptide_count))


    cartesian_product = h.combinatronics.cartesian_product(pots)

    print("\nCyclization Types:")
    print("1) Head to Tail")
    print("2) Huisgen")
    print("3) Macrolactone (Must end in Threonine)")
    print("4) Keep Linear")

    while True:

        try:

            choice = int(input("\nEnter a number: "))

            if choice < 1 or choice > 4:

                print("\nInput was not a valid choice. Try again.")

                continue

            else:

                break

        except ValueError:

            print("\nInput was not a number. Try again.")

            continue

    if choice == 1:

        print("\nCyclizing Head to Tail...")

        cyclization_type = "Head to Tail"

    elif choice == 2:

        print("\nCyclizing Huisgen...")

        cyclization_type = "Huisgen"
        

        while True:

            try:

                choice = input("\nCap exposed carboxyl group? y/n: ")

                if not (choice == "y" or choice == "n"):

                    print("\nInput was not a valid choice. Try again.")

                    continue

                else:

                    break

            except ValueError:

                print("\nInput was not either \"y\" or \"n\". Try again.")

                continue

        while True:

            if choice == "y":

                while True:

                    try:

                        cap = subunit_library[input("\nEnter subunit: ")]

                        print("\n" + str(cap.name) + " cap")

                        cap_bool = True

                        break

                    except KeyError:

                        print("\nInvalid subunit. Check the subunit list and try again.")

                        continue

                break

            elif choice == "n":

                cap_bool = False

                break

            else:

                print("\nSomething went wrong.")

                continue

    elif choice == 3:

        print("\nCyclizing Macrolactone...")

        cyclization_type = "Macrolactone"

        while True:

            try:

                choice = input("\nCap exposed Amino Group? y/n: ")

                if not (choice == "y" or choice == "n"):

                    print("\nInput was not a valid choice. Try again.")

                    continue

                else:

                    break

            except ValueError:

                print("\nInput was not either \"y\" or \"n\". Try again.")

                continue

        while True:

            if choice == "y":

                while True:

                    try:

                        cap = subunit_library[input("\nEnter subunit: ")]

                        print("\n" + str(cap.name) + " cap")

                        cap_bool = True

                        break

                    except KeyError:

                        print("\nInvalid subunit. Check the subunit list and try again.")

                        continue

                break

            elif choice == "n":

                cap_bool = False

                break

            else:

                print("\nSomething went wrong.")

                continue

    elif choice == 4:

        cyclization_type = "None"

    else:

        print("\nSomething went wrong.")

    cyclic_peptides = []

    molecules_list = []

    count = 0
    print("")

    for i in h.tqdm.tqdm(cartesian_product):

        subunits = []

        for j in i:

            try:
                subunit = subunit_library[j]

                subunits.append(subunit)

            except KeyError:

                print(f"\nSubunit \"%s\" is not a valid subunit." % j)
                print("Change the input file and try again.")

                exit()

        peptide = h.bonding.bond_head_to_tail(subunits)

        if cyclization_type == "Head to Tail":

            cyclic_peptide = h.bonding.cyclize_head_to_tail(peptide)

        elif cyclization_type == "Huisgen":

            cyclic_peptide = h.bonding.cyclize_triazole(peptide, subunit_library)

            subunits.append(cap)

            tempPeptide = str(cyclic_peptide[0:-1]) + cap.smiles_string

            cyclic_peptide = tempPeptide

        elif cyclization_type == "Macrolactone":

            cyclic_peptide = h.bonding.cyclize_Macrolactone(peptide, subunit_library)

            if cap_bool:
                subunits.append(cap)

                tempPeptide = str(cyclic_peptide[0:cyclic_peptide.rfind("N") + 1]) + "(" + \
                              cap.smiles_string + ")" + \
                              str(cyclic_peptide[cyclic_peptide.rfind("N") + 1:])

                cyclic_peptide = tempPeptide

        elif cyclization_type == "None":

            cyclic_peptide = peptide

        name = ""


        for i in range(0, len(subunits)):
            name += subunits[i].multiple_letter + " "

        molecule = h.cheminformatics.get_molecule_from_smiles(cyclic_peptide)

        exactMass, tpsa, aLogP = h.cheminformatics.get_chemometrics(molecule)

        cyclic_peptide_object = h.classes.Peptide(name, exactMass, tpsa, aLogP, cyclic_peptide)

        cyclic_peptides.append(cyclic_peptide_object)

        molecules_list.append(molecule)

    if not len(cyclic_peptides) == cyclic_peptide_count:

        print("\nSomething was wrong with one or more of your subunits.")
        print("Check your subunit input file and try again.")

        exit()

    df = h.utilities.peptides_to_dataframe(cyclic_peptides)

    df2 = h.utilities.molecules_to_dataframe(molecules_list)

    h.utilities.print_dataframe(df)

    h.utilities.dataframe_to_csv(df)

    # h.utilities.dataframe_to_sdf(df2)

    # h.utilities.plot_exact_mass_tpsa_a_log_p(df)

    return


def Introduction():
    print("0) Populate Subunit Library")
    print("1) Create Cyclic Peptide in Terminal FIX ME (add Huisgen)")
    print("2) Create Cyclic Peptide Library from .CSV FIX ME (add Huisgen)")
    print("3) Calculate Chemometrics of a SMILES string")
    print("4) Check a Subunit")
    print("5) Quit")

    while True:

        try:

            choice = int(input("\nEnter number: "))

            if choice < 0 or choice > 5:

                print("\nInput was not a valid choice. Try again.")

                continue

            else:

                break

        except ValueError:

            print("\nInput was not a number. Try again.")

            continue

    return choice


def gui_loop():
    h.root.mainloop()
    h.subunit_window.mainloop()
    h.peptide_window.mainloop()


def ui_loop(subunitLibrary):
    print("\n-------------------- Main Menu --------------------")

    choice = Introduction()

    if choice == 0:

        print("\nPopulating Subunit Library...\n")

        subunitLibrary = h.subunit_builder.generate_subunit_library()

        print("\nSubunit Library Populated.")

    elif choice == 1:
        print("\nCreating Cyclic Peptide...")

        cyclic_peptide_in_terminal(subunitLibrary)

        print("\nCyclic Peptide Created.")

    elif choice == 2:
        print("\nCreating Cyclic Peptide Library...")

        peptide_library_from_csv(subunitLibrary)

        print("\nCyclic Peptide Library Created.")

    elif choice == 3:

        while True:

            try:

                cyclic_peptide = input("\nEnter SMILES String: ")

                molecule = h.cheminformatics.get_molecule_from_smiles(cyclic_peptide)

                exact_mass, tpsa, a_log_p = h.cheminformatics.get_chemometrics(molecule)

                cyclic_peptide_object = h.classes.Peptide("", exact_mass, tpsa, a_log_p,
                                                          cyclic_peptide)

                df = h.utilities.peptides_to_dataframe([cyclic_peptide_object])

                h.utilities.print_dataframe(df)

                break

            except TypeError:

                print("\n\nInput was not a valid SMILES string. Try again.")

                continue

    elif choice == 4:

        while True:

            try:
                subunit = subunitLibrary[input("\nEnter subunit: ")]

                print("")

                molecule = h.cheminformatics.get_molecule_from_smiles(subunit.smiles_string)

                exact_mass, tpsa, a_log_p = h.cheminformatics.get_chemometrics(molecule)

                columns = ["Name", "Exact Mass", "TPSA", "ALogP", "SMILES String"]

                df = h.pd.DataFrame(columns=columns)

                df.loc[len(df.index)] = {"Name": subunit.multiple_letter,
                                         "Exact Mass": exact_mass,
                                         "TPSA": tpsa,
                                         "ALogP": a_log_p,
                                         "SMILES String": subunit.smiles_string}

                h.utilities.print_dataframe(df)

            except KeyError:

                print("\nInvalid subunit. Check the subunit list and try again.")

                continue

    elif choice == 5:

        print("\nGoodbye!")

        exit()

    ui_loop(subunitLibrary)
