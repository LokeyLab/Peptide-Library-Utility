#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains all code to run the graphical
part of the GUI. See gui_functions.py for the
functions that are called.
"""

# Internal
from scripts import header as h


def create_root(size_x, size_y, offset_x, offset_y, name):
    window = h.tk.Tk()

    window.grid_columnconfigure(1, weight=1)

    window.title(name)

    main_frame = h.tk.LabelFrame(window, padx=10, pady=10)
    main_frame.pack(fill="both")

    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=1)
    main_frame.rowconfigure(0, weight=1)
    main_frame.rowconfigure(1, weight=1)

    # main_label =h.tk.Label(main_frame, text=name, width=10, pady=10)
    # main_label.grid(column=0, columnspan=1, row=0, rowspan=1, sticky='new')

    progress_bar = h.ttk.Progressbar(main_frame, orient="horizontal", length=500,
                                     mode='indeterminate')
    progress_bar.grid(column=0, columnspan=1, row=0, rowspan=1, padx=10, pady=10, sticky='w')

    # photo = ImageTk.PhotoImage(file="logo.png")
    # window.iconphoto(False, photo)

    window.resizable(False, False)

    # window.attributes('-fullscreen',True)

    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()

    window_width = size_x
    window_height = size_y

    center_x = int(screen_width / 2 - window_width / 2)
    center_y = int(screen_height / 2 - window_height / 2)

    window.geometry(f'{window_width}x{window_height}+{center_x - offset_x}+{center_y + offset_y}')

    return window, main_frame, progress_bar


def create_main_window_main_menu(root, main_frame, console_list, output_list,
                                 subunit_list, subunit_count):
    main_menu_frame = h.tk.LabelFrame(main_frame)
    main_menu_frame.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')
    main_menu_frame.columnconfigure(0, weight=1)

    main_menu_label = h.tk.Label(main_menu_frame, text="Main Menu")
    main_menu_label.grid(column=0, columnspan=2, row=0, rowspan=1, padx=10, pady=10, sticky='n')

    def pep_lib():
        peptide_library_creator(root, console_list, output_list, subunit_list, subunit_count)


    h.tk.Button(main_menu_frame, text="Peptide Library Creator", width=20).grid(column=0,
                                                                                columnspan=1,
                                                                                row=2, rowspan=1,
                                                                                padx=10, pady=10,
                                                                                sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Create a Peptide Library by Uploading a File.")
    main_menu_label.grid(column=0, columnspan=1, row=2, rowspan=1, padx=10, pady=10, sticky='e')

    """
    def pep_chem_info():
        peptide_cheminformatics(root, console_list, output_list)

    h.tk.Button(main_menu_frame, text="Peptide Cheminformatics",
                command=pep_chem_info, width=20).grid(column=0, columnspan=1,
                                    row=3, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Evaluate Cheminformatic Data For a Given SMILES String.")
    main_menu_label.grid(column=0, columnspan=1, row=3, rowspan=1, padx=10, pady=10, sticky='e')

    def pep_lib_chem_info():
        library_cheminformatics(root, console_list, output_list)

    h.tk.Button(main_menu_frame, text="Library Cheminformatics",
                command=pep_lib_chem_info, width=20).grid(column=0, columnspan=1,
                                    row=4, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Evaluate Basic Cheminformatic Data by Uploading a File.")
    main_menu_label.grid(column=0, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='e')
    """

    """
    def plot_data():
        plotting(root, console_list, output_list)
    """

    h.tk.Button(main_menu_frame, text="Plot Data", width=20).grid(column=0, columnspan=1,
                                                                  row=5, rowspan=1, padx=10,
                                                                  pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Plot Cheminformatic Data from .CSV Data File.")

    main_menu_label.grid(column=0, columnspan=1, row=5, rowspan=1, padx=10, pady=10, sticky='e')

    """
    def export_subunit_library():
        subunit_library_export(root, console_list, output_list, subunit_list,
                           subunit_count)
    """

    h.tk.Button(main_menu_frame, text="Export Subunit Library", width=20).grid(column=0,
                                                                               columnspan=1,
                                                                               row=6, rowspan=1,
                                                                               padx=10, pady=10,
                                                                               sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Export Subunit Library to .CSV File.")

    main_menu_label.grid(column=0, columnspan=1, row=6, rowspan=1, padx=10,
                         pady=10, sticky='e')

    def drive_button():
        h.webbrowser.open(
            'https://drive.google.com/drive/folders/1RCPQ5P0D3kjbiswlv-xwSL3Vja2ngBNP')
        console_list.insert("end", "Lokey Lab Google Drive Opened.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(main_menu_frame, text="Subunit Drive",
                command=drive_button, width=20).grid(column=0, columnspan=1,
                                                     row=7, rowspan=1, padx=10, pady=10,
                                                     sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Navigate to Lokey Lab Google "
                                      "Drive Containing Subunit Files.")
    main_menu_label.grid(column=0, columnspan=1, row=7, rowspan=1, padx=10,
                         pady=10, sticky='e')

    def help_button():
        h.webbrowser.open('https://github.com/LokeyLab/Peptide-Library-Utility')
        console_list.insert("end", "Peptide Library Utility Wiki Opened.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(main_menu_frame, text="Help",
                command=help_button, width=20).grid(column=0, columnspan=1,
                                                    row=8, rowspan=1, padx=10, pady=10,
                                                    sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Navigate to the Peptide Library Utility Wiki.")
    main_menu_label.grid(column=0, columnspan=1, row=8, rowspan=1, padx=10,
                         pady=10, sticky='e')

    def quote_button():
        quote = h.gui_functions.quote_generator()
        output_list.insert("end", quote)
        output_list.insert("end", "")
        output_list.see("end")

    def font_button():
        output_list.config(font=h.font.Font(family="Comic Sans MS"), height=11)
        console_list.config(font=h.font.Font(family="Comic Sans MS"), height=11)
        subunit_list.config(font=h.font.Font(family="Comic Sans MS"), height=14)
        main_menu_frame.config(height=14)

    quote_button1 = h.tk.Button(root, text="", command=quote_button)
    quote_button1.place(x=-34, y=702, height=50, width=50)
    quote_button2 = h.tk.Button(root, text="", command=font_button)
    quote_button2.place(x=1062, y=702, height=50, width=50)

    return main_menu_frame


def create_main_window_subunits(main_frame, console_list, main_progress_bar):
    subunit_frame = h.tk.LabelFrame(main_frame)
    subunit_frame.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nsew')
    subunit_frame.columnconfigure(0, weight=1)

    subunit_vertical_scrollbar = h.tk.Scrollbar(subunit_frame)
    subunit_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nse')
    subunit_horizontal_scrollbar = h.tk.Scrollbar(subunit_frame, orient='horizontal')
    subunit_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='sew')

    subunit_list = h.tk.Listbox(subunit_frame, yscrollcommand=subunit_vertical_scrollbar.set,
                                xscrollcommand=subunit_horizontal_scrollbar.set, height=16)

    subunit_list.insert("end", "None")

    h.tk.Label(subunit_frame, text="Available Subunits").grid(column=0, columnspan=2,
                                                              row=0, rowspan=1,
                                                              padx=10, pady=10,
                                                              sticky='')

    subunit_list.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')

    subunit_vertical_scrollbar.config(command=subunit_list.yview)
    subunit_horizontal_scrollbar.config(command=subunit_list.xview)

    subunit_count = h.tk.StringVar()

    subunit_count.set("n = 0")

    def lib_update():
        main_progress_bar.start(10)
        h.gui_functions.subunit_library_update(console_list, subunit_list)
        subunit_count.set("n = " + str(h.SUBUNIT_COUNT))
        main_progress_bar.stop()
        h.SUBUNIT_LIBRARY_POPULATED = True

    populate_subunit_library_button = h.tk.Button(main_frame,
                                                  text="Populate Subunit Library",
                                                  command=lib_update)

    populate_subunit_library_button.grid(column=1, columnspan=1, row=0, rowspan=1, pady=10,
                                         sticky="nw")

    subunit_label = h.tk.Label(subunit_frame, textvariable=subunit_count, width=20)
    subunit_label.grid(column=0, columnspan=2, row=3, rowspan=1, padx=10, pady=10, sticky="w")

    def subunit_scan_key(event):
        subunit_list.delete(0, 'end')
        h.gui_functions.subunit_scan_key(event)

        if len(h.SUBUNIT_LIST) == 0:
            subunit_list.insert('end', "None")

        else:
            for item in h.SUBUNIT_LIST:
                subunit_list.insert('end', item)

        subunit_count.set(h.SUBUNIT_COUNT)

    entry = h.tk.Entry(subunit_frame)
    entry.grid(column=0, columnspan=1,
               row=3, rowspan=1,
               padx=10, pady=10,
               sticky='se')

    entry.bind('<KeyRelease>', subunit_scan_key)

    return subunit_frame, subunit_list, subunit_count


def create_main_window_output(main_frame, console_list):
    peptide_frame = h.tk.LabelFrame(main_frame)
    peptide_frame.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='nesw')
    peptide_frame.columnconfigure(0, weight=1)

    peptide_vertical_scrollbar = h.tk.Scrollbar(peptide_frame)
    peptide_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nse')
    peptide_horizontal_scrollbar = h.tk.Scrollbar(peptide_frame, orient='horizontal')
    peptide_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='ews')

    output_list = h.tk.Listbox(peptide_frame, yscrollcommand=peptide_vertical_scrollbar.set,
                               xscrollcommand=peptide_horizontal_scrollbar.set, height=12)

    h.tk.Label(peptide_frame, text="Output").grid(column=0, columnspan=2,
                                                  row=0, rowspan=1,
                                                  padx=10, pady=10,
                                                  sticky='')

    output_list.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')

    peptide_vertical_scrollbar.config(command=output_list.yview)
    peptide_horizontal_scrollbar.config(command=output_list.xview)

    return peptide_frame, output_list


def create_main_window_console(main_frame):
    console_frame = h.tk.LabelFrame(main_frame)
    console_frame.grid(column=1, columnspan=1, row=2, rowspan=1, sticky='nesw')
    console_frame.columnconfigure(0, weight=1)

    console_vertical_scrollbar = h.tk.Scrollbar(console_frame)
    console_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='ns')
    console_horizontal_scrollbar = h.tk.Scrollbar(console_frame, orient='horizontal')
    console_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='ew')

    console_list = h.tk.Listbox(console_frame, yscrollcommand=console_vertical_scrollbar.set,
                                xscrollcommand=console_horizontal_scrollbar.set, height=12)

    h.tk.Label(console_frame, text="Console").grid(column=0, columnspan=2,
                                                   row=0, rowspan=1,
                                                   padx=10, pady=10,
                                                   sticky='')

    console_list.grid(column=0, columnspan=1,
                      row=1, rowspan=1,
                      sticky='nsew')

    console_vertical_scrollbar.config(command=console_list.yview)
    console_horizontal_scrollbar.config(command=console_list.xview)

    return console_frame, console_list


def create_main_window_subunits(main_frame, console_list, main_progress_bar):
    subunit_frame = h.tk.LabelFrame(main_frame)
    subunit_frame.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nsew')
    subunit_frame.columnconfigure(0, weight=1)

    subunit_vertical_scrollbar = h.tk.Scrollbar(subunit_frame)
    subunit_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nse')
    subunit_horizontal_scrollbar = h.tk.Scrollbar(subunit_frame, orient='horizontal')
    subunit_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='sew')

    subunit_list = h.tk.Listbox(subunit_frame, yscrollcommand=subunit_vertical_scrollbar.set,
                                xscrollcommand=subunit_horizontal_scrollbar.set, height=16)

    subunit_list.insert("end", "None")

    h.tk.Label(subunit_frame, text="Available Subunits").grid(column=0, columnspan=2,
                                                              row=0, rowspan=1,
                                                              padx=10, pady=10,
                                                              sticky='')

    subunit_list.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')

    subunit_vertical_scrollbar.config(command=subunit_list.yview)
    subunit_horizontal_scrollbar.config(command=subunit_list.xview)

    subunit_count = h.tk.StringVar()

    subunit_count.set("n = 0")

    def lib_update():
        main_progress_bar.start(10)
        h.gui_functions.subunit_library_update(console_list, subunit_list)
        subunit_count.set("n = " + str(h.SUBUNIT_COUNT))
        main_progress_bar.stop()
        h.SUBUNIT_LIBRARY_POPULATED = True

    populate_subunit_library_button = h.tk.Button(main_frame,
                                                  text="Populate Subunit Library",
                                                  command=lib_update)

    populate_subunit_library_button.grid(column=1, columnspan=1, row=0, rowspan=1, pady=10,
                                         sticky="nw")

    subunit_label = h.tk.Label(subunit_frame, textvariable=subunit_count, width=20)
    subunit_label.grid(column=0, columnspan=2, row=3, rowspan=1, padx=10, pady=10, sticky="w")

    def subunit_scan_key(event):
        subunit_list.delete(0, 'end')
        h.gui_functions.subunit_scan_key(event)

        if len(h.SUBUNIT_LIST) == 0:
            subunit_list.insert('end', "None")

        else:
            for item in h.SUBUNIT_LIST:
                subunit_list.insert('end', item)

        subunit_count.set(h.SUBUNIT_COUNT)

    entry = h.tk.Entry(subunit_frame)
    entry.grid(column=0, columnspan=1,
               row=3, rowspan=1,
               padx=10, pady=10,
               sticky='se')

    entry.bind('<KeyRelease>', subunit_scan_key)

    return subunit_frame, subunit_list, subunit_count


def create_main_window_output(main_frame, console_list):
    peptide_frame = h.tk.LabelFrame(main_frame)
    peptide_frame.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='nesw')
    peptide_frame.columnconfigure(0, weight=1)

    peptide_vertical_scrollbar = h.tk.Scrollbar(peptide_frame)
    peptide_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nse')
    peptide_horizontal_scrollbar = h.tk.Scrollbar(peptide_frame, orient='horizontal')
    peptide_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='ews')

    output_list = h.tk.Listbox(peptide_frame, yscrollcommand=peptide_vertical_scrollbar.set,
                               xscrollcommand=peptide_horizontal_scrollbar.set, height=12)

    h.tk.Label(peptide_frame, text="Output").grid(column=0, columnspan=2,
                                                  row=0, rowspan=1,
                                                  padx=10, pady=10,
                                                  sticky='')

    output_list.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')

    peptide_vertical_scrollbar.config(command=output_list.yview)
    peptide_horizontal_scrollbar.config(command=output_list.xview)

    return peptide_frame, output_list


def create_main_window_console(main_frame):
    console_frame = h.tk.LabelFrame(main_frame)
    console_frame.grid(column=1, columnspan=1, row=2, rowspan=1, sticky='nesw')
    console_frame.columnconfigure(0, weight=1)

    console_vertical_scrollbar = h.tk.Scrollbar(console_frame)
    console_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='ns')
    console_horizontal_scrollbar = h.tk.Scrollbar(console_frame, orient='horizontal')
    console_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='ew')

    console_list = h.tk.Listbox(console_frame, yscrollcommand=console_vertical_scrollbar.set,
                                xscrollcommand=console_horizontal_scrollbar.set, height=12)

    h.tk.Label(console_frame, text="Console").grid(column=0, columnspan=2,
                                                   row=0, rowspan=1,
                                                   padx=10, pady=10,
                                                   sticky='')

    console_list.grid(column=0, columnspan=1,
                      row=1, rowspan=1,
                      sticky='nsew')

    console_vertical_scrollbar.config(command=console_list.yview)
    console_horizontal_scrollbar.config(command=console_list.xview)

    return console_frame, console_list


def create_main_window_subunits(main_frame, console_list, main_progress_bar):
    subunit_frame = h.tk.LabelFrame(main_frame)
    subunit_frame.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nsew')
    subunit_frame.columnconfigure(0, weight=1)

    subunit_vertical_scrollbar = h.tk.Scrollbar(subunit_frame)
    subunit_vertical_scrollbar.grid(column=1, columnspan=1, row=1, rowspan=1, sticky='nse')
    subunit_horizontal_scrollbar = h.tk.Scrollbar(subunit_frame, orient='horizontal')
    subunit_horizontal_scrollbar.grid(column=0, columnspan=1, row=2, rowspan=1, sticky='sew')

    subunit_list = h.tk.Listbox(subunit_frame, yscrollcommand=subunit_vertical_scrollbar.set,
                                xscrollcommand=subunit_horizontal_scrollbar.set, height=16)

    subunit_list.insert("end", "None")

    h.tk.Label(subunit_frame, text="Available Subunits").grid(column=0, columnspan=2,
                                                              row=0, rowspan=1,
                                                              padx=10, pady=10,
                                                              sticky='')

    subunit_list.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')

    subunit_vertical_scrollbar.config(command=subunit_list.yview)
    subunit_horizontal_scrollbar.config(command=subunit_list.xview)

    subunit_count = h.tk.StringVar()

    subunit_count.set("n = 0")

    def lib_update():
        main_progress_bar.start(10)
        count = h.gui_functions.subunit_library_update(console_list, subunit_list)
        main_progress_bar.stop()
        h.SUBUNIT_LIBRARY_POPULATED = True
        subunit_count.set("n = " + str(count))

    populate_subunit_library_button = h.tk.Button(main_frame,
                                                  text="Populate Subunit Library",
                                                  command=lib_update)

    populate_subunit_library_button.grid(column=1, columnspan=1, row=0, rowspan=1, pady=10,
                                         sticky="nw")

    subunit_label = h.tk.Label(subunit_frame, textvariable=subunit_count)
    subunit_label.grid(column=0, columnspan=1, row=3, rowspan=1, padx=10, pady=10, sticky="ws")

    def subunit_scan_key(event):
        subunit_list.delete(0, 'end')
        h.gui_functions.subunit_scan_key(event)

        if len(h.SUBUNIT_LIST) == 0:
            subunit_list.insert('end', "None")

        else:
            for item in h.SUBUNIT_LIST:
                subunit_list.insert('end', item)

        subunit_count.set(h.SUBUNIT_COUNT)

    entry = h.tk.Entry(subunit_frame, width=10)

    entry.grid(column=0, columnspan=2,
               row=3, rowspan=1,
               padx=10, pady=10,
               sticky='e')

    entry.bind('<KeyRelease>', subunit_scan_key)

    h.tk.Label(subunit_frame, text="Current Subunit").grid(column=2, columnspan=1,
                                                           row=0, rowspan=1,
                                                           padx=10, pady=10,
                                                           sticky='')

    subunit_info_frame = h.tk.LabelFrame(subunit_frame)
    subunit_info_frame.grid(column=2, columnspan=1,
                            row=1, rowspan=4,
                            sticky='nesw')

    fig = h.plt.figure(figsize=(1.5, 1.5))
    ax = fig.add_subplot()
    h.plt.axis('off')
    img = h.FigureCanvasTkAgg(fig, subunit_info_frame)
    img.get_tk_widget().grid(column=0, columnspan=1,
                             row=0, rowspan=1, sticky="new")

    exact_mass = h.tk.StringVar()
    a_log_p_86 = h.tk.StringVar()
    a_log_p = h.tk.StringVar()
    tpsa = h.tk.StringVar()
    hba_hbd = h.tk.StringVar()

    exact_mass.set("Exact Mass = ")
    a_log_p_86.set("AlogP 86 = ")
    a_log_p.set("AlogP 99 = ")
    tpsa.set("TPSA = ")
    hba_hbd.set("HBA : HBD = ")


    h.tk.Label(subunit_info_frame, textvariable=exact_mass).grid(column=0, columnspan=1,
                                                                 row=1, rowspan=1,
                                                                 padx=0, pady=0,
                                                                 sticky='w')

    h.tk.Label(subunit_info_frame, textvariable=a_log_p_86).grid(column=0, columnspan=1,
                                                              row=2, rowspan=1,
                                                              padx=0, pady=0,
                                                              sticky='w')

    h.tk.Label(subunit_info_frame, textvariable=a_log_p).grid(column=0, columnspan=1,
                                                              row=3, rowspan=1,
                                                              padx=0, pady=0,
                                                              sticky='w')

    h.tk.Label(subunit_info_frame, textvariable=tpsa).grid(column=0, columnspan=1,
                                                           row=4, rowspan=1,
                                                           padx=0, pady=0,
                                                           sticky='w')

    h.tk.Label(subunit_info_frame, textvariable=hba_hbd).grid(column=0, columnspan=1,
                                                           row=5, rowspan=1,
                                                           padx=0, pady=0,
                                                           sticky='w')

    def onselect(evt):
        w = evt.widget
        index = int(w.curselection()[0])
        value = w.get(index)
        print('You selected item %d: "%s"' % (index, value))
        key = str(value).split()[-1]

        smiles = h.SUBUNIT_LIBRARY[key].smiles_string
        mol = h.cheminformatics.get_molecule_from_smiles(smiles)
        molecule = h.RemoveHs(mol)
        molecule_img = h.Draw.MolToImage(molecule, size=(256, 256))
        ax.imshow(molecule_img)
        img.draw()

        exact_mass.set("Exact Mass = %.4f" % h.cheminformatics.get_exact_mass(mol))
        a_log_p_86.set("AlogP 86 COMING SOON") # % h.a_log_p_86.get_a_log_p_86(smiles))
        # a_log_p_86.set("AlogP 86 = %.4f" % h.a_log_p_86.get_a_log_p_86(smiles))
        a_log_p.set("AlogP 99 = %.4f" % h.cheminformatics.get_a_log_p_99(mol))
        tpsa.set("TPSA = %.4f" % h.cheminformatics.get_tpsa(mol))
        hba_hbd.set(f"HBA : HBD = {h.cheminformatics.get_num_hba(mol)} : {h.cheminformatics.get_num_hbd(mol)}")

    subunit_list.bind('<<ListboxSelect>>', onselect)

    return subunit_frame, subunit_list, subunit_count


def populate_root(root, main_frame, main_progress_bar):
    h.tk.Button(main_frame, text="Exit",
                command=root.destroy).grid(column=1, columnspan=1,
                                           row=0, rowspan=1,
                                           pady=10, sticky="ne")

    console_frame, console_list = create_main_window_console(main_frame)
    output_frame, output_list = create_main_window_output(main_frame, console_list)
    subunit_frame, subunit_list, subunit_count = create_main_window_subunits(main_frame,
                                                                             console_list,
                                                                             main_progress_bar)
    create_main_window_main_menu(root, main_frame, console_list, output_list,
                                 subunit_list, subunit_count)


def main():
    """main function of gui2 where all other gui related functions are called."""

    root, main_frame, main_progress_bar = create_root(1080, 720, 0, 0, "Peptide Library Utility")

    populate_root(root, main_frame, main_progress_bar)

    root.mainloop()
