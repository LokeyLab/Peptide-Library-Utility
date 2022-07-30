#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script contains all code to run the graphical
part of the GUI. See gui_functions.py for the
functions that are called.
"""

# Internal
import header as h


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
# DONE


def create_window(root, size_x, size_y, offset_x, offset_y, name):

    window = h.tk.Toplevel(root)

    window.grid_columnconfigure(1, weight=1)

    window.title(name)

    main_frame = h.tk.LabelFrame(window, padx=10, pady=10)
    main_frame.pack(fill="both")

    main_frame.grid_columnconfigure(0, weight=1)
    main_frame.grid_columnconfigure(1, weight=1)

    # main_label =h.tk.Label(main_frame, text=name, width=10, pady=10)
    # main_label.grid(column=0, columnspan=1, row=0, rowspan=1, sticky='new')

    progress_bar = h.ttk.Progressbar(main_frame, orient="horizontal", mode='determinate')
    progress_bar.grid(column=0, columnspan=1, row=0, rowspan=1, padx=10, pady=10, sticky='nsew')

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
# DONE


def single_peptide_creator(root, console_list):
    console_list.insert("end", "Peptide Creator Started.")
    console_list.insert("end", "")
    console_list.see("end")

    single_peptide_window, \
        single_peptide_main_frame, \
        single_peptide_progress_bar = create_window(root, 1080, 720, 10, 10, "Peptide Creator")

    def single_peptide_exit_button():
        single_peptide_window.destroy()
        console_list.insert("end", "Peptide Creator Exited.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(single_peptide_main_frame, text="Exit",
                command=single_peptide_exit_button).grid(column=1, columnspan=1,
                                                         row=0, rowspan=1,
                                                         pady=10,
                                                         sticky="ne")

    peptide = []

    h.tk.Label(single_peptide_main_frame, text="Subunit Count").grid(column=0, columnspan=1,
                                                                     row=1, rowspan=1,
                                                                     padx=10, pady=10,
                                                                     sticky='w')

    cap_var = h.tk.StringVar()
    cap_var.set("Cap Subunit")
    h.tk.Label(single_peptide_main_frame, textvariable=cap_var).grid(column=1, columnspan=1,
                                                                     row=1, rowspan=1,
                                                                     padx=10, pady=10,
                                                                     sticky='w')

    subunit_count = h.tk.Entry(single_peptide_main_frame)
    subunit_count.grid(column=1, columnspan=1, row=1, rowspan=1, pady=10, sticky='e')
    # subunit_count.bind('<Return>', Scankey)

    subunit_cyclization_frame = h.tk.LabelFrame(single_peptide_main_frame)

    h.tk.Label(subunit_cyclization_frame, text="Cyclization Type").grid(column=0, columnspan=1,
                                                                     row=3, rowspan=1,
                                                                     padx=10, pady=10,
                                                                     sticky='w')

    subunit_cyclization_frame.grid(column=0, columnspan=2, row=3, rowspan=1, sticky='new')

    subunit_cyclization_frame.grid_columnconfigure(0, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(1, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(2, weight=1)

    radio_var = h.tk.IntVar()

    cyclization_types = ["Keep Linear", "Head to Tail", "Depsipeptide",  "Thioether", "Triazole"]
    R0 = h.tk.Radiobutton(subunit_cyclization_frame, text="Keep Linear", variable=radio_var, value=0)
    R0.grid(column=0, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='new')
    R1 = h.tk.Radiobutton(subunit_cyclization_frame, text="Head to Tail", variable=radio_var, value=1)
    R1.grid(column=1, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='new')
    R2 = h.tk.Radiobutton(subunit_cyclization_frame, text="Depsipeptide", variable=radio_var, value=2)
    R2.grid(column=2, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='new')
    R3 = h.tk.Radiobutton(subunit_cyclization_frame, text="Thioether", variable=radio_var, value=3)
    R3.grid(column=3, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='new')
    R4 = h.tk.Radiobutton(subunit_cyclization_frame, text="Triazole", variable=radio_var, value=4)
    R4.grid(column=4, columnspan=1, row=4, rowspan=1, padx=10, pady=10, sticky='new')

    subunit_positions_frame = h.tk.LabelFrame(single_peptide_main_frame)

    subunit_positions_frame.grid(column=0, columnspan=2, row=4, rowspan=1, sticky='ew')

    subunit_positions_frame.grid_columnconfigure(0, weight=1)
    subunit_positions_frame.grid_columnconfigure(1, weight=1)
    subunit_positions_frame.grid_columnconfigure(2, weight=1)
    subunit_positions_frame.grid_columnconfigure(3, weight=1)
    subunit_positions_frame.grid_columnconfigure(4, weight=1)
    subunit_positions_frame.grid_columnconfigure(5, weight=1)

    def Scankey(event):

        val = event.widget.get()

        if val == '':
            Update(0)

        elif not val == '':
            Update(int(val))

    def Update(val):

        for widget in subunit_positions_frame.winfo_children():
            widget.destroy()

        col = 0
        row = 2

        subunit_entries = []

        for i in range(0, val):
            h.tk.Label(subunit_positions_frame, text=f"Position {i + 1}").grid(column=col, columnspan=1,
                                                                              row=row, rowspan=1,
                                                                              padx=0, pady=5,
                                                                              sticky='new')

            subunit_entry = h.tk.Entry(subunit_positions_frame)
            subunit_entry.grid(column=col, columnspan=1, row=row + 1, rowspan=1, padx=0, pady=5, sticky='new')
            subunit_entries.append(subunit_entry)

            col += 1

            if col == 6:
                col = 0
                row += 2

    subunit_count = h.tk.Entry(single_peptide_main_frame)
    subunit_count.grid(column=0, columnspan=1, row=1, rowspan=1, pady=10, sticky='e')
    subunit_count.bind('<KeyRelease>', Scankey)

    peptide_create_button = h.tk.Button(single_peptide_main_frame, text="Create Peptide",
                                        command=h.gui_functions.single_peptide_create_button(single_peptide_window,
                                                                                             console_list,
                                                                                             h.SUBUNIT_LIBRARY,
                                                                                             h.SUBUNIT_LIST,
                                                                                             single_peptide_progress_bar))

    peptide_create_button.grid(column=0, columnspan=1, row=0, rowspan=1, sticky="e")

    return peptide


def peptide_library_creator(root, console_list, output_list, subunit_list):
    console_list.insert("end", "Peptide Library Creator Started.")
    console_list.insert("end", "")
    console_list.see("end")

    peptide_library_window, \
        peptide_library_main_frame, \
        peptide_library_progress_bar = create_window(root, 1080, 620, 10, -40,
                                                     "Peptide Library Creator")

    peptide_library_main_frame.grid_columnconfigure(0, weight=1)
    peptide_library_main_frame.grid_columnconfigure(1, weight=1)

    def peptide_library_exit_button():
        peptide_library_window.destroy()
        console_list.insert("end", "Peptide Library Creator Exited.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(peptide_library_main_frame, text="Exit",
                command=peptide_library_exit_button).grid(column=1, columnspan=1, row=0, rowspan=1,
                                                          pady=10, sticky="ne")

    subunit_cap_frame = h.tk.LabelFrame(peptide_library_main_frame)
    subunit_cap_frame.grid(column=0, columnspan=2, row=2, rowspan=1, sticky='nsew')

    h.tk.Label(subunit_cap_frame, text="Subunit Capping").grid(column=0, columnspan=1,
                                                                         row=0, rowspan=1,
                                                                         padx=10, pady=10,
                                                                         sticky='w')

    subunit_cap_frame.grid_columnconfigure(0, weight=1)
    subunit_cap_frame.grid_columnconfigure(1, weight=1)

    h.tk.Label(subunit_cap_frame, text="Amine Cap").grid(column=0, columnspan=1,
                                                                     row=1, rowspan=1,
                                                                     padx=10, pady=10,
                                                                     sticky="w")

    amine_cap = h.tk.Entry(subunit_cap_frame, width=20)
    amine_cap.grid(column=0, columnspan=1, row=1, rowspan=1, padx=10, pady=10, sticky='e')

    h.tk.Label(subunit_cap_frame, text="C-Terminus Cap").grid(column=1, columnspan=1,
                                                                  row=1, rowspan=1,
                                                                  padx=10, pady=10,
                                                                  sticky="w")

    c_terminus_cap = h.tk.Entry(subunit_cap_frame, width=20)
    c_terminus_cap.grid(column=1, columnspan=1, row=1, rowspan=1, padx=10, pady=10, sticky='e')

    subunit_cyclization_frame = h.tk.LabelFrame(peptide_library_main_frame)

    subunit_cyclization_frame.grid(column=0, columnspan=2, row=1, rowspan=1, sticky='new')

    subunit_cyclization_frame.grid_columnconfigure(0, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(1, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(2, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(3, weight=1)
    subunit_cyclization_frame.grid_columnconfigure(4, weight=1)

    h.tk.Label(subunit_cyclization_frame, text="Cyclization Type").grid(column=0, columnspan=1,
                                                                        row=0, rowspan=1,
                                                                        padx=10, pady=10,
                                                                        sticky='w')

    cyclization_types = ["Linear", "Head to Tail", "Depsipeptide", "Thioether", "Triazole"]

    radio_var = h.tk.StringVar()
    radio_var.set(value=cyclization_types[0])
    rs = []

    col = 0
    row = 1

    def set_cap_details():
        cyc_type = radio_var.get()
        if cyc_type == cyclization_types[0]:
            amine_cap.configure(state='normal')
            c_terminus_cap.configure(state='normal')

        if cyc_type == cyclization_types[1]:
            amine_cap.configure(state='disabled')
            c_terminus_cap.configure(state='disabled')

        if cyc_type == cyclization_types[2]:
            amine_cap.configure(state='normal')
            c_terminus_cap.configure(state='disabled')

        if cyc_type == cyclization_types[3]:
            amine_cap.configure(state='normal')
            c_terminus_cap.configure(state='disabled')

        if cyc_type == cyclization_types[4]:
            amine_cap.configure(state='disabled')
            c_terminus_cap.configure(state='normal')

    for i in range(0, len(cyclization_types)):
        r = h.tk.Radiobutton(subunit_cyclization_frame, text=cyclization_types[i],
                             value=cyclization_types[i], variable=radio_var,
                             command=set_cap_details)

        r.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        rs.append(r)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1

    cheminformatic_options_frame = h.tk.LabelFrame(peptide_library_main_frame)
    cheminformatic_options_frame.grid(column=0, columnspan=2, row=3, rowspan=1, sticky='nsew')

    h.tk.Label(cheminformatic_options_frame, text="Cheminformatics").grid(column=0, columnspan=1,
                                                                          row=0, rowspan=1,
                                                                          padx=10, pady=10,
                                                                          sticky='w')

    cheminformatic_options_frame.grid_columnconfigure(0, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(1, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(2, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(3, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(4, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(5, weight=1)

    cheminformatic_options = ["Exact Mass", "ALogP", "SLogP", "SASA", "TPSA",
                              "# Atoms", "# Heteroatoms", "# Amide Bonds",
                              "# Rings", "# Aromatic Rings", "# Rotatable Bonds", "# HBA",
                              "# HBD", "Largest Ring Size"]

    file_check_box_vars = []
    file_check_boxes = []

    col = 0
    row = 1

    for i in range(0, len(cheminformatic_options)):
        var = h.tk.StringVar()
        file_check_box_vars.append(var)

        c = h.tk.Checkbutton(cheminformatic_options_frame, text=cheminformatic_options[i],
                             variable=var, onvalue=cheminformatic_options[i], offvalue='')
        c.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        file_check_boxes.append(c)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1

    file_output_type_frame = h.tk.LabelFrame(peptide_library_main_frame)

    file_output_type_frame.grid(column=0, columnspan=2, row=7, rowspan=1, sticky='new')

    file_output_type_frame.grid_columnconfigure(0, weight=1)
    file_output_type_frame.grid_columnconfigure(1, weight=1)
    file_output_type_frame.grid_columnconfigure(2, weight=1)

    h.tk.Label(file_output_type_frame, text="Output Formats").grid(column=0, columnspan=1,
                                                                         row=0, rowspan=1,
                                                                         padx=10, pady=10,
                                                                         sticky='w')

    file_types = [".CSV", ".TXT", ".SDF (Only Mol Data)", ".SMI (Only SMILES)"]

    c_vars = []
    cs = []

    col = 0
    row = 1

    for i in range(0, len(file_types)):
        var = h.tk.StringVar()
        c_vars.append(var)

        c = h.tk.Checkbutton(file_output_type_frame, text=file_types[i],
                             variable=var, onvalue=file_types[i], offvalue='')
        c.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        cs.append(c)

        if col == 3:
            col = 0
            row += 1

        else:
            col += 1

    file_io_frame = h.tk.LabelFrame(peptide_library_main_frame)
    file_io_frame.grid(column=0, columnspan=2, row=4, rowspan=1, sticky='new')
    file_io_frame.grid_columnconfigure(0, weight=1)

    h.tk.Label(file_io_frame, text="Input / Output Files").grid(column=0, columnspan=1,
                                                                row=0, rowspan=1,
                                                                padx=10, pady=10,
                                                                sticky='w')

    def input_file_prompt():
        input_file_name = h.filedialog.askopenfile(initialdir="your directory path",
                                                   title="Input File",
                                                   filetypes=(("csv files", "*.csv"),
                                                              ("all files", "*.*")))
        input_text.set(input_file_name.name)
        h.INPUT_FILE = str(input_file_name.name)

    input_text = h.tk.StringVar()
    input_text.set(value="Select Input File...")

    input_label = h.tk.Label(file_io_frame, textvariable=input_text, text="text")

    input_label.grid(column=0, columnspan=1, row=5, rowspan=1, padx=10, pady=10, sticky='w')

    h.tk.Button(file_io_frame, text="Input File",
                command=input_file_prompt).grid(column=0, columnspan=1,
                                                row=5, rowspan=1,
                                                padx=10, pady=10,
                                                sticky="e")

    output_file_label_text = h.tk.StringVar()

    output_file_label_text.set("Select Output Folder...")

    h.tk.Label(file_io_frame, textvariable=output_file_label_text).grid(column=0, columnspan=1,
                                                                         row=6, rowspan=1,
                                                                         padx=10, pady=10,
                                                                         sticky='w')

    def output_file_prompt():
        output_dir_name = h.filedialog.askdirectory(initialdir="your directory path",
                                                    title="Output Folder")
        output_file_label_text.set(output_dir_name)
        h.OUTPUT_DIR = str(output_dir_name)

    h.tk.Button(file_io_frame, text="Output Folder",
                command=output_file_prompt).grid(column=0, columnspan=1,
                                                 row=6, rowspan=1,
                                                 padx=10, pady=10,
                                                 sticky="e")

    def create_library():
        peptide_library_progress_bar.start(10)
        cyclization_type = radio_var.get()
        amine_cap_subunit = amine_cap.get()

        if not len(amine_cap_subunit) == 0:
            h.CAP_AMINE = True

        c_terminus_cap_subunit = c_terminus_cap.get()

        if not len(amine_cap_subunit) == 0:
            h.CAP_C_TERMINUS = True

        h.OUTPUT_TYPES = []

        for j in range(0, len(file_types)):
            if not c_vars[j] == '':
                h.OUTPUT_TYPES.append(c_vars[j].get())

        h.CHEMINFORMATICS_OPTIONS = []

        for j in range(0, len(cheminformatic_options)):
            if not file_check_box_vars[j] == '':
                h.CHEMINFORMATICS_OPTIONS.append(file_check_box_vars[j].get())

        h.gui_functions.peptide_library_creation(console_list, output_list, subunit_list,
                                                 peptide_library_progress_bar, cyclization_type,
                                                 amine_cap_subunit, c_terminus_cap_subunit,
                                                 peptide_library_window)

        if not cyclization_type == "Linear":
            console_list.insert("end", f"{cyclization_type} Cyclization.")
            console_list.insert("end", "")
            console_list.see("end")

        console_list.insert("end", "Peptide Library Created Successfully!")
        console_list.insert("end", "")
        console_list.see("end")

        peptide_library_progress_bar.stop()

        peptide_library_window.destroy()

        h.INPUT_FILE = ''
        h.OUTPUT_DIR = ''

    h.tk.Button(peptide_library_main_frame, text="Create Library",
                command=create_library).grid(column=0, columnspan=1, row=0, rowspan=1, sticky="e")


def peptide_cheminformatics(root, console_list, output_list):
    console_list.insert("end", "Cheminformatics Started.")
    console_list.insert("end", "")
    console_list.see("end")

    cheminformatics_window, \
        cheminformatics_main_frame, \
        cheminformatics_progress_bar_text = create_window(root, 1080, 720, 10, 10,
                                                          "Cheminformatics")

    def cheminformatics_exit_button():
        cheminformatics_window.destroy()
        console_list.insert("end", "Cheminformatics Exited.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(cheminformatics_main_frame, text="Exit",
                 command=cheminformatics_exit_button).grid(column=1, columnspan=1,
                                                           row=0, rowspan=1,
                                                           pady=10, sticky="ne")

    smiles_string_frame = h.tk.LabelFrame(cheminformatics_main_frame)
    smiles_string_frame.grid_columnconfigure(0, weight=1)
    smiles_string_frame.grid_columnconfigure(1, weight=3)
    smiles_string_frame.grid_columnconfigure(2, weight=1)

    h.tk.Label(smiles_string_frame, text="SMILES String").grid(column=0, columnspan=1,
                                                               row=2, rowspan=1,
                                                               padx=10, pady=10,
                                                               sticky='w')
    def chem_info(events=None):
        smiles_string = smiles_string_entry.get()
        molecule = h.cheminformatics.get_molecule_from_smiles(smiles_string)
        exact_mass, tpsa, a_log_p = h.cheminformatics.get_chemometrics(molecule)
        cyclic_peptide_object = h.classes.Peptide("", exact_mass, tpsa, a_log_p, smiles_string)
        df = h.utilities.peptides_to_dataframe([cyclic_peptide_object])

        line_one = df.to_string().split("\n")[0].split()
        line_one = line_one[1] + " " + line_one[2] + "         " + \
                   line_one[3] + "          " + \
                   line_one[4] + "         " + \
                   line_one[5]

        line_two = df.to_string().split("\n")[1].split()
        line_two = line_two[1] + "        " + \
                   line_two[2] + "        " + \
                   line_two[3] + "        " + \
                   line_two[4] + "        "

        output_list.insert("end", line_one)
        output_list.insert("end", line_two)
        output_list.insert("end", "")
        output_list.see("end")
        console_list.insert("end", "Cheminformatics Evaluated Successfully.")
        console_list.insert("end", "")

    h.tk.Button(cheminformatics_main_frame, text="Evaluate",
                command=chem_info).grid(column=0, columnspan=1,
                                        row=0, rowspan=1,
                                        pady=10, padx=10,
                                        sticky="ne")

    smiles_string_frame.grid(column=0, columnspan=2, row=1, rowspan=1, sticky='new')

    smiles_string_entry = h.tk.Entry(smiles_string_frame)
    smiles_string_entry.grid(column=1, columnspan=2,
                             row=2, rowspan=1,
                             pady=10, padx=10,
                             sticky='ew')

    smiles_string_entry.bind('<Return>', chem_info)

    cheminformatic_options_frame = h.tk.LabelFrame(cheminformatics_main_frame)
    cheminformatic_options_frame.grid(column=0, columnspan=2, row=4, rowspan=1, sticky='nsew')

    cheminformatic_options_frame.grid_columnconfigure(0, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(1, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(2, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(3, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(4, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(5, weight=1)

    cheminformatic_options = ["Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "# Atoms",
                              "# Heteratoms", "# Amide Bonds", "# Rings", "# Aromatic Rings",
                              "# Rotatable Bonds", "# HBA", "# HBD"]

    c_vars = []
    cs = []

    col = 0
    row = 0

    for i in range(0, len(file_types)):
        var = h.tk.StringVar()
        c_vars.append(var)

        c = h.tk.Checkbutton(cheminformatic_options_frame, text=file_types[i],
                             variable=var, onvalue=file_types[i], offvalue='')
        c.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        cs.append(c)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1


def library_cheminformatics(root, console_list, output_list):
    console_list.insert("end", "Cheminformatics Started.")
    console_list.insert("end", "")
    console_list.see("end")

    cheminformatics_window, \
        cheminformatics_main_frame, \
        cheminformatics_progress_bar_text = create_window(root, 1080, 720, 10, 10, "Cheminformatics")

    def cheminformatics_exit_button():
        cheminformatics_window.destroy()
        console_list.insert("end", "Cheminformatics Exited.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(cheminformatics_main_frame, text="Exit",
                 command=cheminformatics_exit_button).grid(column=1, columnspan=1, row=0, rowspan=1, pady=10,
                                                          sticky="ne")

    smiles_string_frame = h.tk.LabelFrame(cheminformatics_main_frame)
    smiles_string_frame.grid_columnconfigure(0, weight=1)
    smiles_string_frame.grid_columnconfigure(1, weight=3)
    smiles_string_frame.grid_columnconfigure(2, weight=1)

    h.tk.Label(smiles_string_frame, text="SMILES String").grid(column=0, columnspan=1,
                                                                         row=2, rowspan=1,
                                                                         padx=10, pady=10,
                                                                         sticky='w')
    def chem_info(events=None):
        smiles_string = smiles_string_entry.get()
        molecule = h.cheminformatics.get_molecule_from_smiles(smiles_string)
        exact_mass, tpsa, a_log_p = h.cheminformatics.get_chemometrics(molecule)
        cyclic_peptide_object = h.classes.Peptide("", exact_mass, tpsa, a_log_p, smiles_string)
        df = h.utilities.peptides_to_dataframe([cyclic_peptide_object])

        line_one = df.to_string().split("\n")[0].split()
        line_one = line_one[1] + " " + line_one[2] + "         " + \
                   line_one[3] + "          " + \
                   line_one[4] + "         " + \
                   line_one[5]

        line_two = df.to_string().split("\n")[1].split()
        line_two = line_two[1] + "        " + \
                   line_two[2] + "        " + \
                   line_two[3] + "        " + \
                   line_two[4] + "        "

        output_list.insert("end", line_one)
        output_list.insert("end", line_two)
        output_list.insert("end", "")
        output_list.see("end")
        console_list.insert("end", "Cheminformatics Evaluated Successfully.")
        console_list.insert("end", "")

    h.tk.Button(cheminformatics_main_frame, text="Evaluate",
                command=chem_info).grid(column=0, columnspan=1,
                                        row=0, rowspan=1,
                                        pady=10, padx=10,
                                        sticky="ne")

    smiles_string_frame.grid(column=0, columnspan=2, row=1, rowspan=1, sticky='new')

    smiles_string_entry = h.tk.Entry(smiles_string_frame)
    smiles_string_entry.grid(column=1, columnspan=2, row=2, rowspan=1, pady=10, padx=10, sticky='ew')
    smiles_string_entry.bind('<Return>', chem_info)

    cheminformatic_options_frame = h.tk.LabelFrame(cheminformatics_main_frame)
    cheminformatic_options_frame.grid(column=0, columnspan=2, row=4, rowspan=1, sticky='nsew')

    cheminformatic_options_frame.grid_columnconfigure(0, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(1, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(2, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(3, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(4, weight=1)
    cheminformatic_options_frame.grid_columnconfigure(5, weight=1)

    file_types = ["Exact MW", "ALogP", "SLogP", "SASA", "TPSA", "# Atoms", "# Heteratoms", "# Amide Bonds",
                  "# Rings", "# Aromatic Rings", "# Rotatable Bonds", "# HBA", "# HBD"]

    check_var1 = h.tk.StringVar()
    check_var2 = h.tk.StringVar()
    check_var3 = h.tk.StringVar()
    check_var4 = h.tk.StringVar()
    check_var5 = h.tk.StringVar()
    check_var6 = h.tk.StringVar()
    check_var7 = h.tk.StringVar()
    check_var8 = h.tk.StringVar()
    check_var9 = h.tk.StringVar()
    check_var10 = h.tk.StringVar()
    check_var11 = h.tk.StringVar()
    check_var12 = h.tk.StringVar()
    check_var13 = h.tk.StringVar()
    check_var14 = h.tk.StringVar()
    check_var15 = h.tk.StringVar()

    c_vars = []
    cs = []

    col = 0
    row = 0

    for i in range(0, len(file_types)):
        var = h.tk.StringVar()

        c_vars.append(var)

        c = h.tk.Checkbutton(cheminformatic_options_frame, text=file_types[i],
                                variable=check_var1, onvalue=file_types[i], offvalue='')
        c.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        cs.append(c)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1


def plotting(root, console_list, output_list):
    console_list.insert("end", "3D Plotting Started.")
    console_list.insert("end", "")
    console_list.see("end")

    plotting_window, \
        plotting_main_frame, \
        plotting_progress_bar_text = create_window(root, 1080, 705, 10, 10, "3D Plotting")

    def plotting_exit_button():
        plotting_window.destroy()
        console_list.insert("end", "3D Ploting Exited.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(plotting_main_frame, text="Exit",
                 command=plotting_exit_button).grid(column=1, columnspan=1, row=0, rowspan=1, pady=10,
                                                          sticky="ne")

    x_options_frame = h.tk.LabelFrame(plotting_main_frame)
    x_options_frame.grid(column=0, columnspan=2, row=3, rowspan=1, sticky='nsew')

    h.tk.Label(x_options_frame, text="X-Axis").grid(column=0, columnspan=1,
                                                                          row=0, rowspan=1,
                                                                          padx=10, pady=10,
                                                                          sticky='w')

    x_options_frame.grid_columnconfigure(0, weight=1)
    x_options_frame.grid_columnconfigure(1, weight=1)
    x_options_frame.grid_columnconfigure(2, weight=1)
    x_options_frame.grid_columnconfigure(3, weight=1)
    x_options_frame.grid_columnconfigure(4, weight=1)
    x_options_frame.grid_columnconfigure(5, weight=1)

    cheminformatic_options = ["Exact Mass", "ALogP", "SASA", "TPSA",
                              "# Atoms", "# Heteroatoms", "# Amide Bonds",
                              "# Rings", "# Aromatic Rings", "# Rotatable Bonds", "# HBA",
                              "# HBD", "Largest Ring Size"]

    x_radio_var = h.tk.StringVar()
    x_radio_var.set(value=cheminformatic_options[0])
    x_radio_buttons = []

    col = 0
    row = 1

    def get_x_var():
        x_radio_var.get()

    for i in range(0, len(cheminformatic_options)):
        r = h.tk.Radiobutton(x_options_frame, text=cheminformatic_options[i],
                             value=cheminformatic_options[i],
                             variable=x_radio_var, command=get_x_var)

        r.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        x_radio_buttons.append(r)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1

    y_options_frame = h.tk.LabelFrame(plotting_main_frame)
    y_options_frame.grid(column=0, columnspan=2, row=4, rowspan=1, sticky='nsew')

    h.tk.Label(y_options_frame, text="Y-Axis").grid(column=0, columnspan=1,
                                                                 row=0, rowspan=1,
                                                                 padx=10, pady=10,
                                                                 sticky='w')

    y_options_frame.grid_columnconfigure(0, weight=1)
    y_options_frame.grid_columnconfigure(1, weight=1)
    y_options_frame.grid_columnconfigure(2, weight=1)
    y_options_frame.grid_columnconfigure(3, weight=1)
    y_options_frame.grid_columnconfigure(4, weight=1)
    y_options_frame.grid_columnconfigure(5, weight=1)

    y_radio_var = h.tk.StringVar()
    y_radio_var.set(value=cheminformatic_options[3])
    y_radio_buttons = []

    col = 0
    row = 1

    def get_y_var():
        y_radio_var.get()

    for i in range(0, len(cheminformatic_options)):
        r = h.tk.Radiobutton(y_options_frame, text=cheminformatic_options[i],
                             value=cheminformatic_options[i],
                             variable=y_radio_var, command=get_y_var)
        r.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        y_radio_buttons.append(r)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1

    z_options_frame = h.tk.LabelFrame(plotting_main_frame)
    z_options_frame.grid(column=0, columnspan=2, row=5, rowspan=1, sticky='nsew')

    h.tk.Label(z_options_frame, text="Z-Axis").grid(column=0, columnspan=1,
                                                                 row=0, rowspan=1,
                                                                 padx=10, pady=10,
                                                                 sticky='w')

    z_options_frame.grid_columnconfigure(0, weight=1)
    z_options_frame.grid_columnconfigure(1, weight=1)
    z_options_frame.grid_columnconfigure(2, weight=1)
    z_options_frame.grid_columnconfigure(3, weight=1)
    z_options_frame.grid_columnconfigure(4, weight=1)
    z_options_frame.grid_columnconfigure(5, weight=1)

    z_radio_var = h.tk.StringVar()
    z_radio_var.set(value=cheminformatic_options[1])
    z_radio_buttons = []

    col = 0
    row = 1

    def get_z_var():
        z_radio_var.get()

    for i in range(0, len(cheminformatic_options)):
        r = h.tk.Radiobutton(z_options_frame, text=cheminformatic_options[i],
                             value=cheminformatic_options[i],
                             variable=z_radio_var,
                             command=get_z_var)

        r.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        z_radio_buttons.append(r)

        if col == 6:
            col = 0
            row += 1

        else:
            col += 1

    file_output_type_frame = h.tk.LabelFrame(plotting_main_frame)

    file_output_type_frame.grid(column=0, columnspan=2, row=7, rowspan=1, sticky='new')

    file_output_type_frame.grid_columnconfigure(0, weight=1)
    file_output_type_frame.grid_columnconfigure(1, weight=1)
    file_output_type_frame.grid_columnconfigure(2, weight=1)
    file_output_type_frame.grid_columnconfigure(3, weight=1)
    file_output_type_frame.grid_columnconfigure(4, weight=1)

    h.tk.Label(file_output_type_frame, text="Output Formats").grid(column=0, columnspan=1,
                                                                   row=0, rowspan=1,
                                                                   padx=10, pady=10,
                                                                   sticky='w')

    file_types = [".CSV", ".TXT", ".SDF (Only Mol Data)", ".SMI (Only SMILES)"]

    c_vars = []
    cs = []

    col = 0
    row = 1

    for i in range(0, len(file_types)):
        var = h.tk.StringVar()
        c_vars.append(var)

        c = h.tk.Checkbutton(file_output_type_frame, text=file_types[i],
                             variable=var, onvalue=file_types[i], offvalue='')
        c.grid(column=col, columnspan=1, row=row, rowspan=1, padx=10, pady=10, sticky='new')

        cs.append(c)

        if col == 4:
            col = 0
            row += 1

        else:
            col += 1

    file_io_frame = h.tk.LabelFrame(plotting_main_frame)
    file_io_frame.grid(column=0, columnspan=2, row=6, rowspan=1, sticky='new')
    file_io_frame.grid_columnconfigure(0, weight=1)

    h.tk.Label(file_io_frame, text="Input / Output Files").grid(column=0, columnspan=1,
                                                                row=0, rowspan=1,
                                                                padx=10, pady=10,
                                                                sticky='w')

    def input_file_prompt():
        input_file_name = h.filedialog.askopenfile(initialdir="your directory path",
                                                   title="Input File",
                                                   filetypes=(("csv files", "*.csv"),
                                                              ("all files", "*.*")))
        input_text.set(input_file_name.name)
        h.INPUT_FILE = str(input_file_name.name)

    input_text = h.tk.StringVar()
    input_text.set(value="Select Input File...")

    input_label = h.tk.Label(file_io_frame, textvariable=input_text, text="text")

    input_label.grid(column=0, columnspan=1, row=5, rowspan=1, padx=10, pady=10, sticky='w')

    h.tk.Button(file_io_frame, text="Input File",
                command=input_file_prompt).grid(column=0, columnspan=1,
                                                row=5, rowspan=1,
                                                padx=10, pady=10,
                                                sticky="e")

    output_file_label_text = h.tk.StringVar()

    output_file_label_text.set("Select Output Folder...")

    h.tk.Label(file_io_frame, textvariable=output_file_label_text).grid(column=0, columnspan=1,
                                                                        row=6, rowspan=1,
                                                                        padx=10, pady=10,
                                                                        sticky='w')

    def output_file_prompt():
        output_dir_name = h.filedialog.askdirectory(initialdir="your directory path",
                                                    title="Output Folder")
        output_file_label_text.set(output_dir_name)
        h.OUTPUT_DIR = str(output_dir_name)

    h.tk.Button(file_io_frame, text="Output Folder",
                command=output_file_prompt).grid(column=0, columnspan=1,
                                                 row=6, rowspan=1,
                                                 padx=10, pady=10,
                                                 sticky="e")

    def plot_data():
        x = x_radio_var.get()
        y = y_radio_var.get()
        z = y_radio_var.get()

        h.gui_functions.plot_generator(x, y, z)

    h.tk.Button(plotting_main_frame, text="Plot", command=plot_data).grid(column=0, columnspan=1,
                                        row=0, rowspan=1,
                                        pady=10, padx=10,
                                        sticky="ne")


def create_main_window_main_menu(root, main_frame, console_list, output_list, subunit_list):
    main_menu_frame = h.tk.LabelFrame(main_frame)
    main_menu_frame.grid(column=0, columnspan=1, row=1, rowspan=1, sticky='nsew')
    main_menu_frame.columnconfigure(0, weight=1)

    main_menu_label =h.tk.Label(main_menu_frame, text="Main Menu")
    main_menu_label.grid(column=0, columnspan=2, row=0, rowspan=1, padx=10, pady=10, sticky='n')

    """
    def single_pep():
        single_peptide_creator(root, console_list)

    h.tk.Button(main_menu_frame, text="Peptide Creator",
                 command=single_pep, width=20).grid(column=0, columnspan=1,
                                              row=1, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame, text=": Create a Single Peptide.")
    main_menu_label.grid(column=0, columnspan=1, row=1, rowspan=1, padx=10, pady=10, sticky='e')
    """

    def pep_lib():
        peptide_library_creator(root, console_list, output_list, subunit_list)

    h.tk.Button(main_menu_frame, text="Peptide Library Creator",
                command=pep_lib, width=20).grid(column=0, columnspan=1,
                                      row=2, rowspan=1,
                                      padx=10, pady=10,
                                      sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame, text=": Create a Peptide Library by Uploading a File.")
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

    def plot_data():
        plotting(root, console_list, output_list)

    h.tk.Button(main_menu_frame, text="Plot Data",
                command=plot_data, width=20).grid(column=0, columnspan=1,
                                    row=5, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Plot Cheminformatic Data from .CSV Data File.")

    main_menu_label.grid(column=0, columnspan=1, row=5, rowspan=1, padx=10, pady=10, sticky='e')

    def drive_button():
        h.webbrowser.open('https://drive.google.com/drive/folders/1RCPQ5P0D3kjbiswlv-xwSL3Vja2ngBNP')
        console_list.insert("end", "Lokey Lab Google Drive Opened.")
        console_list.insert("end", "")
        console_list.see("end")


    h.tk.Button(main_menu_frame, text="Subunit Drive",
                 command=drive_button, width=20).grid(column=0, columnspan=1,
                                         row=6, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                 text=": Navigate to Lokey Lab Google "
                                      "Drive Containing Subunit Files.")
    main_menu_label.grid(column=0, columnspan=1, row=6, rowspan=1, padx=10, pady=10, sticky='e')

    def help_button():
        h.webbrowser.open('https://github.com/LokeyLab/Peptide-Library-Utility')
        console_list.insert("end", "Peptide Library Utility Wiki Opened.")
        console_list.insert("end", "")
        console_list.see("end")

    h.tk.Button(main_menu_frame, text="Help",
                 command=help_button, width=20).grid(column=0, columnspan=1,
                                         row=7, rowspan=1, padx=10, pady=10, sticky="w")

    main_menu_label = h.tk.Label(main_menu_frame,
                                  text=": Navigate to the Peptide Library Utility Wiki.")
    main_menu_label.grid(column=0, columnspan=1, row=7, rowspan=1, padx=10, pady=10, sticky='e')

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

    subunit_frame =h.tk.LabelFrame(main_frame)
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

    populate_subunit_library_button.grid(column=1, columnspan=1, row=0, rowspan=1, pady=10, sticky="nw")

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

    return subunit_frame, subunit_list


def create_main_window_output(main_frame, console_list):
    peptide_frame =h.tk.LabelFrame(main_frame)
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
    console_frame =h.tk.LabelFrame(main_frame)
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


def ui_loop():
    main_window, main_frame, main_progress_bar = create_root(1080, 720, 0, 0,
                                                             "Peptide Library Utility 0.2.0")

    h.tk.Button(main_frame, text="Exit",
                 command=main_window.destroy).grid(column=1, columnspan=1,
                                                   row=0, rowspan=1,
                                                   pady=10, sticky="ne")

    console_frame, console_list = create_main_window_console(main_frame)
    output_frame, output_list = create_main_window_output(main_frame, console_list)
    subunit_frame, subunit_list = create_main_window_subunits(main_frame, console_list,
                                                              main_progress_bar)
    main_menu_frame = create_main_window_main_menu(main_window, main_frame,
                                                   console_list, output_list, subunit_list)

    console_list.insert("end", "Welcome to Peptide Library Utility 0.2.0!")
    console_list.insert("end", "")
    console_list.insert("end", "Get Started by Populating the Subunit Library.")
    console_list.insert("end", "")

    main_window.mainloop()
