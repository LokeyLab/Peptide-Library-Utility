<img align="center" src="https://cpb-us-e1.wpmucdn.com/sites.ucsc.edu/dist/6/1264/files/2022/07/cropped-LOKEYLABLOGO-1.png">
  
# Peptide Library Utility
## Synopsis
Create linear and cyclic peptide libraries *in silico* and evaluate some basic cheminformatic values.

# Table of Contents
- [Features](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#features)
- [Structure](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#structure)
- [Installation](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#installation)
- [Usage](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#usage)
- [Contributors](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#contributors)
- [License](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#license)
- [Funding](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#funding)
- [References](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#references)
- [Change Log](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/README.md#change-log)

## Features
- Create a single peptide using a text-based user inteface. (Outputs cheminformatic data and SMILES string)
- Create a peptide library conbinatorially by inputing a .csv file and following a shortened text-based user inteface. (Writes cheminformatic data and SMILES string to .csv in the output folder.)
- Get basic cheminformatic data from a SMILES string.

## Structure
```
.
├── .idea ... IDE Specific information. 
├── scripts/ ... Directory that contains all of the internal dependencies.
│   └── bonding.py ... Functions for bonding and cyclization.
│   └── cheminformatics.py ... Functions for calculating cheminformatic data.
│   └── classes.py ... Classes used for creating objects that hold data used in the various stages of the peptide creation process.
│   └── combinatronics.py ... Functions for calculating permutations and combinations.
│   └── full_logo.gif
│   └── full_logo.png
│   └── gui.py ... Creates and manages the graphical user interface.
│   └── gui_functions.py ... Functions used by the gui.py script related to library generation and a few other things.
│   └── header.py ... Holds all of the library imports.
│   └── logo.ico
│   └── logo.icns
│   └── logo.png
│   └── peptide_library_utility.py ... The main script from which the gui.py script is called.
│   └── setup.py ... Checks to see if modules are installed and creates desktop shortcut.
│   └── smiles_strings.py ... Holds the SMILE strings for the amino acids from Pub Chem and the reorganized versions. Not used in the program.
│   └── subunit_builder.py ... Functions for generating the subunit library.
│   └── tbui.py ... Deprecated text-based user interface.
│   └── utilities.py ... Various utility functions.
└── subunits/ ... Directory where the .CSVs that hold the data for all of the subunits.
│   └── amines.csv ... Amines .CSV. These will automatically be bonded with a bromoacetic acid.
│   └── canonical_amino_acids.csv ... Canonical Amino Acids .CSV
│   └── noncanonical_amino_acids.csv ... Noncanonical Amino Acids .CSV
│   └── miscellaneous.csv ... Miscellaneous subunits .CSV
└── .gitattribuites
└── .gitignore
└── license ... MIT License.
└── README.md
```

## Installation
0. Have Python 3.10 installed.
1. [Download .ZIP](https://github.com/LokeyLab/Peptide-Library-Utility/archive/refs/heads/main.zip)
2. Run ../scripts/setup.py.
  a. If any python modules are not installed automatically and can not be found by the script, you will need to install them manually.  This can be done using the format:
```
pip install module_name
```
Once you are able to run setup.py with all of the modules reporting as installed, you may proceed to use the program as stated below.


## Usage
Having run setup.py properly, a desktop shortcut will be created with which you may run the program.

For canonical and noncanonical amino acids, the format must cohere with the following examples:
```
(NMe)-ß-D-Leu    ß-Leu-F    H-Leu    Leu    (NMe)-D-Leu-d3
```
- *Note that there is no ß-H-, only either ß- or H-.*

- *Standard practice for amino acids is to include any mass encoding after the amino acid three-letter.*

- For peptoids and miscellaneous subunits, use the multiple letter acronym as defined in the appropriate .CSV file.

### For the input .CSV file, the format must cohere with the following example:

```
Position 1  Position 2  Position 3  Position 4  Position 5  Position 6 <──── Lets the program know how 
                                                                                many pots there are.
Leu         Leu         Leu         Leu         Leu         Leu        <──⌍─ Can be any subunit as long as
                                                                          |     it is on the subunit list.
Pro         D-Pro       Ala         D-Ala       (NMe)-ß-Leu Gly        <──⌏ 

Val                     D-Gly                                          <──── Uneven pot sizes are acceptable.
```
- *Note that subunits should be in reverse synthetic order or N to C.*


### Anatomy of properly formatted Amino Acid SMILES string:
```
N[C@@H](CC(C)C)C(=O)O

Amino Group   Stereocenter   Side Chain   Alpha Carbon    Carboxyl Group
N             [C@@H]         (CC(C)C)     C               (=O)O
                             ^      ^
                             (Some side chains may or may not require surrounding parentheses)
```

### Adding deuterated amino acids to noncanonical-amino-acids.csv:
```
Improperly formatted L-Leucine:
CC(C)C[C@H](N)C(O)=O

Improperly formatted L-Leucine-d3:
[2H]C([2H])([2H])C(C)C[C@H](N)C(O)=O

Properly formatted L-Leucine: 
N[C@@H](CC(C)C)C(=O)O

Properly formatted L-Leucine-d3: 
N[C@@H](CC(C)C([2H])([2H])([2H]))C(=O)O
```

## Contributors
- Code written by [Adam Murray](https://github.com/Adiaslow)
- Coding help from [Akshar Lohith](https://github.com/alohith)

In collaboration with Lokey Lab.

## License
[MIT License](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/license)

## Funding
[National Institutes of Health](https://www.nih.gov/)

[Bridges to the Baccalaureate Program](https://access.ucsc.edu/)

- Grant #R25GM51765

## References
### Formating of amino acids and peptides:
- [Annotation of Peptide Structures Using SMILES and Other Chemical Codes–Practical Solutions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149970/)

### SMILES string reference:
- [Simplified molecular-input line-entry system](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
- [3. SMILES - A Simplified Chemical Language](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
- [Appendix F. SMILES Notation Tutorial](https://www.epa.gov/sites/default/files/2015-05/documents/appendf.pdf)

## Change Log
All notable changes to this project will be documented in this file. This project adheres to [Semantic Versioning](https://semver.org/).

## 1.0.1 - 2022-08-5
### Added
- Ala-d4 and variants

### Changed

### Fixed
- Val and Bip

### Removed
- input and output folder

## 1.0.0 - 2022-08-1
### Added
- GUI added
- setup.py script for the initialization of the program for first time users.
- Desktop shortcut to run program without using an IDE.

### Changed
- __main__.py renamed to peptide_library_utility.py
- Input files are now taken in reverse syntheic order.

### Fixed
- Minor code clean up
- subunit_builder .CSV file importing is now fixed.

### Removed
- input and output folder

## 0.2.0 - 2022-07-23
### Added
- Beta and Homo versions of amino acids (still working on this)

### Changed
- Made the alpha carbon editable (useful for deuterated amino acids)
- Moved __main__.py to scripts folder

### Fixed
- Minor code clean up

### Removed
- Classes from subunits.py

## 0.1.0 - 2022-07-18
### Added
- Initial Commit

### Changed

### Fixed

### Removed
