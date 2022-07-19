# Peptide Library Utility
## Synopsis
Create linear and cyclic peptide libraries *in silico* and evaluate some basic cheminformatic values.

## Features
- Create a single peptide using a text-based user inteface. (Outputs cheminformatic data and SMILES string)
- Create a peptide library conbinatorially by inputing a .csv file and following a shortened text-based user inteface. (Writes cheminformatic data and SMILES string to .csv in the output folder.)
- Get basic cheminformatic data from a SMILES string.

## Directory Structure
```
.
├── input/ ... Directory for .CSV files containing the subunits to be made into peptides.
├── output/ ... Directory where the peptide library will be output as a .CSV file.
├── scripts/ ... Directory that contains all of the internal dependencies.
│   └── bonding.py ... Functions for bonding and cyclization.
│   └── cheminformatics.py ... Functions for calculating cheminformatic data.
│   └── classes.py ... Classes used for creating objects that hold data used in the various stages of the peptide creation process.
│   └── combinatronics.py ... Functions for calculating permutations and combinations.
│   └── header.py ... Holds all of the library imports.
│   └── smiles_strings.py ... Holds the SMILE strings for the amino acids from Pub Chem and the reorganized versions. Not used in the program.
│   └── subunit_builder.py ... Functions for generating the subunit library.
│   └── user_interface.py ... Where most of functions are called as the user navigates the text-based interface.
│   └── utilities.py ... Various utility functions.
└── subunits/ ... Directory where the .CSVs that hold the data for all of the subunits.
│   └── amines.csv ... Amines .CSV. These will automatically be bonded with a bromoacetic acid.
│   └── canonical_amino_acids.csv ... Canonical Amino Acids .CSV
│   └── noncanonical_amino_acids.csv ... Noncanonical Amino Acids .CSV
│   └── miscellaneous.csv ... Miscellaneous subunits .CSV
└── .gitattribuites
└── .gitignore
└── __main__.py ... The main script.
└── README.md
└── license ... MIT License.
```

## Installation
Below are the required libraries.
```
pip install itertools
pip install matplotlib
pip install pandas
pip install rdkit
```

## Usage
Run the __main__.py script and follow the text-based user interface to create peptides, peptide libraries, and calculate cheminformatic information from SMILES strings.

For canonical and noncanonical amino acids, the format must comply to the following examples:

```
(NMe)-ß-H-D-Leu    ß-Leu    H-Leu    Leu    (NMe)-Leu
```

For peptoids and miscellaneous subunits, use the multiple letter acronym as defined in the appropriate .CSV file.

For the input .CSV file, the format must comply to the following example:

```
Pot 1       Pot 2       Pot 3       Pot 4       Pot 5       Pot 6       <──── Lets the program know how 
                                                                                many pots there are.
Leu         Leu         Leu         Leu         Leu         Leu         <──⌍─ Can be any subunit as long as
                                                                           |    it is on the subunit list.
Pro         D-Pro       Ala         D-Ala       (NMe)-ß-Leu Gly         <──⌏ 

Val                     D-Gly                                           <──── Uneven pot sizes are acceptable.
```

## Contributors
Code written by [Adam Murray](https://github.com/Adiaslow)

In collaboration with Lokey Lab.

## License
[MIT License](https://github.com/LokeyLab/Peptide-Library-Utility/blob/main/license)

## Funding
[National Institutes of Health](https://www.nih.gov/)

[Bridges to the Baccalaureate Program](https://access.ucsc.edu/)

Grant #R25GM51765

## References
### Formating of amino acids and peptides:
- [Annotation of Peptide Structures Using SMILES and Other Chemical Codes–Practical Solutions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149970/)

### SMILES string reference:
- [Simplified molecular-input line-entry system](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
- [3. SMILES - A Simplified Chemical Language](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
- [Appendix F. SMILES Notation Tutorial](https://www.epa.gov/sites/default/files/2015-05/documents/appendf.pdf)

## Change Log
All notable changes to this project will be documented in this file. This project adheres to [Semantic Versioning](https://semver.org/).

## 0.1.0 - 2022-07-18
### Added
Initial Commit

### Changed

### Fixed

### Removed
