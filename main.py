# Internal
import bonding as b
import combinatronics as c
import subunitBuilder as s
import peptideBuilder as pb
import userInterface as ui
import utilities as u


def Main():
    cyclicPeptides.append(ui.MainUILoop())
    u.PrintPeptides(cyclicPeptides[0])
    u.PeptidesToCSV(cyclicPeptides[0])

if __name__ == '__main__':
    cyclicPeptides = []

    while True:
        Main()

# References
# Peptide Annotation
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149970/

# SMILES Guide
# https://www.epa.gov/sites/production/files/2015-05/documents/appendf.pdf

# Amino Acids
# https://www.peptide.com/resources/solid-phase-peptide-synthesis/amino-acid-derivatives-for-peptide-synthesis/amino-acid-abbreviations/
