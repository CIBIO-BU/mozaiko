# IUPAC Dictionary: Based on Johnson, A. D. (2010).
# An extended IUPAC nomenclature code for polymorphic nucleic acids.
# Bioinformatics, 26(10), 1386-1389.
# https://doi.org/10.1093/bioinformatics/btq098

IUPAC = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def calculate_iupac_mismatches(seq1, seq2):
    """
    Calculate the number of mismatches between two sequences, according to IUPAC ambiguity codes.
    """
    mismatches = 0
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 != base2:
            if not set(IUPAC[base1]).intersection(set(IUPAC[base2])):
                mismatches += 1
    return mismatches