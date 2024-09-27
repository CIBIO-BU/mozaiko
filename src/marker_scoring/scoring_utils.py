from collections import defaultdict

from Bio import SeqIO

# IUPAC Dictionary: Based on Johnson, A. D. (2010).
# An extended IUPAC nomenclature code for polymorphic nucleic acids.
# Bioinformatics, 26(10), 1386-1389.
# https://doi.org/10.1093/bioinformatics/btq098

IUPAC = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


def calculate_iupac_mismatches(sequence1, sequence2):
    """
    Calculate the number of mismatches between two sequences, according to IUPAC ambiguity codes.

    Parameters:
    sequence1 (str): First sequence to be analysed.
    sequence2 (str): Second sequence to be analysed against.

    Returns:
    int: Number of mismatches found between two sequences.
    """
    mismatches = 0

    for base1, base2 in zip(sequence1.upper(), sequence2.upper()):
        if base1 != base2:
            # check if the bases are compatible according to IUPAC
            # by taking the union of the sets of compatible bases
            # and checking if the intersection is empty
            if not set(IUPAC.get(base1, [base1])).intersection(
                set(IUPAC.get(base2, [base2]))
            ):
                mismatches += 1

    return mismatches


def detect_fwd_rev_primer_len(sequence):
    """
    Calculate the length of the forward and reverse primer.

    Parameters:
    sequence (str): DNA sequence.

    Returns:
    int: Length of the forward and reverse primer.
    """

    fwd_len = len(sequence) // 2
    rev_len = len(sequence) - fwd_len

    return fwd_len, rev_len


def calculate_ambiguous_percentage(sequence):
    """
    Calculate the percentage of ambiguous bases in a DNA sequence.

    Parameters:
    sequence (str): The DNA sequence.

    Returns:
    float: The percentage of ambiguous bases in the sequence.
    """
    ambiguous_bases = set("RYWSMKHBVDN")

    return sum(base in ambiguous_bases for base in sequence) / len(sequence)


def write_filtered_sequence(output_handle, record):
    """
    Write a filtered sequence to the output file.

    Parameters:
    - output_handle:
    -record:

    Returns:
    - None
    """
    sequence = str(record.seq)
    output_handle.write(f">{record.description}\n{sequence}\n")


def filter_sequences_by_ambiguity(input_file, out_file, max_ambiguous_percentage=0.05):
    """
    Filter DNA sequences based on the maximum allowed percentage of ambiguous bases.

    Parameters:
    - input_file (str): Path to the input file containing DNA sequences in FASTA format.
    - out_file (str): Path to the output file to write the filtered sequences.
    - max_ambiguous_percentage (float): Maximum allowed percentage of ambiguous bases.

    Returns:
    - None
    """

    with open(out_file, "w", encoding="UTF-8") as output_handle:

        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)
            ambiguous_percentage = calculate_ambiguous_percentage(sequence)

            if ambiguous_percentage <= max_ambiguous_percentage:
                write_filtered_sequence(output_handle, record)


def extract_primers(input_file):
    primers = defaultdict(lambda: ["", "", 0])

    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)

        fwd_len, rev_len = detect_fwd_rev_primer_len(sequence)
        accession_number = record.id
        # stores forward primer
        primers[accession_number][0] = sequence[:fwd_len]
        # stores reverse primer
        primers[accession_number][1] = sequence[-rev_len:]
        # stores sequence lenght (total lenght minus lenght of both primers)
        primers[accession_number][2] = len(sequence) - fwd_len - rev_len
    return primers
