import glob
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# IUPAC Dictionary: Based on Johnson, A. D. (2010).
# An extended IUPAC nomenclature code for polymorphic nucleic acids.
# Bioinformatics, 26(10), 1386-1389.
# https://doi.org/10.1093/bioinformatics/btq098


def calculate_iupac_mismatches(sequence1, sequence2, search_gc_clamp: bool = False):
    """
    Calculate the number of mismatches between two sequences, according to IUPAC ambiguity codes.

    Parameters:
    sequence1 (str): First sequence to be analysed.
    sequence2 (str): Second sequence to be analysed against.

    Returns:
    int: Number of mismatches found between two sequences.
    """
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

    mismatches = 0
    gc_matches = 0

    for base1, base2 in zip(sequence1.upper(), sequence2.upper()):
        # check if the bases are compatible according to IUPAC
        # by taking the union of the sets of compatible bases
        # and checking if the intersection is empty
        if not set(IUPAC.get(base1, {base1})).intersection(
            set(IUPAC.get(base2, {base2}))
        ):
            mismatches += 1

        if search_gc_clamp and (
            (base1 == "G" and base2 == "C") or (base1 == "C" and base2 == "G")
        ):
            gc_matches += 1

    return (mismatches, gc_matches) if search_gc_clamp else mismatches


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
    - output_handle: File handle for writing
    - record: SeqIO record object
    """
    output_handle.write(f">{record.description}\n{record.seq}\n")


def filter_sequences_by_ambiguity(
    input_path, output_dir=None, max_ambiguous_percentage=0.05
):
    """
    Filter DNA sequences based on the maximum allowed percentage of ambiguous bases.
    Can process either a single file or all FASTA files in a directory.

    Parameters:
    - input_path (str or Path): Path to input file or directory containing FASTA files
    - output_dir (str or Path): Directory to write the filtered files (optional)
    - max_ambiguous_percentage (float): Maximum allowed percentage of ambiguous bases (0.0 to 1.0)

    Returns:
    - dict: Mapping of input files to their corresponding output files
    """
    input_path = Path(input_path)

    if input_path.is_file():
        input_files = [input_path]
    elif input_path.is_dir():
        input_files = list(input_path.glob("*.fasta"))
    else:
        raise ValueError(f"mozaiko ERROR: Input path {input_path} does not exist")

    if not input_files:
        raise ValueError(f"mozaiko ERROR: No FASTA files found in {input_path}")

    if not (0.0 <= max_ambiguous_percentage <= 1.0):
        raise ValueError("mozaiko ERROR: max_ambiguous_percentage must be between 0.0 and 1.0")

    if output_dir is None:
        if input_path.is_dir():
            output_dir = input_path / "filtered"
        else:
            output_dir = input_path.parent / "filtered"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    consecutive_Ns = re.compile(r"N{4,}")

    # processed_files = {}
    for input_file in input_files:
        output_filename = f"{input_file.name}"
        output_path = output_dir / output_filename

        with open(output_path, "w", encoding="UTF-8") as output_handle:
            for record in SeqIO.parse(input_file, "fasta"):
                upper_sequence = str(record.seq.upper())
                ambiguous_percentage = calculate_ambiguous_percentage(upper_sequence)

                has_consecutive_Ns = consecutive_Ns.search(upper_sequence) is not None

                if (
                    ambiguous_percentage <= max_ambiguous_percentage
                    and not has_consecutive_Ns
                ):
                    write_filtered_sequence(output_handle, record)

    # print(f"mozaiko INFO: Successfully filtered ambiguous sequences in {input_path} to {output_dir}.")


def read_fasta(file):
    """
    Helper function to read a FASTA file and return a dictionary of sequence.
    """
    sequences = {}

    with open(file, "r") as f:
        header = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # New record
                if header is not None:
                    sequences[header] = sequence
                header = line
                sequence = ""
            else:
                sequence += line
        if header is not None:
            sequences[header] = sequence
    return sequences


def extract_primer_binding_sites(amplicon_file, insert_file):
    """
    This method overlaps the insert sequence with the amplicon sequence to extract the primer
    binding sites. It does so by defining the foward adpter sequence as the sequence before the
    insert, and the reverse adapter sequence as the sequence after the insert.

    The function assumes that sequence headers in both files contain matching identifiers
    that can be used to pair the correct sequences together.
    """
    amplicon_basename = str(amplicon_file).split("/")[-1]

    amplicon_data = read_fasta(amplicon_file)
    insert_data = read_fasta(insert_file)

    results = []

    def get_seq_id(header):
        clean_header = header.replace(">", "")
        return clean_header.split("|")[0]

    amplicon_dict = {
        get_seq_id(header): (header, seq) for header, seq in amplicon_data.items()
    }
    insert_dict = {
        get_seq_id(header): (header, seq) for header, seq in insert_data.items()
    }

    if not amplicon_dict or not insert_dict:
        print(
            f"mozaiko WARNING: Empty data in amplicon or insert file for {amplicon_basename}"
        )
        return pd.DataFrame(
            columns=["header", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"]
        )

    for seq_id in insert_dict:
        if seq_id in amplicon_dict:
            amplicon_header, amplicon_sequence = amplicon_dict[seq_id]
            insert_header, insert_sequence = insert_dict[seq_id]

            if not amplicon_sequence or not insert_sequence:
                print(f"mozaiko WARNING: Empty sequence found for ID {seq_id}")
                continue

            start_index = amplicon_sequence.find(insert_sequence)

            if start_index != -1:
                fwd_seq = amplicon_sequence[:start_index]
                rev_seq = amplicon_sequence[start_index + len(insert_sequence) :]

                # append if both sequences are non-empty
                if fwd_seq and rev_seq:
                    fwd_seq_len = len(fwd_seq)
                    rev_seq_len = len(rev_seq)
                    results.append(
                        [amplicon_header, fwd_seq, rev_seq, fwd_seq_len, rev_seq_len]
                    )
                else:
                    print(
                        f"mozaiko WARNING: Empty binding site extracted for ID {seq_id}"
                    )
            else:
                print(
                    f"mozaiko WARNING: Insert sequence not found within amplicon for ID {seq_id} in {amplicon_basename}"
                )
        else:
            print(
                f"mozaiko WARNING: No matching amplicon found for insert ID {seq_id} in {amplicon_basename}"
            )

    if not results:
        print(
            "mozaiko WARNING: No matching sequences found between amplicon and insert files."
        )
        return pd.DataFrame(
            columns=["header", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"]
        )

    primer_dataframe = pd.DataFrame(
        results, columns=["header", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"]
    )
    return primer_dataframe

def remove_rc_suffix_from_fasta_files(results_directory): # make this static?
    """
    Recursively processes all FASTA files in the given directory and its subdirectories,
    removing ' rc' suffix from the end of headers.

    Parameters:
        results_directory (str): Path to the results directory to search for FASTA files

    Returns:
        dict: Statistics about the processing operation
    """
    stats = {"files_processed": 0, "headers_modified": 0, "errors": []}

    # Pattern to match FASTA headers ending with ' rc'
    # This ensures we only match ' rc' when it's at the very end of the header
    header_pattern = re.compile(r"^(>.+)\s+rc$")

    for dirpath, dirnames, filenames in os.walk(results_directory):
        for filename in filenames:
            if filename.endswith((".fasta")):
                file_path = os.path.join(dirpath, filename)
                temp_file_path = file_path + ".temp"

                try:
                    headers_changed = 0

                    with open(file_path, "r") as input_file, open(
                        temp_file_path, "w"
                    ) as output_file:
                        for line in input_file:
                            if line.startswith(">"):
                                # Check if the header ends with ' rc'
                                match = header_pattern.match(line.strip())
                                if match:
                                    # Replace the header without the ' rc' suffix
                                    new_header = match.group(1) + "\n"
                                    output_file.write(new_header)
                                    headers_changed += 1
                                else:
                                    output_file.write(line)
                            else:
                                output_file.write(line)

                    if headers_changed > 0:
                        os.replace(temp_file_path, file_path)
                        stats["headers_modified"] += headers_changed
                    else:
                        os.remove(temp_file_path)

                    stats["files_processed"] += 1

                except Exception as e:
                    stats["errors"].append(f"Error processing {file_path}: {str(e)}")
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)

    return stats
