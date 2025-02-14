import glob
import os
import re
import subprocess
from pathlib import Path
from typing import Union
from collections import defaultdict

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
        print(f"mozaiko WARNING: Empty data in amplicon or insert file for {amplicon_basename}")
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
                    print(f"mozaiko WARNING: Empty binding site extracted for ID {seq_id}")
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


def sequence_count_tracking(
    original_database, analysis_folder, save_results: bool = True
):
    """
    This method tracks sequence count per each analysis step.

    Parameters:
    - original_database: path to the FASTA file containing the original inputted database.
    - analysis_folder: path to the folder contaning the analysis outcomes.
    - save_results: bool
        A boolean parameter to save the results to an CSV file, if the parameter is set to True.

    Output:
    - sequence_count_track (Dataframe): TSV file containing the number of sequenc21h45 	22h03es in the original
    database and the number of sequences considered after each analysis step.
    """
    try:
        file_list = []

        for root, dirs, files in os.walk(analysis_folder):
            for file in files:
                file_list.append(os.path.join(root, file))

        if os.path.exists(original_database):
            file_list.append(original_database)
            original_dir = "original_database"
        else:
            print(f"mozaiko WARNING: Original database file not found.")
            original_dir = "original_database_not_found"

        sequence_counts = {}

        for file_path in file_list:
            try:
                count_result = subprocess.run(
                    ["grep", "-c", "^>", file_path],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                count = int(count_result.stdout.strip())
                sequence_counts[file_path] = count
            except subprocess.CalledProcessError as e:
                print(f"mozaiko ERROR: Error processing file {file_path}: {e}")
                sequence_counts[file_path] = np.na
            except ValueError as e:
                print(
                    f"mozaiko ERROR: Error counting sequences for file {file_path}: {e}"
                )
                sequence_counts[file_path] = np.na

        records = []

        for path, count in sequence_counts.items():
            if path == original_database:
                primer_name = "original_database"
                directory = original_dir

            else:
                primer_name = ((os.path.basename(path)).split("."))[0]
                directory = path.split("/")[-2]
                if directory == "filtered":
                    directory = (path.split("/")[-3]) + "-" + directory

            sequence_count = count

            records.append(
                {
                    "primer_name": primer_name,
                    "analysis_step": directory,
                    "number_of_sequences": sequence_count,
                }
            )

        df = pd.DataFrame(records)
        df_pivoted = df.pivot(
            index="primer_name", columns="analysis_step", values="number_of_sequences"
        )

        analysis_order = [
            "original_database",
            "amplicon",
            "amplicon-filtered",
            "insert",
            "insert-filtered",
            "all_complete_pbs",
            "all_complete_pbs-filtered",
            "all_inserts",
            "all_inserts_filtered",
        ]
        existing_cols = [col for col in analysis_order if col in df_pivoted.columns]

        df_pivoted = df_pivoted[existing_cols].fillna("NA")

        if save_results == True:
            run_name = os.path.basename(analysis_folder)
            output_name = "sequence_count_track.tsv"
            output_path = analysis_folder + "/" + run_name + "-" + output_name

            df_pivoted.to_csv(
                output_path, sep="\t", index=True, header=True, index_label=""
            )

            print(
                f"mozaiko INFO: Sequence count tracking file successfully saved to {output_path}."
            )

        return df_pivoted

    except Exception as e:
        print(f"mozaiko ERROR: Unexpected error occurred: {e}")
        return None


# def create_MultiBarcodeTools_input(insert_folder, pbs_incomplete_folder, output_file):
#     """
#     Process all FASTA files in a given folder and write extracted information to a TSV file.

#     Parameters:
#     - folder_path:Path to the folder containing FASTA files
#     - output_file: Path to the output TSV file
#     """
#     try:
#         fasta_files = glob.glob(os.path.join(insert_folder, "*.fasta")) + glob.glob(
#             os.path.join(pbs_incomplete_folder, "*.fasta")
#         )

#         if not fasta_files:
#             raise FileNotFoundError(
#                 "mozaiko ERROR: No FASTA files found in the specified folders."
#             )

#         with open(output_file, "w") as tsv_file:

#             for fasta_path in fasta_files:
#                 primer_name = os.path.splitext(os.path.basename(fasta_path))[0]

#                 with open(fasta_path, "r") as fasta:

#                     current_header = None
#                     current_sequence = []

#                     for line in fasta:
#                         line = line.strip()

#                         if line.startswith(">"):
#                             if current_header and current_sequence:
#                                 process_sequence(
#                                     current_header,
#                                     current_sequence,
#                                     primer_name,
#                                     tsv_file,
#                                 )

#                             current_header = line[1:]
#                             current_sequence = []

#                         elif line:
#                             current_sequence.append(line)

#                     if current_header and current_sequence:
#                         process_sequence(
#                             current_header, current_sequence, primer_name, tsv_file
#                         )

#         print(f"mozaiko INFO: MultiBarcodeTools file created to {output_file}")
#         return output_file

#     except Exception as e:
#         print(f"mozaiko ERROR: Error creating MultiBarcodeTools input - {str(e)}")
#         return None

def split_fasta_by_family(fasta_path, output_dir):
    """
    Split a FASTA file into multiple files based on family information in the header.

    Parameters:
    - fasta_path: Path to input FASTA file
    - output_dir: Directory to store family-specific FASTA files

    Returns:
    - dict: Mapping of family names to their output file paths
    """
    family_sequences = defaultdict(list)
    current_header = None
    current_sequence = []

    # Read and sort sequences by family
    with open(fasta_path, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_sequence:
                    family = current_header.split('|')[-3].strip()
                    family_sequences[family].extend([current_header, ''.join(current_sequence)])
                current_header = line
                current_sequence = []
            elif line:
                current_sequence.append(line)

        # Don't forget the last sequence
        if current_header and current_sequence:
            family = current_header.split('|')[-3].strip()
            family_sequences[family].extend([current_header, ''.join(current_sequence)])

    # Write family-specific FASTA files
    family_files = {}
    primer_name = os.path.splitext(os.path.basename(fasta_path))[0]

    for family, sequences in family_sequences.items():
        family_safe_name = family.replace(' ', '_')
        output_path = os.path.join(output_dir, f"{family_safe_name}_{primer_name}.fasta")

        with open(output_path, 'w') as f:
            for i in range(0, len(sequences), 2):
                f.write(f"{sequences[i]}\n{sequences[i+1]}\n")

        family_files[family] = output_path

    return family_files

def create_family_multibarcode_input(family_fasta_files, output_dir):
    """
    Create MultiBarcodeTools input files for each family.

    Parameters:
    - family_fasta_files: Dict mapping families to their FASTA file paths
    - output_dir: Directory to store the output files

    Returns:
    - dict: Mapping of family names to their MultiBarcodeTools input files
    """
    family_inputs = {}

    for family, fasta_files in family_fasta_files.items():
        family_safe_name = family.replace(' ', '_')
        output_file = os.path.join(output_dir, f"multibarcode_input_{family_safe_name}.tsv")

        with open(output_file, 'w') as tsv_file:
            for fasta_path in fasta_files:
                primer_name = os.path.splitext(os.path.basename(fasta_path))[0]
                with open(fasta_path, 'r') as fasta:
                    current_header = None
                    current_sequence = []
                    for line in fasta:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_header and current_sequence:
                                process_sequence(
                                    current_header,
                                    current_sequence,
                                    primer_name,
                                    tsv_file
                                )
                            current_header = line[1:]
                            current_sequence = []
                        elif line:
                            current_sequence.append(line)
                    if current_header and current_sequence:
                        process_sequence(
                            current_header,
                            current_sequence,
                            primer_name,
                            tsv_file
                        )

        family_inputs[family] = output_file

    return family_inputs

def process_sequence(header, sequence_lines, primer_name, tsv_file):
    """
    Process a single FASTA sequence and write to TSV.

    Parameters:
    - header : FASTA header line
    - sequence_lines : ist of sequence lines
    - barcode_name : Name of the barcode/primer
    - tsv_out : Output TSV file handle
    """
    full_sequence = "".join(sequence_lines)

    if "|" in header:
        parts = header.split("|")

        if len(parts) >= 10:
            seq_ID = parts[0].strip()
            species_name = parts[2].strip()

            tsv_file.write(
                f"{seq_ID}\t{primer_name}\t{species_name}\t{full_sequence}\n"
            )
        else:
            print(f"mozaico WARNING: Incorrect header format: {header}. Check if the header has" +
                  "been through harmonization.")
    else:
        print(f"mozaico WARNING: Unexpected header format: {header}")
