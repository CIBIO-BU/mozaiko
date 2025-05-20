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


# def sequence_count_tracking(
#     original_database, analysis_folder, save_results: bool = True
# ):
#     """
#     This method tracks sequence count per each analysis step.
#     Parameters:
#     - original_database: path to the FASTA file containing the original inputted database.
#     - analysis_folder: path to the folder contaning the analysis outcomes.
#     - save_results: bool
#     A boolean parameter to save the results to an CSV file, if the parameter is set to True.
#     Output:
#     - sequence_count_track (Dataframe): TSV file containing the number of sequences in the original
#     database and the number of sequences considered after each analysis step.
#     """
#     try:
#         # Define the specific folders to analyze and their clean names for columns
#         target_folders = ["all_complete_pbs", "amplicon", "incomplete_pbs", "insert"]

#         file_list = []
#         # Walk through the analysis folder
#         for root, dirs, files in os.walk(analysis_folder):
#             # Check if this directory is or is a subdirectory of any target folder
#             path_parts = root.split(os.path.sep)

#             # Only process files if they're in target folder or subfolders
#             is_target_or_subfolder = any(folder in path_parts for folder in target_folders)

#             if is_target_or_subfolder:
#                 for file in files:
#                     file_list.append(os.path.join(root, file))

#         # Add the original database if it exists
#         if os.path.exists(original_database):
#             file_list.append(original_database)
#             original_dir = "original_database"
#         else:
#             print(f"mozaiko WARNING: Original database file not found.")
#             original_dir = "original_database_not_found"

#         # Dictionary to track sequence counts across folders and primers
#         primer_step_counts = {}

#         # Process each file to count sequences
#         for file_path in file_list:
#             try:
#                 # Count sequences in the file
#                 count_result = subprocess.run(
#                     ["grep", "-c", "^>", file_path],
#                     capture_output=True,
#                     text=True,
#                     check=True,
#                 )
#                 count = int(count_result.stdout.strip())

#                 # Extract primer name from filename
#                 filename = os.path.basename(file_path)
#                 primer_name = filename.split('.')[0]

#                 # Determine the analysis step (folder)
#                 if file_path == original_database:
#                     primer_name = "original_database"
#                     step = "original_database"
#                 else:
#                     # Determine which target folder this belongs to
#                     step = None
#                     path_parts = file_path.split(os.path.sep)

#                     for folder in target_folders:
#                         if folder in path_parts:
#                             step = folder
#                             # Check if it's in a filtered subfolder
#                             folder_index = path_parts.index(folder)
#                             if folder_index + 1 < len(path_parts) and path_parts[folder_index + 1] == "filtered":
#                                 step = f"{folder}-filtered"
#                             break

#                 # Skip if we couldn't determine the step
#                 if step is None:
#                     continue

#                 # Create a simplified primer name by removing any folder prefix
#                 if primer_name != "original_database":
#                     # Remove folder prefixes if present (like "insert-" or "all_complete_pbs-")
#                     for folder in target_folders:
#                         prefix = f"{folder}-"
#                         if primer_name.startswith(prefix):
#                             primer_name = primer_name[len(prefix):]

#                 # Store in our dictionary, keeping the highest count if duplicates exist
#                 key = (primer_name, step)
#                 if key not in primer_step_counts or primer_step_counts[key] < count:
#                     primer_step_counts[key] = count

#             except subprocess.CalledProcessError as e:
#                 print(f"mozaiko ERROR: Error processing file {file_path}: {e}")
#             except ValueError as e:
#                 print(f"mozaiko ERROR: Error counting sequences for file {file_path}: {e}")

#         # Convert to records for DataFrame creation
#         records = []
#         for (primer_name, step), count in primer_step_counts.items():
#             records.append({
#                 "primer_name": primer_name,
#                 "analysis_step": step,
#                 "number_of_sequences": count
#             })

#         # Create DataFrame and pivot
#         df = pd.DataFrame(records)

#         # Handle empty DataFrame case
#         if df.empty:
#             print("mozaiko WARNING: No data found to create sequence count tracking.")
#             return None

#         df_pivoted = df.pivot(
#             index="primer_name",
#             columns="analysis_step",
#             values="number_of_sequences"
#         )

#         # Define column order with original_database first
#         column_order = ["original_database"]
#         # Add target folders in specified order
#         column_order.extend(target_folders)
#         # Add filtered versions
#         column_order.extend([f"{folder}-filtered" for folder in target_folders])

#         # Keep only columns that exist in our data
#         existing_cols = [col for col in column_order if col in df_pivoted.columns]

#         # Reorder columns and fill NaN values with "NA" string
#         df_pivoted = df_pivoted[existing_cols].fillna("NA")

#         if save_results:
#             run_name = os.path.basename(analysis_folder)
#             output_name = "sequence_count_track.tsv"
#             output_path = os.path.join(analysis_folder, f"{run_name}-{output_name}")

#             df_pivoted.to_csv(
#                 output_path, sep="\t", index=True, header=True, index_label=""
#             )
#             print(
#                 f"mozaiko INFO: Sequence count tracking file successfully saved to {output_path}."
#             )

#         return df_pivoted

#     except Exception as e:
#         print(f"mozaiko ERROR: Unexpected error occurred: {e}")
#         import traceback
#         traceback.print_exc()
#         return None

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
    - sequence_count_track (Dataframe): TSV file containing the number of sequences in the original
    database and the number of sequences considered after each analysis step.
    """
    try:
        # Define the specific target folders to analyze
        target_folders = ["all_complete_pbs", "amplicon", "incomplete_pbs", "insert", "input_B"]

        # Dictionary to store file paths and their analysis steps
        step_files = {}

        # Walk through the analysis folder to find all relevant files
        for root, dirs, files in os.walk(analysis_folder):
            # Get path components
            path_parts = root.split(os.path.sep)

            # Check if any target folder is in the path
            target_folder_present = False
            target_folder_path = None

            for target in target_folders:
                if target in path_parts:
                    target_folder_present = True
                    # Get the index of the target folder in the path
                    target_index = path_parts.index(target)
                    # Extract the target and all subsequent folders to create the step name
                    target_folder_path = path_parts[target_index:]
                    break

            # Skip if no target folder found
            if not target_folder_present:
                continue

            # Process files in this directory
            for file in files:
                # Only process files with extensions that could be sequence files
                if file.endswith((".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn")):
                    file_path = os.path.join(root, file)

                    # Create a standardized analysis step name based on the subfolder path
                    # Join all path components after the target folder with hyphens
                    step_name = "-".join(target_folder_path)

                    # Add to our collection
                    if step_name not in step_files:
                        step_files[step_name] = []
                    step_files[step_name].append(file_path)

        # Add the original database to our analysis steps
        if os.path.exists(original_database):
            if "original_database" not in step_files:
                step_files["original_database"] = []
            step_files["original_database"].append(original_database)
        else:
            print(f"mozaiko WARNING: Original database file not found.")

        # Dictionary to store sequence counts by primer and step
        primer_step_counts = {}

        # Process each analysis step and count sequences in each file
        for step_name, file_paths in step_files.items():
            for file_path in file_paths:
                try:
                    # Count sequences in the file
                    count_result = subprocess.run(
                        ["grep", "-c", "^>", file_path],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    count = int(count_result.stdout.strip())

                    # Extract primer name from filename
                    filename = os.path.basename(file_path)
                    file_base_name = filename.split('.')[0]

                    # Handle the original database case
                    if step_name == "original_database":
                        primer_name = "original_database"
                    else:
                        # Extract primer name, removing any folder prefixes
                        primer_name = file_base_name
                        for prefix in [f"{folder}-" for folder in target_folders]:
                            if primer_name.startswith(prefix):
                                primer_name = primer_name[len(prefix):]

                    # Use the primer name and step as a key
                    key = (primer_name, step_name)

                    # Store the highest count if there are duplicates
                    if key not in primer_step_counts or primer_step_counts[key] < count:
                        primer_step_counts[key] = count

                except subprocess.CalledProcessError as e:
                    print(f"mozaiko ERROR: Error processing file {file_path}: {e}")
                except ValueError as e:
                    print(f"mozaiko ERROR: Error counting sequences for file {file_path}: {e}")

        # Convert to records for DataFrame creation
        records = []
        for (primer_name, step), count in primer_step_counts.items():
            records.append({
                "primer_name": primer_name,
                "analysis_step": step,
                "number_of_sequences": count
            })

        # Create DataFrame and pivot
        df = pd.DataFrame(records)

        # Handle empty DataFrame case
        if df.empty:
            print("mozaiko WARNING: No data found to create sequence count tracking.")
            return None

        df_pivoted = df.pivot(
            index="primer_name",
            columns="analysis_step",
            values="number_of_sequences"
        )

        # Define a clean column order with original_database first
        preferred_columns = ["original_database"]

        # Add base target folders
        preferred_columns.extend(target_folders)

        # Add known filtered versions and other common patterns
        for folder in target_folders:
            preferred_columns.append(f"{folder}-filtered")

        # Get the remaining columns that weren't specified in our preferred order
        remaining_columns = [col for col in df_pivoted.columns if col not in preferred_columns]

        # Combine preferred columns (that exist) with remaining columns (sorted)
        final_columns = [col for col in preferred_columns if col in df_pivoted.columns]
        final_columns.extend(sorted(remaining_columns))

        # Reorder columns and fill NaN values with "NA" string
        df_pivoted = df_pivoted[final_columns].fillna("NA")

        if save_results:
            run_name = os.path.basename(analysis_folder)
            output_name = "sequence_count_track.tsv"
            output_path = os.path.join(analysis_folder, f"{run_name}-{output_name}")

            df_pivoted.to_csv(
                output_path, sep="\t", index=True, header=True, index_label=""
            )
            print(
                f"mozaiko INFO: Sequence count tracking file successfully saved to {output_path}."
            )

        return df_pivoted

    except Exception as e:
        print(f"mozaiko ERROR: Unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        return None


def remove_rc_suffix_from_fasta_files(results_directory):
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


def create_MultiBarcodeTools_input(insert_folder, pbs_incomplete_folder, output_file):
    """
    Process all FASTA files in a given folder and write extracted information to a TSV file.

    Parameters:
    - folder_path:Path to the folder containing FASTA files
    - output_file: Path to the output TSV file
    """
    try:
        fasta_files = glob.glob(os.path.join(insert_folder, "*.fasta")) + glob.glob(
            os.path.join(pbs_incomplete_folder, "*.fasta")
        )

        if not fasta_files:
            raise FileNotFoundError(
                "mozaiko ERROR: No FASTA files found in the specified folders."
            )

        with open(output_file, "w") as tsv_file:

            for fasta_path in fasta_files:
                primer_name = os.path.splitext(os.path.basename(fasta_path))[0]

                with open(fasta_path, "r") as fasta:

                    current_header = None
                    current_sequence = []

                    for line in fasta:
                        line = line.strip()

                        if line.startswith(">"):
                            if current_header and current_sequence:
                                process_sequence(
                                    current_header,
                                    current_sequence,
                                    primer_name,
                                    tsv_file,
                                )

                            current_header = line[1:]
                            current_sequence = []

                        elif line:
                            current_sequence.append(line)

                    if current_header and current_sequence:
                        process_sequence(
                            current_header, current_sequence, primer_name, tsv_file
                        )

        print(f"mozaiko INFO: MultiBarcodeTools file created to {output_file}")
        return output_file

    except Exception as e:
        print(f"mozaiko ERROR: Error creating MultiBarcodeTools input - {str(e)}")
        return None


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

        if len(parts) == 11:
            seq_ID = parts[0].strip()
            species_name = parts[2].strip()

            tsv_file.write(
                f"{seq_ID}\t{primer_name}\t{species_name}\t{full_sequence}\n"
            )
        else:
            print(
                f"mozaico WARNING: Incorrect header format: {header}. Check if the header has"
                + "been through harmonization."
            )
    else:
        print(f"mozaico WARNING: Unexpected header format: {header}")


def process_bmi_sequences(
    header, sequence_lines, primer_name, tsv_file, skipped_counter
):
    full_sequence = "".join(sequence_lines)
    if "|" in header:
        parts = header.split("|")
        original_length = len(parts)
        if original_length == 11:
            seq_ID = parts[0].strip()
            species_name = parts[2].strip()
            kingdom = parts[4].strip()
            phylum = parts[5].strip()
            classt = parts[6].strip()
            order = parts[7].strip()
            family = parts[8].strip()
            genus = parts[9].strip()
            species = parts[10].strip()
            tsv_file.write(
                f"{seq_ID}\t{primer_name}\t{species_name}\t{kingdom}\t{phylum}\t{classt}\t{order}\t{family}\t{genus}\t{species}\t{full_sequence}\n"
            )
        else:
            skipped_counter["insufficient_columns"] += 1
            print(
                f"mozaico WARNING: Incorrect header format: {header}. Check if the header has been through harmonization."
            )
    else:
        skipped_counter["unexpected_format"] += 1
        print(f"mozaico WARNING: Unexpected header format: {header}")

    return skipped_counter


def create_multibarcode_input_for_bmi(results_folder, output_file):
    insert_folder = results_folder + "/insert/filtered"
    pbs_incomplete_folder = (
        results_folder + "/incomplete_pbs/filtered/filtered_intersection"
    )

    try:
        fasta_files = glob.glob(os.path.join(insert_folder, "*.fasta")) + glob.glob(
            os.path.join(pbs_incomplete_folder, "*.fasta")
        )
        if not fasta_files:
            raise FileNotFoundError(
                "mozaiko ERROR: No FASTA files found in the specified folders."
            )

        # Initialize counters
        skipped_counter = {"insufficient_columns": 0, "unexpected_format": 0}
        processed_counter = 0

        with open(output_file, "w") as tsv_file:
            for fasta_path in fasta_files:
                primer_name = os.path.splitext(os.path.basename(fasta_path))[0]
                with open(fasta_path, "r") as fasta:
                    current_header = None
                    current_sequence = []
                    for line in fasta:
                        line = line.strip()
                        if line.startswith(">"):
                            if current_header and current_sequence:
                                skipped_counter = process_bmi_sequences(
                                    current_header,
                                    current_sequence,
                                    primer_name,
                                    tsv_file,
                                    skipped_counter,
                                )
                                processed_counter += 1
                            current_header = line[1:]
                            current_sequence = []
                        elif line:
                            current_sequence.append(line)
                    if current_header and current_sequence:
                        skipped_counter = process_bmi_sequences(
                            current_header,
                            current_sequence,
                            primer_name,
                            tsv_file,
                            skipped_counter,
                        )
                        processed_counter += 1

        # Print summary at the end
        print(f"mozaiko INFO: MultiBarcodeTools file created to {output_file}")
        print(f"mozaiko INFO: Total entries processed: {processed_counter}")
        print(
            f"mozaiko WARNING: {skipped_counter['insufficient_columns']} entries were not written because they did not meet the minimum length of columns. Check if headers have been through harmonization."
        )
        print(
            f"mozaiko WARNING: {skipped_counter['unexpected_format']} entries were not written because of unexpected header format."
        )

    except Exception as e:
        print(f"mozaiko ERROR: Error creating MultiBarcodeTools input - {str(e)}")
