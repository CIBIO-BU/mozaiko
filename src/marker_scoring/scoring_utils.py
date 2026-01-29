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

from src.reference_database.sequence_import import *
from src.marker_scoring.metrics_system import *
import pandas as pd

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


def sequence_count_tracking(original_database, analysis_folder, save_results: bool = True):
    """
    Tracks sequence count per analysis step.

    Parameters:
    - original_database: Path to the original input FASTA database.
    - analysis_folder: Path to folder containing the analysis outcomes.
    - save_results: If True, saves the tracking results to a TSV file.

    Returns:
    - A pivoted DataFrame showing sequence counts per primer and analysis step.
    """

    try:
        analysis_folder = Path(analysis_folder)
        rename_analysis_folder(analysis_folder)

        # Define analysis folders
        target_folders = ["all_complete_pbs", "amplicon", "incomplete_pbs", "insert", "input_B"]
        renamed_folders = ["input_ABC", "input_AC", "input_A"]

        step_files = {}

        # Walk through analysis folder to find files
        for root, dirs, files in os.walk(analysis_folder):
            path_parts = root.split(os.path.sep)
            for target in target_folders + renamed_folders:
                if target in path_parts:
                    target_index = path_parts.index(target)
                    step_name = "-".join(path_parts[target_index:])
                    for file in files:
                        if file.endswith((".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn")):
                            file_path = os.path.join(root, file)
                            step_files.setdefault(step_name, []).append(file_path)
                    break

        # Count original database sequences and print
        original_count = 0
        if os.path.exists(original_database):
            try:
                count_result = subprocess.run(
                    ["grep", "-c", "^>", original_database],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                original_count = int(count_result.stdout.strip())
                print(f"mozaiko INFO: Number of sequences in original database: {original_count}")
            except (subprocess.CalledProcessError, ValueError) as e:
                print(f"mozaiko WARNING: Could not count sequences in original database: {e}")
        else:
            print("mozaiko WARNING: Original database file not found.")

        primer_step_counts = {}

        # Count sequences using grep for analysis steps only
        for step_name, file_paths in step_files.items():
            for file_path in file_paths:
                try:
                    count_result = subprocess.run(
                        ["grep", "-c", "^>", file_path],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    count = int(count_result.stdout.strip())
                    file_base_name = os.path.basename(file_path).split('.')[0]

                    # Determine primer name by removing folder prefixes
                    primer_name = file_base_name
                    for prefix in [f"{folder}-" for folder in target_folders + renamed_folders]:
                        if primer_name.startswith(prefix):
                            primer_name = primer_name[len(prefix):]
                            break

                    key = (primer_name, step_name)
                    if key not in primer_step_counts or primer_step_counts[key] < count:
                        primer_step_counts[key] = count

                except subprocess.CalledProcessError as e:
                    print(f"mozaiko ERROR: Error processing file {file_path}: {e}")
                except ValueError as e:
                    print(f"mozaiko ERROR: Invalid sequence count for file {file_path}: {e}")

        # Convert to DataFrame
        records = [
            {"primer_name": primer, "analysis_step": step, "number_of_sequences": count}
            for (primer, step), count in primer_step_counts.items()
        ]
        df = pd.DataFrame(records)

        if df.empty:
            print("mozaiko WARNING: No data found to create sequence count tracking.")
            return None

        df_pivoted = df.pivot(index="primer_name", columns="analysis_step", values="number_of_sequences")

        # Column order: base folders, renamed folders, then filtered versions
        preferred_columns = target_folders + renamed_folders
        preferred_columns += [f"{folder}-filtered" for folder in target_folders + renamed_folders]

        # Add remaining columns alphabetically
        remaining_columns = [col for col in df_pivoted.columns if col not in preferred_columns]
        final_columns = [col for col in preferred_columns if col in df_pivoted.columns] + sorted(remaining_columns)

        df_pivoted = df_pivoted[final_columns].fillna(0)

        # Remove 'analysis_step' from column names (pivot artifact)
        df_pivoted.columns.name = None

        # Compute percentage metrics with corrected column names
        percentage_columns = []
        try:
            # Check for the actual column names that should exist
            required_cols = ["all_complete_pbs-input_ABC", "insert-input_A", "input_B"]
            missing_cols = [col for col in required_cols if col not in df_pivoted.columns]

            if not missing_cols:
                # Calculate percentage of sequences with PBS
                # (sequences with incomplete PBS + sequences with complete PBS) / total sequences with PBS
                total_pbs_sequences = df_pivoted["input_B"] + df_pivoted["insert-input_A"]
                total_sequences_with_pbs = df_pivoted["all_complete_pbs-input_ABC"]

                df_pivoted["percentage_of_sequences_with_PBS"] = (
                    total_pbs_sequences / total_sequences_with_pbs * 100
                ).round(2)

                # Calculate percentage of amplified sequences among those with PBS
                # sequences with complete PBS / (sequences with incomplete PBS + sequences with complete PBS)
                df_pivoted["percentage_of_amplified_sequences_with_PBS"] = (
                    df_pivoted["insert-input_A"] / total_pbs_sequences * 100
                ).round(2)

                percentage_columns = ["percentage_of_sequences_with_PBS", "percentage_of_amplified_sequences_with_PBS"]
                print("mozaiko INFO: Percentage calculations completed successfully.")

            else:
                print(f"mozaiko WARNING: Cannot compute percentages. Missing columns: {missing_cols}")

        except Exception as e:
            print(f"mozaiko WARNING: Failed to compute percentages: {e}")
            import traceback
            traceback.print_exc()

        # Insert percentage columns at the end
        if percentage_columns:
            final_columns.extend(percentage_columns)

        # Apply the final column order
        available_columns = [col for col in final_columns if col in df_pivoted.columns]
        final_df = df_pivoted[available_columns]

        # Save result if requested
        if save_results:
            run_name = analysis_folder.name
            output_path = analysis_folder / f"{run_name}-sequence_count_track.tsv"
            final_df.to_csv(output_path, sep="\t", index=True, header=True, index_label="primer_name")
            print(f"mozaiko INFO: Sequence count tracking file saved to {output_path}.")

        return final_df

    except Exception as e:
        print(f"mozaiko ERROR: Unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        return None


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

def rename_analysis_folder(analysis_folder):
    """
    Rename the in silico analysis folder to respective input names.

    Parameters:
    - analysis_folder: Path to the analysis folder
    """
    from pathlib import Path

    analysis_folder = Path(analysis_folder)

    # Input ABC
    filtered_folder1 = analysis_folder/ "all_complete_pbs" / "filtered"
    renamed_folder1 = analysis_folder / "all_complete_pbs" / "input_ABC"

    if filtered_folder1.exists() and not renamed_folder1.exists():
        filtered_folder1.rename(renamed_folder1)

    # Input AC
    filtered_folder2 = analysis_folder/ "incomplete_pbs" / "filtered"
    renamed_folder2 = analysis_folder / "incomplete_pbs" / "input_AC"

    if filtered_folder2.exists() and not renamed_folder2.exists():
        filtered_folder2.rename(renamed_folder2)

    # Input A
    filtered_folder3 = analysis_folder/ "insert" / "filtered"
    renamed_folder3 = analysis_folder / "insert" / "input_A"

    if filtered_folder3.exists() and not renamed_folder3.exists():
        filtered_folder3.rename(renamed_folder3)


def compute_pbs_stats(input_A_file, input_B_file, input_ABC_file, otl, output_dir):
    """
    Process input files and compute taxonomic coverage statistics.
    """
    input_A_df = _process_fasta_file(input_A_file, output_dir)
    input_B_df = _process_fasta_file(input_B_file, output_dir)
    input_ABC_df = _process_fasta_file(input_ABC_file, output_dir)

    # Create A+B input
    merged_AB = pd.concat([input_A_df, input_B_df], ignore_index=True).drop_duplicates()

    # Calculate coverage metrics
    AB_tax_coverage = _calculate_tax_coverage(merged_AB, otl)
    A_tax_coverage = _calculate_tax_coverage(input_A_df, otl)
    ABC_tax_coverage = _calculate_tax_coverage(input_ABC_df, otl)

    # Compute percentages
    prcnt_pbs = round((AB_tax_coverage / ABC_tax_coverage) * 100, 2) if ABC_tax_coverage > 0 else 0
    prcnt_amplified_with_pbs = round((A_tax_coverage / AB_tax_coverage) * 100, 2) if AB_tax_coverage > 0 else 0

    return {
        'AB_tax_coverage': float(AB_tax_coverage),
        'A_tax_coverage': float(A_tax_coverage),
        'ABC_tax_coverage': float(ABC_tax_coverage),
        'prcnt_pbs': prcnt_pbs,
        'prcnt_amplified_with_pbs': prcnt_amplified_with_pbs
    }


def _process_fasta_file(input_file, output_dir):
    """
    Convert FASTA file to DataFrame with taxonomic information.
    """
    file_name = input_file.split('/')[-1].strip('.fasta')
    step_name = input_file.split('/')[-2]
    output_file_name = file_name + '_' + step_name + '.tsv'
    output_file = os.path.join(output_dir, output_file_name)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.path.exists(output_file):
        df = pd.read_csv(output_file, sep='\t', usecols=['taxa_info'])
        df[['family', 'genus', 'species']] = df['taxa_info'].str.split('|', expand=True).iloc[:, -3:]
        return df[['family', 'genus', 'species']]

    custom_fasta_import = CustomFastaImport(input_file)
    custom_fasta_import.read_fasta(input_file, check_taxid=False)
    custom_fasta_import.df2csv(output_file)

    df = pd.read_csv(output_file, sep='\t', usecols=['taxa_info'])
    df[['family', 'genus', 'species']] = df['taxa_info'].str.split('|', expand=True).iloc[:, -3:]

    return df[['family', 'genus', 'species']]


def _calculate_tax_coverage(input_df, otl, cutff_val: int = 1):
    """
    Calculate taxonomic coverage percentage for a given dataset.
    """
    otl_handler = OtlHandler(otl)
    otl_handler.import_otl()
    total_taxa = otl_handler.total_taxa

    otl_df = otl_handler.otl[['family', 'genus', 'species', 'rank']].drop_duplicates().reset_index(drop=True)

    # Count matching taxa
    def count_taxa_matches(row):
        rank = row['rank']
        entry = row[rank]
        return (input_df[rank] == entry).sum()

    otl_df['tax_count'] = otl_df.apply(count_taxa_matches, axis=1)

    cutoff_num = (otl_df['tax_count'] >= cutff_val).sum()
    tax_coverage = round((cutoff_num / total_taxa) * 100, 2)
    tax_coverage = float(tax_coverage)

    return tax_coverage

def compute_pbs_stats_multiple_otls(input_A_file, input_B_file, input_ABC_file, otl_path, output_dir):
    for otl_file in os.listdir(otl_path):
        otl_file_p = os.path.join(otl_path, otl_file)
        country = otl_file.split('_')[0]
        country_stats = compute_pbs_stats(input_A_file, input_B_file, input_ABC_file, otl_file_p, output_dir)

        print(f"{country} stats: {country_stats}")

def barcoded_taxa(input_ABC_file, otl, output_dir):
    input_ABC_df = _process_fasta_file(input_ABC_file, output_dir)

    _calculate_tax_coverage(input_ABC_df, otl, cutff_val=1)