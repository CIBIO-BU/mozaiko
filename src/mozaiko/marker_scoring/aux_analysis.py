import numpy as np
import pandas as pd
import os
import subprocess
from pathlib import Path

from mozaiko.reference_database.sequence_import import CustomFastaImport
from mozaiko.marker_scoring.metrics_system import OtlHandler

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

        # add type annotation to step_files
        step_files: dict[str, list[str]] = {}

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

        primer_step_counts: dict[tuple[str, str], int] = {}

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

    return _calculate_tax_coverage(input_ABC_df, otl, cutff_val=1)