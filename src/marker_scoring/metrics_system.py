import json
import os
import sys
from collections import defaultdict
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction

from src.in_silico_analysis.amplification import InSilicoAmplification
from src.marker_scoring.scoring_utils import *
from src.reference_database.sequence_import import CustomFastaImport


class OtlHandler:
    def __init__(self, otl=None, fasta=None):
        self.otl = otl
        self.fasta = fasta

    def validate_otl(self, otl=None):
        """
        This method validates the inputed Operational Taxonomic List (OTL).
        """
        if otl is not None:
            self.otl = otl

        if not os.path.exists(self.otl):
            print("mozaiko INFO: The OTL does not exist. Exiting...")
            sys.exit(1)

        _, file_extension = os.path.splitext(self.otl)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "tsv":
            print("mozaiko ERROR: The OTL must be a TSV file. Exiting...")
            sys.exit(1)

        otl_table = pd.read_csv(self.otl, sep="\t", header=0)

        otl_table_fields = otl_table.columns.tolist()

        if "taxa" not in otl_table_fields:
            print(f"mozaiko ERROR: The OTL must contain a column labeled 'taxa'.")
            sys.exit(1)

        otl_table = otl_table.dropna(subset=["taxa"])
        print(otl_table)

        self.otl = otl_table

    def import_otl(self):
        """
        This method retrieves a set of unique taxa from an OTL.
        """
        if self.otl is None:
            print(
                "mozaiko INFO: To continue the evaluation, a Operational Taxonomic List (OTL) is \
                    required. An OTL is a list contaning information on the taxonomic numenclature \
                        of all identifiable taxa in routine biomonitoring initiatives."
            )
            self.otl = input("Please enter the path to the OTL: ")

        self.validate_otl()

        otl = self.otl
        unique_otl_taxa = set()

        for entry in otl["taxa"]:
            unique_otl_taxa.add(entry)

        total_taxa_count = len(unique_otl_taxa)

        self.total_taxa = total_taxa_count
        self.otl_taxa_set = unique_otl_taxa

        return total_taxa_count, unique_otl_taxa

    def filter_fasta_for_species_not_in_otl(
        self, fasta_file, otl_taxa_set: set, overwrite: bool = True
    ):
        """
        This method filter a FASTA based on the set of unique taxa retrieved from the OTL. Entries
        whose taxon is not present in the OTL are removed.

        Parameters:
        input_file: Path to input FASTA file
        taxa_set: set
            Set of taxa names to keep
        overwrite: bool, default=True
            If True, overwrites the input file.
            If False, creates new file with '_filtered' suffix

        Returns:
        tuple
            (number of sequences kept, path to output file)
        """
        input_folder = os.path.dirname(fasta_file)
        file_name = os.path.basename(fasta_file)

        if overwrite:
            output_file = fasta_file
        else:
            filtered_folder = os.path.join(
                input_folder, f"{os.path.basename(input_folder)}_otl_filtered"
            )
            os.makedirs(filtered_folder, exist_ok=True)
            output_file = os.path.join(filtered_folder, file_name)

        total_seq_count: int = 0
        kept_seq_count: int = 0
        current_header: str = ""
        keep_sequence: bool = False
        sequences_to_write: List[str] = []
        current_sequence: List[str] = []

        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header and keep_sequence:
                        sequences_to_write.extend(
                            [current_header, "".join(current_sequence)]
                        )

                    current_header = line
                    current_sequence = []
                    total_seq_count += 1

                    try:
                        if "|" not in line:
                            print(f"mozaico WARNING: No '|' found in header - {line}")
                            keep_sequence = False
                            continue

                        header_parts = line.split(" | ")
                        print(header_parts)
                        if len(header_parts) < 2 or not header_parts[1].strip():
                            print(
                                f"mozaico WARNING: Taxonomy seems to not be present for - {line}"
                            )

                        taxa = header_parts[1].strip()
                        print(taxa)
                        print(otl_taxa_set)

                        if len(header_parts) > 1:
                            taxa = header_parts[1].strip()
                            keep_sequence = taxa in otl_taxa_set
                            print(keep_sequence)
                            if keep_sequence:
                                kept_seq_count += 1
                        else:
                            print(f"mozaico WARNING: Header parsing error - {line}")
                            keep_sequence = False
                    except Exception as e:
                        print(f"mozaico WARNING: Header parsing error - {line}")
                        print(f"Error details: {str(e)}")
                        keep_sequence = False
                elif line and keep_sequence:
                    current_sequence.append(line)

        if current_header and keep_sequence:
            sequences_to_write.extend([current_header, "".join(current_sequence)])

        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        if overwrite:
            temp_file = output_file + ".temp"
            with open(temp_file, "w") as f:
                for line in sequences_to_write:
                    f.write(f"{line}\n")
            os.replace(temp_file, output_file)
        else:
            with open(output_file, "w") as f:
                for line in sequences_to_write:
                    f.write(f"{line}\n")

        return total_seq_count, kept_seq_count, output_file

    def apply_fasta_filtering_for_all_fasta_files_in_folder(
        self, fasta_folder, otl_taxa_set, overwrite: bool = False
    ):
        """
        Applies the "filter_fasta_for_species_not_in_otl" method to all FASTA files in a folder.

        Parameters:
        fasta_folder: str
            Path to the folder containing FASTA files.
        otl_taxa_set: set
            Set of taxa names to keep.
        overwrite: bool, default=True
            If True, overwrites the input files.
            If False, creates new files with '_filtered' suffix.

        Returns:
        List[Tuple[str, int, int]]
            A list of tuples containing the file name, number of sequences processed, and number of sequences kept.
        """
        results = []

        for file_name in os.listdir(fasta_folder):
            file_path = os.path.join(fasta_folder, file_name)

            if file_name.endswith((".fasta")):
                print(f"Processing file: {file_name}")
                try:
                    total_seq_count, kept_seq_count, output_file = (
                        self.filter_fasta_for_species_not_in_otl(
                            fasta_file=file_path,
                            otl_taxa_set=otl_taxa_set,
                            overwrite=overwrite,
                        )
                    )
                    results.append((output_file, total_seq_count, kept_seq_count))
                except Exception as e:
                    print(
                        f"mozaico ERROR: Failed to process {file_name} due to {str(e)}"
                    )

        return results


class ReferenceDatabaseQuality:
    def __init__(self, all_inserts_path=None, otl=None):
        self.custom_fasta_import = CustomFastaImport()
        self.all_inserts_path = all_inserts_path

    def parse_fasta_files(self):
        """
        Parse all FASTA files in the given folder.
        Returns a dictionary with primer pairs as keys and file paths as values.
        """
        if self.all_inserts_path is None:
            print("mozaiko INFO: Please provide a folder containing FASTA files.")
            all_inserts_path = input("Please enter the path to the folder: ")
            self.all_inserts_path = str(all_inserts_path)

        fasta_files = {}
        for f in os.listdir(self.all_inserts_path):
            if f.endswith((".fasta")):
                primer_pair = os.path.splitext(f)[0]
                fasta_files[primer_pair] = os.path.join(self.all_inserts_path, f)

        return fasta_files

    def calculate_number_of_barcodes_per_taxon(self):
        """
        This method calculates the total number of barcodes that exists per taxon.

        Parameters:
        - all_inserts_file: File containing all inserts present in the reference database, either
        successfully amplified or not.

        Ouput:
        - barcodes_per_species: Dictionary containing information on how many barcodes (values)
        exists per species (keys).
        """
        fasta_files = self.parse_fasta_files()
        barcodes_per_species = {}

        for primer_pair, file_path in fasta_files.items():
            fasta_data = self.custom_fasta_import.read_fasta(
                file_path, check_taxid=False
            )

            barcodes_per_species[primer_pair] = {}

            for taxa in fasta_data["taxa_info"].unique():
                unique_sequences = fasta_data[fasta_data["taxa_info"] == taxa][
                    "sequence"
                ].unique()

                barcodes_per_species[primer_pair][taxa.strip()] = len(unique_sequences)

        return barcodes_per_species

    def calculate_percentage_of_taxa_w_x_barcodes(
        self,
        total_taxa_count: int,
        barcode_threshold=1,
    ):
        """
        This method calculates the percentage of taxa with more than X barcodes.

        Parameters:
        - otl: Operational Taxonomy List. List of taxons considered in a country for biomonitoring
        purposes.
        - barcode_threshold: Number of barcodes to consider as a threshold. Only taxons with more
        than this value will be considered.

        Output:
        - percentage_of_taxa: float
            Decimal percentage of taxa with more than X barcodes
        """
        barcodes_per_species = self.calculate_number_of_barcodes_per_taxon()
        results = {}

        for primer_pair, taxa_data in barcodes_per_species.items():
            taxa_meeting_threshold = sum(
                1 for count in taxa_data.values() if count > barcode_threshold
            )

            percentage = (taxa_meeting_threshold * 100) / total_taxa_count
            results[primer_pair] = round(percentage, 2)

        return results

    def barcoded_taxa_ratio(self, total_taxa_count: int):
        """
        This method calculates the Ratio of Barcoded Taxa (RBT).

        Parameters:
        - barcoded_taxa_five_plus: float
            Percentage of taxa with more than five barcodes.
        - barcoded_taxa_one_plus: float
            Percentage of taxa with more than one barcode.

        Output:
        - results: Dict
            A dictionary containing the raatio of barcoded taxa and the percentage of taxa with more than five barcodes per primer pair.
        """

        barcoded_taxa_five_plus = self.calculate_percentage_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=5
        )
        barcoded_taxa_one_plus = self.calculate_percentage_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=1
        )

        results = {}

        for primer_pair in barcoded_taxa_five_plus.keys():
            percent_5plus = barcoded_taxa_five_plus[primer_pair]
            percent_1plus = barcoded_taxa_one_plus[primer_pair]

            ratio_barcoded_taxa = (
                percent_5plus / percent_1plus if percent_5plus > 0 else 0
            )

            results[primer_pair] = {
                "barcoded_taxa_five_plus": percent_5plus,
                "ratio_barcoded_taxa": round(ratio_barcoded_taxa, 2),
            }

        return results


class Binding:
    def __init__(self, number_of_mismatches=None):
        self.amplification_instance = InSilicoAmplification()
        if number_of_mismatches is None:
            number_of_mismatches = (
                self.amplification_instance.get_number_of_mismatches()
            )
        self.number_of_mismatches = number_of_mismatches
        self.processed_primers = {}

    def parse_files_with_same_extension_in_folders(self, folder_path_A, folder_path_B):
        """
        This method checks if there are any files with the same name in two different folders
        and, if present, returns the path to this set of files.
        """
        files_A = [
            f
            for f in os.listdir(folder_path_A)
            if f.endswith((".fasta", ".fa", ".txt"))
        ]
        files_B = [
            f
            for f in os.listdir(folder_path_B)
            if f.endswith((".fasta", ".fa", ".txt"))
        ]

        base_files_A = {os.path.splitext(f)[0]: f for f in files_A}
        base_files_B = {os.path.splitext(f)[0]: f for f in files_B}

        files_same_base_name = set(base_files_A.keys()) & set(base_files_B.keys())

        if not files_same_base_name:
            print(
                "mozaiko ERROR: No matching files found between the two folders. Check folder directory and content before re-running."
            )
            return

        matching_files = []
        for file_name in files_same_base_name:
            file_A = base_files_A[file_name]
            file_B = base_files_B[file_name]
            file_path_A = os.path.join(folder_path_A, file_A)
            file_path_B = os.path.join(folder_path_B, file_B)
            matching_files.append((file_path_A, file_path_B))

        return matching_files

    def get_primer_table(self, primer_table):
        """
        This method reads the primer table into memory for downstream tasks.
        """
        self.amplification_instance.validate_primer_table(primer_table)

        self.primer_table = pd.read_csv(primer_table, sep="\t", header=0)

        return self.primer_table

    def get_pbs_table(self, amplicon_file, insert_file):
        """
        This method reads the PBS table into memory for downstream tasks.
        """

        pbs_table = extract_primer_binding_sites(amplicon_file, insert_file)

        pbs_table["primer_name"] = insert_file

        return pbs_table

    def primer_pbs_analysis(
        self,
        amplicon_folder,
        insert_folder,
        primer_table,
        save_results: bool = False,
    ) -> Tuple[Optional[Dict[str, pd.DataFrame]], Optional[pd.DataFrame]]:
        """
        This method analyzes Primer and PBS sequences to compute mismatch metrics, melting temperatures,
        and GC fractions for both the forward and reverse sequences.

        The implementation retrieves Primer-PBS information, processes amplicon and insert data files,
        and calculates various metrics related to primer-template binding. It generates a comprehensive
        DataFrame for each primer pair, including full-length mismatches, three-end mismatches, GC
        clamp matches, and melting temperature averages. It also computes GC fractions for forward and
        reverse primers.

        Parameters:
        - amplicon_folder: str
            Path to the folder containing amplicon sequence files from the in-silico amplification process.
        - insert_folder: str
             Path to the folder containing insert sequence files from the in-silico amplification process.
        - primer_table: str
            Path to the table containing the primer pair sequences and details.
        - save_results: bool
            If True, saves the comprehensive DataFrame for each primer pair and the primer GC fraction
        DataFrame as CSV files (default is False).

        Output:
        Tuple[Dict[str, pd.DataFrame], pd.DataFrame]
        - primer_pbs_df : dict
            A dictionary where keys are primer pair identifiers (based on barcode region and assay name),
            and values are DataFrames containing:
            - `seq_id`: Sequence identifier.
            - `taxon`: Taxonomic information.
            - `full_len_mismatch_sum`: Sum of mismatches for the forward and reverse primers across the
              entire primer-binding site.
            - `three_end_mismatch_sum`: Sum of mismatches for the last 5 bases of the primer-binding site
              (three-end region).
            - `min_tm`: Min. melting temperature between the forward and reverse primer-binding sites.
            - `gc_matches_sum`: Total GC matches in the three-end region for the forward and reverse primers.
        - primer_gc_df : pd.DataFrame
            A DataFrame summarizing GC fractions for each primer, with columns:
            - `barcode_region`: Barcode region of the primer.
            - `assay_name`: Assay name for the primer.
            - `forward_primer_gc_fraction`: GC fraction of the forward primer.
            - `reverse_primer_gc_fraction`: GC fraction of the reverse primer.
        """
        self.primer_table = self.get_primer_table(primer_table)
        matching_files = self.parse_files_with_same_extension_in_folders(
            amplicon_folder, insert_folder
        )

        if not matching_files:
            print("mozaiko ERROR: No matching files found in the provided folders.")
            return None, None

        primer_pbs_df = {}
        primer_gc_fractions = []

        for _primer_ind, primer_row in self.primer_table.iterrows():
            barcode_region = primer_row["barcode_region"]
            assay_name = primer_row["assay_name"]
            pbs_filename = f"{barcode_region}_{assay_name}"
            primer_seq_fwd = primer_row["fwd_seq"]
            primer_seq_rev = primer_row["rev_seq"]
            rev_comp_primer_seq_rev = str(Seq(primer_seq_rev).reverse_complement())

            # Compute GC fraction for Primer
            primer_gc_fractions.append(
                {
                    "barcode_region": barcode_region,
                    "assay_name": assay_name,
                    "forward_primer_gc_fraction": gc_fraction(primer_seq_fwd),
                    "reverse_primer_gc_fraction": gc_fraction(primer_seq_rev),
                }
            )

            full_mismatches_data = []
            three_end_mismatches_data = []
            pbs_melting_temperature_data = []
            three_end_gc_matches_data = []

            # Process matching files (same primer name)
            matching_files_found = False

            for amplicon_file, insert_file in matching_files:
                amplicon_filename = os.path.splitext(os.path.basename(amplicon_file))[0]

                if pbs_filename == amplicon_filename:
                    matching_files_found = True
                    pbs_table = self.get_pbs_table(amplicon_file, insert_file)

                    for _, pbs_row in pbs_table.iterrows():
                        seq_header = pbs_row["header"].replace(">", "")
                        seq_id = seq_header.split("|")[0].replace(" ", "")
                        taxon = seq_header.split("|")[1]
                        pbs_fwd_seq = pbs_row["fwd_seq"]
                        pbs_rev_seq = pbs_row["rev_seq"]

                        # Melting Temperature
                        pbs_fwd_tm = MeltingTemp.Tm_GC(
                            pbs_fwd_seq, valueset=7, strict=False
                        )
                        pbs_rev_tm = MeltingTemp.Tm_GC(
                            pbs_rev_seq, valueset=7, strict=False
                        )

                        min_tm = min(pbs_fwd_tm, pbs_rev_tm)
                        min_tm = float(round(min_tm, 2))
                        substraction = pbs_fwd_tm - pbs_rev_tm
                        delta_tm = abs(substraction)
                        delta_tm = round(delta_tm, 2)

                        pbs_melting_temperature_data.append(
                            {
                                "seq_id": seq_id,
                                "taxon": taxon,
                                "min_tm": min_tm,
                                "delta_tm": delta_tm,
                            }
                        )

                        # Compute full Primer-Template mismatches
                        full_fwd_mismatches = calculate_iupac_mismatches(
                            primer_seq_fwd, pbs_fwd_seq
                        )
                        full_rev_mismatches = calculate_iupac_mismatches(
                            rev_comp_primer_seq_rev, pbs_rev_seq
                        )

                        full_mismatches_data.append(
                            {
                                "seq_id": seq_id,
                                "taxon": taxon,
                                "full_len_mismatch_sum": full_fwd_mismatches
                                + full_rev_mismatches,
                            }
                        )

                        # Compute three-end mismatches and GC clamp
                        three_end_fwd_seq = pbs_fwd_seq[-5:]
                        three_end_rev_seq = pbs_rev_seq[-5:]
                        three_end_fwd_primers = primer_seq_fwd[-5:]
                        three_end_rev_primers = rev_comp_primer_seq_rev[-5:]

                        three_end_fwd_mismatches = calculate_iupac_mismatches(
                            three_end_fwd_primers, three_end_fwd_seq
                        )
                        three_end_rev_mismatches = calculate_iupac_mismatches(
                            three_end_rev_primers, three_end_rev_seq
                        )

                        three_end_mismatches_data.append(
                            {
                                "seq_id": seq_id,
                                "taxon": taxon,
                                "three_end_mismatch_sum": three_end_fwd_mismatches
                                + three_end_rev_mismatches,
                            }
                        )

                        # GC matches at three-end
                        three_end_fwd_gc_matches = calculate_iupac_mismatches(
                            three_end_fwd_primers,
                            three_end_fwd_seq,
                            search_gc_clamp=True,
                        )[1]
                        three_end_rev_gc_matches = calculate_iupac_mismatches(
                            three_end_rev_primers,
                            three_end_rev_seq,
                            search_gc_clamp=True,
                        )[1]

                        three_end_gc_matches_data.append(
                            {
                                "seq_id": seq_id,
                                "taxon": taxon,
                                "gc_matches_fwd": three_end_fwd_gc_matches,
                                "gc_matches_rev": three_end_rev_gc_matches,
                            }
                        )

            if matching_files_found:
                # Start by adding one df
                comprehensive_df = pd.DataFrame(full_mismatches_data)
                # Merge all other df's by seq_id & taxon
                comprehensive_df = comprehensive_df.merge(
                    pd.DataFrame(three_end_mismatches_data)[
                        ["seq_id", "taxon", "three_end_mismatch_sum"]
                    ],
                    on=["seq_id", "taxon"],
                )

                comprehensive_df = comprehensive_df.merge(
                    pd.DataFrame(pbs_melting_temperature_data)[
                        ["seq_id", "taxon", "min_tm", "delta_tm"]
                    ],
                    on=["seq_id", "taxon"],
                )

                comprehensive_df = comprehensive_df.merge(
                    pd.DataFrame(three_end_gc_matches_data)[
                        ["seq_id", "taxon", "gc_matches_fwd", "gc_matches_rev"]
                    ],
                    on=["seq_id", "taxon"],
                )

                column_order = [
                    "seq_id",
                    "taxon",
                    "full_len_mismatch_sum",
                    "three_end_mismatch_sum",
                    "gc_matches_fwd",
                    "gc_matches_rev",
                    "min_tm",
                    "delta_tm",
                ]

                comprehensive_df = comprehensive_df[column_order]

                primer_pbs_df[pbs_filename] = comprehensive_df

                if save_results:
                    comprehensive_df.to_csv(
                        f"{pbs_filename}_comprehensive.csv", index=False
                    )

        primer_gc_df = pd.DataFrame(primer_gc_fractions)

        if save_results:
            primer_gc_df.to_csv("primer_gc_fractions.csv", index=False)

        return primer_pbs_df, primer_gc_df

    def add_missing_otl_taxa_to_df_with_values_of_zero(self, primer_df, otl_taxa_set):
        """
        This methid adds entries for taxons that are not present in the in-silico analysis
        but need to be present for downstream analysis regarding the OTL.

        Parameters:
        - dataframe: Dataframe
            A dataframe containing the primer-PBS analysis.
        - otl_taxa_set: set
            A set containing all unique values of taxa present in the OTL.

        Returns:
        - otl_populated_df: Dataframe
            A dataframe equal to input but with additional entries of taxa that was in the OTL but
            not on the in-silico files. These entries will have 'otl-import' for 'seq_id' and
            a value of zero for all other analysis.
        """
        taxon_on_df = set(primer_df["taxon"].unique())
        missing_taxon = otl_taxa_set - taxon_on_df
        intersection = taxon_on_df & missing_taxon
        len_intersection = len(intersection)

        if missing_taxon:
            new_entries = pd.DataFrame(
                {
                    "taxon": list(missing_taxon),
                    "seq_id": ["otl-import"] * len(missing_taxon),
                    **{
                        col: 0
                        for col in primer_df.columns
                        if col not in ["taxon", "seq_id"]
                    },
                }
            )

        otl_populated_df = pd.concat([primer_df, new_entries], ignore_index=True)

        return otl_populated_df

    def iterate_over_primer_pbs_df(
        self,
        primer_pbs_df,
        add_otl_taxa: bool = True,
        otl_taxa_set: Optional[set] = None,
    ):
        """
        Method to iterate over the comprehensive primer Dataframe to extract a Dataframe for each
        primer and attribute it to a variable.

        Parameter:
        - primer_pbs_df: Dict
            A dictionary containing the primers as keys and Dataframes as values.

        Return:
        - processed_primers: List
            A list containing the primer's dataframes.
        """
        processed_primers = []

        for key, dataframe in primer_pbs_df.items():
            variable_name = key.replace("-", "_").replace(" ", "_")
            self.processed_primers[variable_name] = dataframe
            processed_primers.append(variable_name)

        if add_otl_taxa == True:
            otl_populated_dfs = [
                self.add_missing_otl_taxa_to_df_with_values_of_zero(
                    self.processed_primers[primer_df], otl_taxa_set
                )
                for primer_df in processed_primers
            ]
            self.processed_primers = {
                key: df for key, df in zip(processed_primers, otl_populated_dfs)
            }

        return self.processed_primers

    def get_primer_df(self, primer_name):
        """
        This method retrieves a DataFrame by its processed primer name.

        Parameter:
        - primer_name: str
            The name of the primer to retrieve.

        Returns:
        - DataFrame corresponding to the primer.
        """
        return self.processed_primers.get(primer_name, None)

    def process_analysis_per_taxon(
        self,
        primer_df: pd.DataFrame,
        operation: Literal["min", "max", "sum", "mean", "coef_var"],
        analysis_name: str,
    ) -> Union[pd.DataFrame, pd.Series]:
        """
        This method performs user-inputed operations on a groupby of 'taxon'.

        Parameters:
        - operation: Literal["min", "max", "sum", "mean", "coef_var"]
            The operation to perform.
        - analysis_name: str
            The column on which to perform the operation.

        Returns:
        - pd.DataFrame
            The grouped and processed DataFrame.
        """
        if "taxon" not in primer_df.columns:
            raise ValueError(
                "mozaiko ERROR: The input DataFrame must contain a 'taxon' column."
            )

        if analysis_name not in primer_df.columns:
            raise ValueError(
                f"mozaiko ERROR: The specified analysis '{analysis_name}' does not exist in the DataFrame."
            )

        grouped_taxa = primer_df.groupby("taxon")

        if operation in {"min", "max", "sum", "mean"}:
            result = getattr(grouped_taxa[analysis_name], operation)().astype(float)
        elif operation == "coef_var":
            mean = grouped_taxa[analysis_name].mean()
            print(mean)
            std = grouped_taxa[analysis_name].std()
            print(std)
            result = (std / mean.replace(0, pd.NA)) * 100

        else:
            raise ValueError(
                f"mozaiko ERROR: Unrecognized operation: '{operation}'. "
                "Please choose from 'min', 'max', 'sum', 'mean', or 'coef_var'."
            )

        result = pd.DataFrame(result)

        return result

    def process_analysis_across_taxon(
        self,
        tax_grouped_df: pd.DataFrame,
        operation: Literal["min", "max", "sum", "mean", "coef_var"],
    ):
        """
        This method will perform operation across the taxon to determine a value for a primer pair.

        Parameters:
        - tax_grouped_df: pd.DataFrame
            DataFrame containing grouped taxa data.
        - operation: Literal["min", "max", "sum", "mean", "coef_var"]
            The operation to apply.

        Returns:
        - value: float
            A value that reflects a characteristic of the primer set.
        """
        if tax_grouped_df.empty:
            raise ValueError("mozaiko ERROR: Input DataFrame is empty.")

        if operation in {"min", "max", "sum", "mean"}:
            result = getattr(tax_grouped_df, operation)().astype(float)
        elif operation == "coef_var":
            mean = tax_grouped_df.mean()
            std = tax_grouped_df.std()
            result = (std / mean.replace(0, pd.NA)) * 100
        else:
            raise ValueError(
                f"mozaiko ERROR: Unrecognized operation: '{operation}'. "
                "Please choose from 'min', 'max', 'sum', 'mean', or 'coef_var'."
            )

        return round(float(result.iloc[0]), 2)

    def get_priming_ratio(
        self, max_mismatch_full_len: pd.DataFrame, max_mismatch_three_end: pd.DataFrame
    ):
        """
        This method computes the priming ratio between the maximum number of mismatches per taxon
        and the maximum number of mismatches at the 3'end per taxon. This metric is performed at the
        primer-level, thefore, all of the ratio are then summed.

        Parameters:
        - max_mismatch_full_len: DataFrame
            A DataFrame containing the max number of mismatches attained for each taxon,
            across the total sequence length.
        - max_mismatch_three_end: DataFrame
            A DataFrame containing the max number of mismatches attained for each taxon,
            across the 3'end length.

        Output:
        - priming_ratio: float
            The value for the priming ratio.
        """
        merged_df = max_mismatch_three_end.join(max_mismatch_full_len, how="inner")
        merged_df["priming_ratio"] = (
            merged_df["three_end_mismatch_sum"] / merged_df["full_len_mismatch_sum"]
        )

        ratio_col = merged_df["priming_ratio"]
        ratio_df = pd.DataFrame(ratio_col)
        priming_ratio = self.process_analysis_across_taxon(ratio_df, operation="sum")

        priming_ratio = float(round(priming_ratio, 2))

        return priming_ratio

    def get_total_gc_matches(self, primer_pbs_df: pd.DataFrame):
        """
        This method retrived a score for the number of GCC matches between the PBS-Primer. The
        presence of G and C bases at the 3' end of both forward and reverse primers (GC clamp)
        is an indicator of higher binding stability. However, more than 3 G's or C's should be
        avoided in the last 5 bases at the 3' end of the primer. Given this, a score system is
        applied to the number of GC matches between the Primer-PBS:
        - 1-3 matches -> 1 point
        - 0 matches -> 0 points
        - >3 matches -> 0 points

        The method translates the GC matches to points for both the forward and reverse sequences
        and sums the result for each Seq Id.

        Parameters:
        - primer_pbs_df: DataFrame
            A dictionary containing the primers as keys and Dataframes as values.

        Return:
        -
        """

        def gc_scoring(value):
            if value == 0:
                value = 0
            elif value > 3:
                value = 0
            elif 1 <= value <= 3:
                value = 1
            return value

        primer_pbs_df["gc_matches_fwd"] = primer_pbs_df["gc_matches_fwd"].apply(
            gc_scoring
        )
        primer_pbs_df["gc_matches_rev"] = primer_pbs_df["gc_matches_rev"].apply(
            gc_scoring
        )

        primer_pbs_df["gc_matches_score"] = (
            primer_pbs_df["gc_matches_fwd"] + primer_pbs_df["gc_matches_rev"]
        )

        return primer_pbs_df

    def tm_score(self, primer_pbs_df: pd.DataFrame, temp_threshold: float = 2.0):
        """
        This method retrieves the percentage of entries whose variation of temperature between
        the forward and reverse PBS sequence is lower than the temp_threshold.

        Parameters:
        - primer_pbs_df: Dataframe
            A dataframe containing the results from the comprehensive Primer-PBS analysis.
        - temp_theresold: float
            The temperature threshold to look for.

        Return:
        - tm_score: int
        """
        delta_col = primer_pbs_df["delta_tm"]
        number_of_entries_passing_threshold = (delta_col < temp_threshold).sum()
        total_count = delta_col.count()

        tm_score = (number_of_entries_passing_threshold / total_count) * 100

        tm_score = float(round(tm_score, 2))

        return tm_score

    def count_unique_taxa(self, fasta_file):
        unique_taxa = set()

        for record in SeqIO.parse(fasta_file, "fasta"):
            parts = record.description.split("|")
            if len(parts) > 1:
                taxonomy = parts[1].strip()
                unique_taxa.add(taxonomy)

        return len(unique_taxa)

    def get_outputs_taxa_counts(self, results_folder):
        """
        This method computes the total number of taxa that were successfuly amplified in-silico,
        the number of taxa that contain PBS, even if amplification was not successful, and the number
        of taxa that did not contain PBS.
        A taxa is amplified if at least one sequence identified with the nomenclature is kept after
        the in-silico process.

        Parameters:
        - results_folder: str
            Path to the folder containing the results from the amplification process and its
             subdirectories. Subdirectories names should not be changed.

        Return:
        - amplified_taxa_count: int
            Total number of taxa that were successfuly amplified
        """
        # Input A
        in_silico_amplified_inserts = results_folder + "/all_inserts"
        # Input B
        all_inserts_with_pbs = results_folder + "/all_complete_pbs/filtered"
        # Input C
        inserts_with_incomplete_pbs = results_folder + "/incomplete_pbs/filtered"

        folder_list = [
            ("taxa_in_silico_amplified", in_silico_amplified_inserts),
            ("taxa_with_pbs", all_inserts_with_pbs),
            ("taxa_with_incomplete_pbs", inserts_with_incomplete_pbs),
        ]
        data = {}

        for folder_name, folder_path in folder_list:
            data[folder_name] = {}
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    if file.endswith(".fasta"):
                        primer_name = os.path.splitext(file)[0]
                        file_path = os.path.join(root, file)

                        unique_taxa_count = self.count_unique_taxa(file_path)

                        data[folder_name][primer_name] = unique_taxa_count

        insert_taxa_counts_df = pd.DataFrame(data).fillna(0).astype(int)

        insert_taxa_counts_df = insert_taxa_counts_df

        insert_taxa_counts_df["taxa_with_pbs"] = (
            insert_taxa_counts_df["taxa_in_silico_amplified"]
            + insert_taxa_counts_df["taxa_with_pbs"]
        )

        return insert_taxa_counts_df

    def calculate_amplification_success_score(self, results_folder):
        """
        This method calculates the percentage of successfully amplified taxa, computed by divinding
        the number of in-silico amplified taxa by taxa with PBS.
        """
        insert_taxa_counts_df = self.get_outputs_taxa_counts(results_folder)

        amplification_score = (
            insert_taxa_counts_df["taxa_in_silico_amplified"]
            / insert_taxa_counts_df["taxa_with_pbs"]
            * 100
        )
        insert_taxa_counts_df["amplification_sucess_percent"] = round(
            amplification_score, 2
        )

        return insert_taxa_counts_df


class MetricsSystemExecutor:
    """
    This class orchestrates the entire Metrics System, coordinating the execution of all evaluation
    categories by calling the necessary classes and methods.
    """

    def __init__(self, all_inserts_folder=None, otl=None):
        # Initialize Reference Database Quality Category
        self.ref_db = ReferenceDatabaseQuality()
        # Initialize OTL and related variables
        self.otl = otl
        self.all_inserts_folder = all_inserts_folder
        self.otl_handler = OtlHandler(self.otl)
        self.otl_handler.import_otl()
        self.total_otl_taxa_count = self.otl_handler.total_taxa
        self.otl_unique_taxa_set = self.otl_handler.otl_taxa_set
        # Filter FASTA files per OTL species
        self.otl_handler.apply_fasta_filtering_for_all_fasta_files_in_folder(
            fasta_folder=self.all_inserts_folder,
            otl_taxa_set=self.otl_unique_taxa_set,
            overwrite=False,
        )
        self.filtered_inserts_folder = os.path.join(
            self.all_inserts_folder,
            os.path.basename(self.all_inserts_folder) + "_otl_filtered",
        )

    def calculate_reference_database_quality(self):
        """
        This method processed the Reference Database Quality evaluation.
        It calculates the Barcode Coverage Score (BCS) for each primer and stores the results.

        Output:
        - ref_bd_scores: dict
            A dictionary where the key is the primer name and the value is a tuple containing
            the percentage of taxa with more than five barcodes and the rounded Ratio of Barcoded
            Taxa.
        """

        cls = ReferenceDatabaseQuality(self.filtered_inserts_folder, self.otl)
        reference_db_quality = cls.barcoded_taxa_ratio(
            total_taxa_count=self.total_otl_taxa_count
        )

        return reference_db_quality

    # TODO: check how to best do folder and file intake, otl handling
