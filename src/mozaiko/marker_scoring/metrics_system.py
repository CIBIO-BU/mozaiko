"""
This module implements the metrics system for the mozaiko pipeline.

It includes classes and methods for handling Operational Taxonomic Lists (OTL),
calculating primer binding sites (PBS), and evaluating the quality of reference databases.

"""



import copy
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp, gc_fraction

from mozaiko.in_silico_analysis.amplification import InSilicoAmplification
from mozaiko.marker_scoring.scoring_utils import *
from mozaiko.reference_database.sequence_import import CustomFastaImport


class OtlHandler:
    def __init__(self, otl=None, fasta=None):
        self.otl = otl
        self.fasta = fasta
        self.taxa_hierarchy = {}
        self.fasta_handler = CustomFastaImport()

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

        self.otl = otl_table

    def pre_process_otl(self):
        """
        Pre-processed OTL to:
        1) Transform entries with '-' to NA.
        2) Tranform entries to lower case.
        3) Remove entries with with 'kingdom', 'phylum', 'class', 'order' in 'ranks'.
        5) Clean ASCII characters from the 'scientificName' column.

        Returns:
        - otl: DataFrame
            The pre-processed OTL.
        """
        # 1) Transform entries with '-' to NA.
        pd.set_option('future.no_silent_downcasting', True)
        self.otl = self.otl.replace("-", np.nan)

        # 2) Tranform entries to lower case.
        self.otl["rank"] = self.otl["rank"].str.lower()

        # 3) Normalize level above species ('form', 'variety', 'subspecies') to 'species'
        self.otl["rank"] = self.otl["rank"].replace(
            {"form": "species", "variety": "species", "subspecies": "species"}
        )

        # 4) Keep only entries with rank 'family', 'genus', or 'species'
        ranks_to_keep = ["family", "genus", "species"]
        self.otl = self.otl[self.otl["rank"].isin(ranks_to_keep)]

        # 4) Remove entries where 'scientificName' is NA.
        self.otl = self.otl.dropna(subset=["scientificName"])

        # 4) Clean ASCII characters from the 'scientificName' column.
        self.otl["scientificName"] = self.otl["scientificName"].apply(
            lambda x: self.fasta_handler.clean_header(x) if pd.notnull(x) else x
        )

        # 5) Dereplicate entries based on 'scientificName'
        self.otl = self.otl.drop_duplicates(subset=["scientificName"])

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
        self.pre_process_otl()
        self.create_taxonomic_hierarchy()

        otl = self.otl

        otl_taxa_count = otl[['family', 'genus', 'species']]
        otl_taxa_count = otl_taxa_count.replace(np.nan, "nan")
        otl_taxa_count = otl_taxa_count.drop_duplicates()

        total_taxa_count = otl_taxa_count.shape[0]

        self.total_taxa = total_taxa_count

        return total_taxa_count

    def create_taxonomic_hierarchy(self):
        """
        Creates a nested dictionary representing the taxonomic hierarchy from the OTL.
        """
        self.taxa_hierarchy = {}

        otl = self.otl

        # ensure we have the minimum-required columns for taxonomy
        required_cols = ["family", "genus", "species"]
        if not all(col in otl.columns for col in required_cols):
            raise ValueError(
                "mozaiko INFO: OTL must contain family, genus, and species columns."
            )

        # group by family and genus to create the hierarchy
        for _, row in otl.iterrows():
            family = row["family"] if pd.notna(row["family"]) else None
            genus = row["genus"] if pd.notna(row["genus"]) else None
            species = row["species"] if pd.notna(row["species"]) else None

            if family not in self.taxa_hierarchy:
                self.taxa_hierarchy[family] = {"genera": {}, "count": 0}

            if genus:
                if genus not in self.taxa_hierarchy[family]["genera"]:
                    self.taxa_hierarchy[family]["genera"][genus] = {
                        "species": {},
                        "count": 0,
                    }

                if species:
                    self.taxa_hierarchy[family]["genera"][genus]["species"][species] = {
                        "count": 0
                    }

        return self.taxa_hierarchy


class ReferenceDatabaseQuality:
    def __init__(self, otl, all_inserts_path=None):
        self.custom_fasta_import = CustomFastaImport()
        self.all_inserts_path = all_inserts_path
        self.otl_handler = OtlHandler(otl)
        self.otl_handler.import_otl()
        self.taxa_hierarchy = self.otl_handler.taxa_hierarchy

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

    def calculate_number_of_barcodes_per_fasta_entry(self):
        """
        Calculates the total number of barcodes per unique taxonomy entry.
        Uses pandas groupby for faster processing and handles non-string values.

        Returns:
            dict: Dictionary containing number of barcodes (values) per taxonomy entry (keys),
                organized by primer pair.
        """
        # print("mozaiko INFO: Calculating the number of barcodes per each taxon per FASTA file.")

        # Get all FASTA files
        fasta_files = self.parse_fasta_files()
        barcodes_per_entry = {}

        # Process each primer pair and its corresponding file
        for primer_pair, file_path in fasta_files.items():
            # Read FASTA file with minimal validation for speed
            fasta_data = self.custom_fasta_import.read_fasta(
                file_path,
                check_taxid=False,
                overide_validation=True,
                taxa_column_end=11,
            )

            # Convert taxa_info to string type and handle NaN values
            fasta_data["taxa_info"] = (
                fasta_data["taxa_info"].astype(str).replace("nan", "")
            )

            # Remove any leading/trailing whitespace
            fasta_data["taxa_info"] = fasta_data["taxa_info"].str.strip()

            # Filter out empty taxa_info entries
            fasta_data = fasta_data[fasta_data["taxa_info"] != ""]

            # Use groupby for counting
            counts = fasta_data.groupby("taxa_info")["sequence"].size()

            # Convert Series to dictionary
            barcodes_per_entry[primer_pair] = counts.to_dict()

        return barcodes_per_entry

    def calculate_number_of_barcodes_per_otl_taxonomy(
        self, barcodes_per_entry: dict, otl_hierarchical_taxonomy: dict
    ):
        """
        Calculate barcode counts per taxonomic level, handling missing species/genus data.
        """
        # Replicate OTL taxonomic hierarchy for each primer
        # Allows us to add barcode counts based on OTL structure
        primer_taxa_data = {}
        for primer in barcodes_per_entry.keys():
            primer_taxa_data[primer] = copy.deepcopy(otl_hierarchical_taxonomy)
        # print(primer_taxa_data)

        def normalize_taxon(value):
            if value in ["nan"]:
                return None
            return value

        # Create mapping between parsed taxonomy in headers and their occurences (counts)
        taxa_counts_mapping = {}
        for primer, unique_header_counts in barcodes_per_entry.items():
            for taxa_header, count in unique_header_counts.items():
                taxa_parts = taxa_header.split("|")
                family = normalize_taxon(taxa_parts[7])
                genus = normalize_taxon(taxa_parts[8] if len(taxa_parts) > 8 else None)
                species = normalize_taxon(taxa_parts[9] if len(taxa_parts) > 9 else None)

                mapping_key = (primer, family, genus, species)
                taxa_counts_mapping[mapping_key] = count

        def update_counts(primer_otl_hierarchy, taxa_counts_mapping, primer):
            """
            Updates counts for each taxonomic level in the OTL-based hierarchy according to the
            number of barcodes.
            """
            # print("mozaiko INFO: Updating counts for each taxonomic level.")
            for family, family_data in primer_otl_hierarchy.items():
                family_total = 0

                # count entries that only have family level information
                family_only_count = sum(
                    count
                    for (p, f, g, s), count in taxa_counts_mapping.items()
                    if p == primer
                    and f == family
                    and (g == "nan" or g is None)
                    and (s == "nan" or s is None)
                )

                if "genera" in family_data:
                    for genus, genus_data in family_data["genera"].items():
                        genus_total = 0

                        # count entries that only have genus level information
                        genus_only_count = sum(
                            count
                            for (p, f, g, s), count in taxa_counts_mapping.items()
                            if p == primer
                            and f == family
                            and g == genus
                            and (s == "nan" or s is None)
                        )

                        if "species" in genus_data:
                            for species_name, species_data in genus_data[
                                "species"
                            ].items():
                                # count entries with species level information
                                # excludes 'nan' values from count
                                species_total = sum(
                                    count
                                    for (
                                        p,
                                        f,
                                        g,
                                        s,
                                    ), count in taxa_counts_mapping.items()
                                    if p == primer
                                    and f == family
                                    and g == genus
                                    and s
                                    and s != "nan"
                                    and s == species_name
                                )

                                # add species-level counts to genus count
                                species_data["count"] = species_total
                                genus_total += species_total

                        # add genus-only counts to genus total
                        genus_total += genus_only_count
                        genus_data["count"] = genus_total

                        # add all genus counts (including species) to family total
                        family_total += genus_total

                # Add family-only counts to family total
                family_total += family_only_count
                family_data["count"] = family_total

        # update counts for each primer according to barcode numbers
        for primer in barcodes_per_entry.keys():
            update_counts(primer_taxa_data[primer], taxa_counts_mapping, primer)
            # print(f"mozaiko INFO: Counts updated for {primer}.")

        return primer_taxa_data

    def calculate_proportion_of_taxa_w_x_barcodes(
        self,
        total_taxa_count: int,
        barcode_threshold: int = 1,
    ):
        """
        Calculates the proportion of taxa with more than X barcodes.
        """
        barcodes_per_entry = self.calculate_number_of_barcodes_per_fasta_entry()
        barcodes_per_species = self.calculate_number_of_barcodes_per_otl_taxonomy(
            barcodes_per_entry, self.taxa_hierarchy
        )

        def count_qualifying_taxa(taxa_data, threshold):
            # print(f"mozaiko INFO: Counting taxa with more than {barcode_threshold} barcodes.")
            qualifying_taxa = 0

            for family, family_data in taxa_data.items():
                if not isinstance(family_data, dict):
                    continue

                family_has_qualifying_descendant = False

                if "genera" in family_data:
                    for genus, genus_data in family_data["genera"].items():
                        if not isinstance(genus_data, dict):
                            continue

                        genus_has_qualifying_descendant = False

                        if "species" in genus_data:
                            for species_name, species_data in genus_data[
                                "species"
                            ].items():
                                if (
                                    isinstance(species_data, dict)
                                    and species_data.get("count", 0) >= threshold
                                ):
                                    qualifying_taxa += 1
                                    genus_has_qualifying_descendant = True
                                    family_has_qualifying_descendant = True

                        if (
                            not genus_has_qualifying_descendant
                            and isinstance(genus_data, dict)
                            and genus_data.get("count", 0) >= threshold
                        ):
                            qualifying_taxa += 1
                            family_has_qualifying_descendant = True

                if (
                    not family_has_qualifying_descendant
                    and isinstance(family_data, dict)
                    and family_data.get("count", 0) >= threshold
                ):
                    qualifying_taxa += 1

            return qualifying_taxa

        results = {}
        for primer_pair, taxa_data in barcodes_per_species.items():
            taxa_meeting_threshold = count_qualifying_taxa(taxa_data, barcode_threshold)
            proportion = round(taxa_meeting_threshold / total_taxa_count, 2)
            results[primer_pair] = proportion
            # print(f"mozaiko INFO: Calculated proportion for {primer_pair}.")

        return results

    def barcoded_taxa_ratio(self, total_taxa_count: int):
        """
        This method calculates the Ratio of Barcoded Taxa (RBT).

        Parameters:
        - barcoded_taxa_five_plus: float
            Percentage of taxa with more than five barcodes.
        - barcoded_taxa: float
            Percentage of taxa with more than one barcode.

        Output:
        - barcoded_taxa_ratio: Dict
            A dictionary containing the ratio of barcoded taxa and the proportion of taxa with more than five barcodes per primer pair.
        """
        barcoded_taxa_five_plus = self.calculate_proportion_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=5
        )
        barcoded_taxa = self.calculate_proportion_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=1
        )

        barcoded_taxa_ratio = {}

        for primer_pair in barcoded_taxa_five_plus.keys():
            percent_5plus = barcoded_taxa_five_plus[primer_pair]
            percent_1plus = barcoded_taxa[primer_pair]

            ratio_barcoded_taxa = (
                percent_5plus / percent_1plus if percent_5plus > 0 else 0
            )

            barcoded_taxa_ratio[primer_pair] = {
                "barcoded_taxa": percent_1plus,
                "ratio_barcoded_taxa": round(ratio_barcoded_taxa, 2),
            }

        barcoded_taxa_ratio_df = pd.DataFrame(barcoded_taxa_ratio)

        return barcoded_taxa_ratio_df.T


class Binding:
    """
    This class implements the method related to evaluating the binding efficency and performance
    between the primer-sets and the PBS.
    """

    def __init__(self, otl, number_of_mismatches=None):
        self.amplification_instance = InSilicoAmplification()
        if number_of_mismatches is None:
            number_of_mismatches = (
                self.amplification_instance.get_number_of_mismatches()
            )
        self.number_of_mismatches = number_of_mismatches
        self.processed_primers = {}
        self.otl_handler = OtlHandler(otl)
        self.otl_handler.import_otl()
        self.otl = self.otl_handler.otl

    @staticmethod
    def parse_files_with_same_extension_in_folders(folder_path_A, folder_path_B):
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


    @staticmethod
    def get_pbs_table(amplicon_file, insert_file):
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

        The method processes amplicon and insert files from the in-silico amplification step,
        calculating metrics such as full-length mismatches, 3′-end mismatches, GC clamp matches,
        minimum melting temperature, and primer GC fractions for each primer pair.

        Parameters:
        amplicon_folder : str
            Path to the folder containing amplicon sequence files.
        insert_folder : str
            Path to the folder containing insert sequence files.
        primer_table : str
            Path to the table containing primer pair sequences and metadata.
        save_results : bool, optional
            If True, saves output DataFrames as CSV files. Default is False.

        Output:
        Tuple[Dict[str, pd.DataFrame], pd.DataFrame]
        - primer_pbs_df : dict
            Dictionary of DataFrames indexed by primer pair identifier, containing:
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
        print("mozaiko INFO: Retriving information on the Binding Efficiency.")

        self.primer_table = self.get_primer_table(primer_table)
        matching_files = Binding.parse_files_with_same_extension_in_folders(
            amplicon_folder, insert_folder
        )

        if not matching_files:
            print("mozaiko ERROR: No matching primer files found between the insert and amplicon folders.")
            return None, None

        primer_pbs_df = {}
        primer_gc_fractions = []

        for _primer_ind, primer_row in self.primer_table.iterrows():
            barcode_region = primer_row["barcode_region"]
            assay_name = primer_row["assay_name"]
            primer_name = barcode_region + "_" + assay_name
            pbs_filename = f"{barcode_region}_{assay_name}"
            primer_seq_fwd = primer_row["fwd_seq"]
            primer_seq_rev = primer_row["rev_seq"]
            rev_comp_primer_seq_rev = str(Seq(primer_seq_rev).reverse_complement())

            # Compute GC fraction for Primer
            primer_gc_fractions.append(
                {
                    "primer_name": primer_name,
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
                    pbs_table = Binding.get_pbs_table(amplicon_file, insert_file)

                    for _, pbs_row in pbs_table.iterrows():
                        pbs_fwd_seq = pbs_row["fwd_seq"]
                        pbs_rev_seq = pbs_row["rev_seq"]

                        if not pbs_fwd_seq or not pbs_rev_seq:
                            seq_header = pbs_row["header"].replace(">", "")
                            seq_id = seq_header.split("|")[0].replace(" ", "")
                            print(
                                f"mozaiko WARNING: Skipping entry with empty PBS sequence(s): {seq_id}"
                            )
                            continue

                        seq_header = pbs_row["header"].replace(">", "")
                        seq_id = seq_header.split("|")[0].replace(" ", "")
                        taxon = seq_header.split("|")[2]
                        rank = seq_header.split("|")[3]
                        family = seq_header.split("|")[8]
                        genus = seq_header.split("|")[9]
                        species = seq_header.split("|")[10]

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
                        # delta_tm = abs(substraction)
                        # delta_tm = round(delta_tm, 2)

                        pbs_melting_temperature_data.append(
                            {
                                "seq_id": seq_id,
                                "taxon": taxon,
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
                                "min_tm": min_tm
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
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
                                "full_len_mismatch_sum": full_fwd_mismatches
                                + full_rev_mismatches,
                            }
                        )

                        # Compute three-end mismatches and GC clamp
                        three_end_fwd_seq = pbs_fwd_seq[-5:]  # last 5 bases
                        three_end_rev_seq = pbs_rev_seq[:5]  # first 5 bases
                        three_end_fwd_primers = primer_seq_fwd[-5:]  # last 5 bases
                        three_end_rev_primers = rev_comp_primer_seq_rev[
                            :5
                        ]  # first 5 bases (reverse complement)

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
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
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
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
                                "gc_matches_fwd": three_end_fwd_gc_matches,
                                "gc_matches_rev": three_end_rev_gc_matches,
                            }
                        )

            if matching_files_found:
                if (
                    not full_mismatches_data
                    or not three_end_mismatches_data
                    or not pbs_melting_temperature_data
                    or not three_end_gc_matches_data
                ):
                    print(
                        f"mozaiko WARNING: Missing data for {pbs_filename}, skipping..."
                    )
                    continue

                try:
                    # Start by adding one df
                    comprehensive_df = pd.DataFrame(full_mismatches_data)

                    if comprehensive_df.empty:
                        print(
                            f"mozaiko WARNING: Empty data for {pbs_filename}, skipping..."
                        )
                        continue

                    # Create and verify other DataFrames before merging
                    three_end_df = pd.DataFrame(three_end_mismatches_data)
                    temp_df = pd.DataFrame(pbs_melting_temperature_data)
                    gc_matches_df = pd.DataFrame(three_end_gc_matches_data)

                    required_columns = {
                        "three_end_df": ["seq_id", "taxon", "three_end_mismatch_sum"],
                        "temp_df": ["seq_id", "taxon", "min_tm"],
                        "gc_matches_df": [
                            "seq_id",
                            "taxon",
                            "gc_matches_fwd",
                            "gc_matches_rev",
                        ],
                    }

                    for df_name, columns in required_columns.items():
                        df = locals()[df_name]
                        missing_cols = [col for col in columns if col not in df.columns]
                        if missing_cols:
                            print(
                                f"mozaiko WARNING: Missing columns {missing_cols} in {df_name} for {pbs_filename}"
                            )
                            continue

                    #  Merge all other df's by seq_id & taxon
                    comprehensive_df = comprehensive_df.merge(
                        pd.DataFrame(three_end_mismatches_data)[
                            ["seq_id", "taxon", "three_end_mismatch_sum"]
                        ],
                        on=["seq_id", "taxon"],
                    )

                    comprehensive_df = comprehensive_df.merge(
                        pd.DataFrame(pbs_melting_temperature_data)[
                            ["seq_id", "taxon", "min_tm"]
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
                        "family",
                        "genus",
                        "species",
                        "rank",
                        "full_len_mismatch_sum",
                        "three_end_mismatch_sum",
                        "gc_matches_fwd",
                        "gc_matches_rev",
                        "min_tm"
                    ]

                    comprehensive_df = comprehensive_df[column_order]

                    primer_pbs_df[pbs_filename] = comprehensive_df

                    if save_results:
                        comprehensive_df.to_csv(
                            f"{pbs_filename}_comprehensive.csv", index=False
                        )

                except Exception as e:
                    print(f"mozaiko ERROR: Failed to process {pbs_filename}: {str(e)}")
                    continue

        primer_gc_df = pd.DataFrame(primer_gc_fractions)
        primer_gc_df.set_index("primer_name", inplace=True)

        if save_results:
            primer_gc_df.to_csv("primer_gc_fractions.csv", index=False)

        return primer_pbs_df, primer_gc_df

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

    def fill_hierarchical_values(self, df, operation_name, analysis_name=None):
        """
        Fill missing values in a hierarchical dataset based on genus-level and family-level
        aggregations.

        Parameters:
        - df (pd.DataFrame): The input DataFrame with columns ['family', 'genus', 'species'] and
            numerical data.
        - operation_name (str): The operation to perform ('min', 'max', 'sum', 'mean', or
        'coef_var').
        - analysis_name (str or None): The specific column to analyze. If None, applies to all
        numeric columns.

        Returns:
            pd.DataFrame: A DataFrame with missing values filled hierarchically.
        """
        # Get OTL reference with rank information
        otl = self.otl_handler.otl[["family", "genus", "species", "rank"]].drop_duplicates().copy()

        # If analysis_name is None, apply the logic to all numeric columns
        if analysis_name is None:
            numeric_cols = df.select_dtypes(include="number").columns
        else:
            numeric_cols = [analysis_name]

        # Initialize result with OTL taxonomy and empty numeric columns
        result = otl.copy()
        for col in numeric_cols:
            result[col] = np.nan

        # For each OTL entry, aggregate values from df based on rank
        for idx, row in otl.iterrows():
            rank = row['rank']

            # Build filter based on rank - filter df by the taxonomic name at this rank
            if rank == 'species':
                subset = df[
                    (df['family'] == row['family']) &
                    (df['genus'] == row['genus']) &
                    (df['species'] == row['species'])
                ]
            elif rank == 'genus':
                subset = df[
                    (df['family'] == row['family']) &
                    (df['genus'] == row['genus'])
                ]
            elif rank == 'family':
                subset = df[df['family'] == row['family']]
            else:
                subset = pd.DataFrame()  # Empty if rank is unknown

            # Aggregate each numeric column
            if not subset.empty:
                for col in numeric_cols:
                    values = subset[col].dropna()  # Remove NaN values before aggregation

                    if not values.empty:
                        agg = self._compute_aggregation(values, operation_name)
                    else:
                        agg = np.nan

                    result.loc[idx, col] = agg

        # Drop rank column before returning
        result.drop(columns=["rank"], inplace=True)

        return result

    def _compute_aggregation(self, series, operation_name):
        """
        Helper method to compute the aggregation based on operation_name.

        Parameters:
        - series (pd.Series): The data to aggregate
        - operation_name (str): The operation to perform

        Returns:
            Aggregated value
        """
        if operation_name == "min":
            return series.min()
        elif operation_name == "max":
            return series.max()
        elif operation_name == "sum":
            return series.sum()
        elif operation_name == "mean":
            return series.mean()
        elif operation_name == "coef_var":
            mean = series.mean()
            std = series.std()
            if mean == 0:
                return np.nan
            return (std / mean) * 100
        else:
            raise ValueError(f"Unknown operation: {operation_name}")

    def process_analysis_per_taxon(
        self,
        primer_df: pd.DataFrame,
        operation: Literal["min", "max", "sum", "mean", "coef_var"],
        analysis_name: str,
        ) -> Union[pd.DataFrame, pd.Series]:
        """
        This method performs user-inputed operations on a groupby of 'taxon'.
        Modified to handle empty groups consistently across operations.

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

        grouped_taxa = primer_df.groupby(["family", "genus", "species"])

        if operation == "sum":
            result = pd.DataFrame({analysis_name: grouped_taxa[analysis_name].sum()})
        elif operation in {"min", "max", "mean"}:
            result = pd.DataFrame(
                {analysis_name: getattr(grouped_taxa[analysis_name], operation)()}
            )
        elif operation == "coef_var":
            result = primer_df[
                ["family", "genus", "species", analysis_name]
            ].copy()
        else:
            raise ValueError(
                f"mozaiko ERROR: Unrecognized operation: '{operation}'. "
                "Please choose from 'min', 'max', 'sum', 'mean', or 'coef_var'."
            )

        # Reset index to get family, genus, species as columns
        if "family" not in result.columns:
            result.reset_index(inplace=True)
        result = result.replace("nan", np.nan)
        hierarchical_result = self.fill_hierarchical_values(result, operation)

        return hierarchical_result

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

        numeric_df = tax_grouped_df.select_dtypes(include=["number"])

        if operation in {"min", "max", "sum", "mean"}:
            result = getattr(numeric_df, operation)().astype(float)
        elif operation == "coef_var":
            mean = numeric_df.mean()
            std = numeric_df.std()
            result = (std / mean) * 100
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
        - ratio_df: DataFrame
            A DataFrame containing the priming ratio for each taxon, calculated as the sum of
        """
        merged_df = pd.merge(
            max_mismatch_three_end,
            max_mismatch_full_len,
            on=["family", "genus", "species"],
            suffixes=("_three_end", "_full_len"),
        )
        merged_df["priming_ratio"] = (
            merged_df["three_end_mismatch_sum"] / merged_df["full_len_mismatch_sum"]
        )

        ratio_col = merged_df["priming_ratio"]
        ratio_df = pd.DataFrame(ratio_col)

        return ratio_df

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


class TraitsAndResolution:
    def __init__(
        self,
        otl,
        results_folder: Optional[str] = None,
        insert_folder_path: Optional[str] = None,
        amplicon_folder_path: Optional[str] = None,
        incomplete_pbs_folder_path: Optional[str] = None,
    ):
        # print(
        #     f"mozaiko INFO: Provided output folder: {results_folder}."
        # )
        if results_folder is not None:
            self.results_folder = results_folder
            self.insert_folder_path = os.path.join(results_folder, "insert/filtered")
            self.amplicon_folder_path = os.path.join(results_folder, "amplicon")
            self.incomplete_pbs_path = os.path.join(
                results_folder, "incomplete_pbs/filtered"
            )
            # print(
            #     f"mozaiko INFO: Setting insert folder path to {self.insert_folder_path} and amplicon folder path to {self.amplicon_folder_path} and incomplete_pbs_path to {self.incomplete_pbs_path}."
            # )
        elif (
            insert_folder_path is not None
            and amplicon_folder_path is not None
            and incomplete_pbs_folder_path is not None
        ):
            self.insert_folder_path = insert_folder_path
            self.amplicon_folder_path = amplicon_folder_path
            self.incomplete_pbs_path = incomplete_pbs_folder_path
            self.results_folder = os.path.dirname(os.path.dirname(insert_folder_path))
            # print(
            #     "mozaiko INFO: Using provided insert_folder_path, amplicon_folder_path and incomplete_pbs_folder_path."
            # )
        else:
            raise ValueError(
                "mozaiko ERROR: Either provide a path to the in-silico amplification results folder "
                "('results_folder') or paths to the insert ('insert_folder_path'), amplicon "
                "('amplicon_folder_path') and incomplete PBS ('incomplete_pbs_path') results folders."
            )
        self.otl_handler = OtlHandler(otl)
        self.otl_handler.import_otl()
        self.otl = self.otl_handler.otl
        self.otl = self.otl[['family', 'genus', 'species', 'rank']]
        self.otl = self.otl.replace(np.nan, "nan")
        self.otl = self.otl.drop_duplicates()
        self.otl = self.otl.reset_index(drop=True)
        self.country_name = Path(otl).stem.split('_')[0]
        self.catnip_dir = os.path.join(self.results_folder, "catnip")

    def run_catnip(self, threshold: float | list = 10.0):
        clustering_threshold = None
        if isinstance(threshold, float):
            clustering_threshold = threshold
        elif isinstance(threshold, list):
            if len(threshold) != 3:
                raise TypeError("mozaiko ERROR: threhsolds must be a list of three float numbers to be used in family, genus and species filtering, respectively.")
            clustering_threshold = max(threshold)
        else:
            raise TypeError("mozaiko ERROR: threhsolds must be float or a list of three numbers.")

        # create directory for catnip analysis
        os.makedirs(self.catnip_dir, exist_ok=True)

        print("mozaiko INFO: Starting catnip to retrieve nucleotide divergence across taxa levels...")

        # get catnip script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        catnip_script = os.path.join(script_dir, "catnip.sh")
        os.chmod(catnip_script, 0o755)

        # run catnip for each primer
        processed_files = 0
        failed_files = []

        for file in os.listdir(self.incomplete_pbs_path):
            if file.endswith(".fasta"):
                original_fasta_path = os.path.join(self.incomplete_pbs_path, file)
                primer_name = os.path.splitext(file)[0]

                # Create a subdirectory for each primer to keep outputs organized
                primer_output_dir = os.path.join(self.catnip_dir, primer_name)

                os.makedirs(primer_output_dir, exist_ok=True)

                # Copy FASTA to the primer's output directory
                catnip_fasta_path = os.path.join(primer_output_dir, file)
                shutil.copy2(original_fasta_path, catnip_fasta_path)

                mapping_file_name = f"mapping_{primer_name}.tsv"
                cols = "0,8,9,10"  # seq_id, family, genus, species

                try:
                    result = subprocess.run(
                        [
                            catnip_script,
                            file,
                            mapping_file_name,
                            cols,
                            str(clustering_threshold)
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                        cwd=primer_output_dir  # Run in primer-specific directory
                    )

                    processed_files += 1
                    print(f"mozaiko INFO: catnip completed successfully for {primer_name}.")

                except subprocess.CalledProcessError as e:
                    failed_files.append(primer_name)
                    print(f"mozaiko ERROR: catnip failed for {primer_name}.")
                    print(f"mozaiko ERROR: Error code: {e.returncode}")
                    if e.stderr:
                        print(f"mozaiko ERROR: {e.stderr.strip()}")
                    if e.stdout:
                        print(f"mozaiko ERROR: {e.stdout.strip()}")

                except Exception as e:
                    failed_files.append(primer_name)
                    print(f"mozaiko ERROR: Unexpected error processing {primer_name}: {str(e)}")

        print(f"\nmozaiko INFO: catnip processing complete.")
        # print(f"mozaiko INFO: Output directory: {self.catnip_dir}")
        # print(f"mozaiko INFO: Successfully processed: {processed_files} primers(s)")

        if failed_files:
            print(f"mozaiko WARNING: Failed to process {len(failed_files)} file(s): {', '.join(failed_files)}")

    def clean_and_split(self, df, col_prefix):
        """
        Cleans up taxa columns and splits them into family, genus, species.
        Keeps NaN as np.nan.

        Args:
            df: DataFrame to process
            col_prefix: Either "query" or "target"

        Returns:
            DataFrame with added family, genus, species columns
        """
        # col = f"{col_prefix}_cat"

        # # Remove brackets and parentheses of tuple
        # df[col] = df[col].str.replace(r"[\'\(\)]", "", regex=True)

        # new_cols = [f"{col_prefix}_family", f"{col_prefix}_genus", f"{col_prefix}_species"]
        # df[new_cols] = df[col].str.split(',', expand=True)

        # df[new_cols] = df[new_cols].apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        # return df


        col = f"{col_prefix}_cat"

        new_cols = [
            f"{col_prefix}_family",
            f"{col_prefix}_genus",
            f"{col_prefix}_species",
        ]

        # Work on a safe string version but preserve NaNs
        s = df[col].astype("string")

        # Remove brackets and parentheses of tuple
        s = s.str.replace(r"[\'()]", "", regex=True)
        split = s.str.split(",", expand=True)
        split = split.reindex(columns=range(3))
        split.columns = new_cols

        # Strip whitespace
        split = split.map(
        lambda x: x.strip() if isinstance(x, str) else x
        )

        # Keep NaNs as NaNs
        split = split.where(df[col].notna())

        df = df.join(split)

        return df


    def make_dataframe_symmetric(self, df):
        """
        Creates symmetric version of dataframe by swapping query and target columns.
        This ensures that both query->target and target->query relationships are present.

        Args:
            df: DataFrame to make symmetric

        Returns:
            Concatenated DataFrame with original + swapped rows
        """
        df_copy = df.copy()

        rename_map = {
            'query': 'target',
            'target': 'query',
            'query_family': 'target_family',
            'target_family': 'query_family',
            'query_genus': 'target_genus',
            'target_genus': 'query_genus',
            'query_species': 'target_species',
            'target_species': 'query_species'
        }

        # Swap query and target columns by renaming
        df_copy = df_copy.rename(columns=rename_map)

        first_cols = [
            'query', 'query_family', 'query_genus', 'query_species',
            'target', 'target_family', 'target_genus', 'target_species'
        ]
        rest_cols = [col for col in df_copy.columns if col not in first_cols]
        df_copy = df_copy[first_cols + rest_cols]

        # Concatenate original and swapped dataframes
        symmetric_df = pd.concat([df, df_copy], ignore_index=True)

        return symmetric_df


    def add_taxa_from_mapping(self, df, folder_path, primer_name):
        """
        Add taxa present in mapping_*.tsv file but not in the dataframe.
        These are taxa that were in homogeneous clusters and didn't get aligned by Bowtie.

        For these taxa:
        - Add them as query taxa
        - Set target to 'nan'
        - Set divergence_prct to 'inf' (couldn't be calculated due to clustering threshold)

        Args:
            df: DataFrame to augment
            folder_path: Path to the primer folder
            primer_name: Name of the primer (for logging/debugging)

        Returns:
            DataFrame with additional rows for missing taxa
        """
        mapping_cols = [
            'seq_id', 'og_taxa', 'scientificName', 'rank',
            'kingdom', 'phylum', 'class', 'order',
            'family', 'genus', 'species'
        ]

        # Read mapping file for THIS primer only
        mapping_files = list(folder_path.glob("mapping_*.tsv"))

        if not mapping_files:
            print(f"mozaiko WARNING: No mapping file found for {primer_name}")
            return df

        mapping_file = pd.read_csv(mapping_files[0], sep="\t", names=mapping_cols, header=0, index_col=False)

        tax_cols = ['family', 'genus', 'species', 'seq_id']
        mapping_df = mapping_file[tax_cols].copy()
        mapping_df = mapping_df.fillna('nan')

        # Create a map from taxa tuple to seq_id
        seqid_map = {}
        for _, row in mapping_df.iterrows():
            taxa_tuple = (row['family'], row['genus'], row['species'])
            seqid_map[taxa_tuple] = row['seq_id']

        # Get all unique taxa from mapping file
        all_taxa_in_mapping = mapping_df[['family', 'genus', 'species']].drop_duplicates()

        # Get existing taxa from current dataframe (both query and target sides)
        existing_query = df[['query_family', 'query_genus', 'query_species']].rename(
            columns={'query_family': 'family', 'query_genus': 'genus', 'query_species': 'species'}
        )
        existing_target = df[['target_family', 'target_genus', 'target_species']].rename(
            columns={'target_family': 'family', 'target_genus': 'genus', 'target_species': 'species'}
        )

        # Combine all existing taxa
        existing_taxa = pd.concat([existing_query, existing_target], ignore_index=True).drop_duplicates()

        # Find taxa that's in mapping but not in dataframe
        merge_cols = ['family', 'genus', 'species']
        new_taxa = pd.merge(all_taxa_in_mapping, existing_taxa, on=merge_cols, how='left', indicator=True)
        new_taxa = new_taxa[new_taxa['_merge'] == 'left_only'].drop(columns=['_merge'])

        # Create rows for new taxa
        new_rows = []
        for _, row in new_taxa.iterrows():
            taxa_tuple = (row['family'], row['genus'], row['species'])
            seq_id = seqid_map.get(taxa_tuple, pd.NA)

            new_row = {
                'query': seq_id,
                'query_family': row['family'],
                'query_genus': row['genus'],
                'query_species': row['species'],
                'target': 'nan',
                'target_family': 'nan',
                'target_genus': 'nan',
                'target_species': 'nan',
                'divergence_prct': np.inf  # Couldn't be calculated (homogeneous cluster)
            }

            new_rows.append(new_row)

        # Add new taxa rows to dataframe
        if new_rows:
            new_taxa_df = pd.DataFrame(new_rows)
            df = pd.concat([df, new_taxa_df], ignore_index=True)

        return df

    def process_single_primer(self, folder_name, df, folder_path, thresholds: float | list = [10.0, 5.0, 2.0]):
        """
        Process a single primer's catnip output:
        1. Clean and split taxonomy columns
        2. Make dataframe symmetric (swap query<->target)
        3. Add taxa from mapping file that weren't in heterogeneous clusters (as inf values)

        Args:
            folder_name: Name of the primer folder
            df: DataFrame with catnip results
            folder_path: Path to the primer folder

        Returns:
            Processed DataFrame
        """
        print(f"mozaiko INFO: Processing primer {folder_name}...")

        clustering_threshold = None
        if isinstance(thresholds, float):
            clustering_threshold = thresholds
        elif isinstance(thresholds, list):
            if len(thresholds) != 3:
                raise TypeError("mozaiko ERROR: threhsolds must be a list of three float numbers to be used in family, genus and species filtering, respectively.")
            clustering_threshold = max(thresholds)
        else:
            raise TypeError("mozaiko ERROR: threhsolds must be float or a list of three numbers.")

        # Step 1: Clean and split taxonomy columns
        df = self.clean_and_split(df, "query")
        df = self.clean_and_split(df, "target")

        ordered_cols = [
            "query", "query_family", "query_genus", "query_species",
            "target", "target_family", "target_genus", "target_species",
            "divergence_prct"
        ]
        df = df[ordered_cols]

        # Fill NaN with 'nan' string for taxonomy match
        df = df.fillna('nan')

        # Step 2: Replace entries that are above the clustering threshold with 'inf'
        df_above_clustering_th = df[df['divergence_prct'] > clustering_threshold].copy()
        df = df[df['divergence_prct'] <= clustering_threshold]
        df_above_clustering_th.loc[:, 'divergence_prct'] = np.inf
        df_above_clustering_th.loc[
            :, ['target_family', 'target_genus', 'target_species']
        ] = 'nan'

        # Step 3: Make symmetric (add swapped query/target rows)
        df_symm = self.make_dataframe_symmetric(df)

        df_symm = pd.concat([df_symm, df_above_clustering_th], ignore_index=True)

        # Step 4: Add taxa from mapping file that were not processed with bowtie (infs)
        df_inf = self.add_taxa_from_mapping(df_symm, folder_path, folder_name)
        df_inf = df_inf.drop(['query', 'target'], axis=1)
        self.df_inf = df_inf.copy()
        processed_file_path = folder_path / f"df_inf_{folder_name}_{self.country_name}.tsv"
        df_inf.to_csv(processed_file_path, sep='\t', index=False, header=True)

        # # Step 5: Find minimum divergence_prct for each OTL entry in query and all target taxa
        df_catnipt_all_on_target = self.get_values_otl(df_inf)
        processed_file_path = folder_path / f"catnip_target_prefilt_{folder_name}_{self.country_name}.tsv"
        df_catnipt_all_on_target.to_csv(processed_file_path, sep='\t', index=False, header=True)

        # Step 6: Filter target taxa by OTL & find minimum divergence_prct
        df_otl_on_target = self.filter_results_by_otl(df_inf)
        df_otl_on_target = self.get_values_otl(df_otl_on_target)

        # # Step 7: Filter both dfs by thresholds
        df_otl_on_target = self.filter_divergence_threshold(df_otl_on_target, thresholds)
        processed_file_path = folder_path / f"otl_target_{folder_name}_{self.country_name}.tsv"
        df_otl_on_target.to_csv(processed_file_path, sep='\t', index=False, header=True)

        df_catnipt_all_on_target = self.filter_divergence_threshold(df_catnipt_all_on_target, thresholds)
        processed_file_path = folder_path / f"catnip_target_{folder_name}_{self.country_name}.tsv"
        df_catnipt_all_on_target.to_csv(processed_file_path, sep='\t', index=False, header=True)

        return df_otl_on_target, df_catnipt_all_on_target

    def post_process_catnip_primer_results(self, thresholds: float | list = [10.0, 5.0, 2.0]):
        """
        Process each primer's final_output_interclst.tsv file individually.
        Saves each processed file as processed_<primer_name>.tsv in the same folder.~

        Returns: None. Saves processed files to disk.
        """
        # print(f"mozaiko INFO: Processing {len(list(os.listdir(self.catnip_dir)))} primers...")

        for folder in os.listdir(self.catnip_dir):
            folder_path = Path(self.catnip_dir) / folder

            if not folder_path.is_dir():
                continue

            output_file_path = folder_path / "final_output_interclst.tsv"

            if not output_file_path.exists():
                print(f"mozaiko WARNING: No final_output_interclst.tsv found for {folder}")
                continue

            df = pd.read_csv(output_file_path, sep="\t")

            if 'edit_distance' in df.columns:
                df = df.drop('edit_distance', axis=1)

            self.process_single_primer(folder, df, folder_path, thresholds)

    def get_values_otl(self, df):
        """
        For each entry in OTL, find the minimum divergence_prct in df.
        Removes entries where query and target taxa are the same.
        If multiple entries have the same minimum divergence_prct, set target taxa to 'nan'.
        Args:
            df: DataFrame with query/target taxa and divergence_prct. Should be symmetric and
                contain all taxa from OTL as query taxa, with 'inf' for missing target taxa.
        Returns:
            DataFrame filtered by OTL with minimum divergence_prct per OTL entry
        """
        otl = self.otl.copy()

        otl = otl.rename(columns={'family':'query_family', 'genus': 'query_genus', 'species': 'query_species'})
        new_cols = ['target_family', 'target_genus', 'target_species', 'divergence_prct']

        for col in new_cols:
            if col != 'divergence_prct':
                otl[col] = 'nan'
            else:
                otl[col] = np.nan

        otl_filtered_df = otl.copy()
        for index, row in otl.iterrows():
            # Check highest rank and select taxonomy of that rank
            rank = row['rank']
            query_col = 'query_' + rank
            target_col = 'target_' + rank
            entry = row[query_col]

            # Filter df to find minimum divergence for this entry
            search_entry = df[df[query_col] == entry]
            search_entry = search_entry[search_entry[target_col] != entry]

            if search_entry.empty:
                if rank == 'species':
                    otl_filtered_df.loc[index, ['target_family', 'target_genus', 'target_species',
                                            'divergence_prct']] = ['nan', 'nan', 'nan', np.nan]
                    continue

                if rank == 'genus' or rank == 'family':
                    # Filter df_inf for this entry
                    search_entry_inf = self.df_inf[self.df_inf[query_col] == entry]
                    # print(search_entry_inf)

                    if not search_entry_inf.empty:
                        otl_filtered_df.loc[index, ['target_family', 'target_genus', 'target_species',
                                                'divergence_prct']] = ['nan', 'nan', 'nan', np.inf]
                    elif search_entry_inf.empty:
                        otl_filtered_df.loc[index, ['target_family', 'target_genus', 'target_species',
                                                'divergence_prct']] = ['nan', 'nan', 'nan', np.nan]

                    continue

            # Extract values as numpy array
            vals = search_entry['divergence_prct'].to_numpy()

            # Find minimum values
            finite_vals = vals[np.isfinite(vals)]      # only finite numbers (no inf, no NaN)

            if len(finite_vals) > 0:
                # Case 1: finite numbers exist, just use the smallest
                min_value = finite_vals.min()
            else:
                # Case 2: only inf or NaN exist
                if np.isinf(vals).any():
                    min_value = np.inf                 # at least one infinite, use it
                else:
                    min_value = np.nan                 # all NaN

            # check if there are multiple equal mins
            same_min_div = search_entry[search_entry['divergence_prct'] == min_value]

            if len(same_min_div) > 1:
                otl_filtered_df.loc[index, ['target_family', 'target_genus', 'target_species',
                                            'divergence_prct']] = [
                    'nan', 'nan', 'nan', min_value
                ]
            else:
                best_value = same_min_div.iloc[0]
                otl_filtered_df.loc[index,
                    ['target_family', 'target_genus', 'target_species', 'divergence_prct']
                ] = [
                    best_value['target_family'],
                    best_value['target_genus'],
                    best_value['target_species'],
                    best_value['divergence_prct']
                ]

        otl_filtered_df.drop(columns=['rank'], inplace=True)

        otl_filtered_df = otl_filtered_df.replace(
        {'nan': np.nan}
        )

        return otl_filtered_df


    def filter_results_by_otl(self, df):
        """
        Filter target taxa by taxa that's only in OTL.
        """
        df_filtered = self.otl.merge(
            df,
            how='inner',
            left_on=['family', 'genus', 'species'],
            right_on=['target_family', 'target_genus', 'target_species']
            )

        df_inf = df[df['divergence_prct'] == np.inf]

        # Combine both dataframes
        df_filtered = pd.concat([df_filtered, df_inf], ignore_index=True).drop_duplicates()

        df_filtered = df_filtered.replace(
        {'nan': np.nan, 'inf': np.inf}
        )

        return df_filtered

    def exclude_common_ancestry(self, otl_filtered_df):
        only_family = (
            otl_filtered_df['query_family'].notna() &
            otl_filtered_df['query_genus'].isna() &
            otl_filtered_df['query_species'].isna() &
            (otl_filtered_df['query_family'] != otl_filtered_df['target_family'])
        )

        family_genus = (
            otl_filtered_df['query_genus'].notna() &
            otl_filtered_df['query_species'].isna() &
            (
                (otl_filtered_df['query_family'] != otl_filtered_df['target_family']) |
                (otl_filtered_df['query_genus'] != otl_filtered_df['target_genus'])
            )
        )

        family_genus_species = (
            otl_filtered_df['query_species'].notna() &
            (
                otl_filtered_df['query_species'] != otl_filtered_df['target_species']
            )
        )

        different_values = only_family | family_genus | family_genus_species

        mask_keep = (
            different_values |
            otl_filtered_df['divergence_prct'].isin([np.inf, -np.inf]) |
            otl_filtered_df['divergence_prct'].isna()
        )

        otl_filtered_df = otl_filtered_df[mask_keep].copy()

        return otl_filtered_df

    def filter_divergence_threshold(self, df, thresholds: list | float = None):
        if thresholds is None:
            thresholds = [10.0, 5.0, 2.0]
        single_threhold = None
        if isinstance(thresholds, float):
            single_threhold = thresholds
            print(f"mozaiko INFO: Using single divergence threshold of {single_threhold}%.")
        elif isinstance(thresholds, list):
            th_family, th_genus, th_species = thresholds
            print(f"mozaiko INFO: Using divergence thresholds of {th_family}% for family, {th_genus}% for genus and {th_species}% for species.")
        else:
            raise TypeError("mozaiko ERROR: threhsolds must be float or a list of three numbers.")

        otl_filtered = df.copy()

        # Step 1: exclude rows with common ancestry
        otl_filtered = self.exclude_common_ancestry(otl_filtered)

        # Step 2: Apply rank-dependent divergence thresholds
        if single_threhold is not None:
            threshold_mask = (
                (otl_filtered['divergence_prct'] > single_threhold) |
                (otl_filtered['divergence_prct'].isin([np.inf, -np.inf])) |
                (otl_filtered['divergence_prct'].isna())
            )
        else:
            cond_family = (
                otl_filtered['query_family'].notna() &
                otl_filtered['query_genus'].isna() &
                otl_filtered['query_species'].isna() &
                (otl_filtered['divergence_prct'] > th_family)
            )

            cond_genus = (
                otl_filtered['query_genus'].notna() &
                otl_filtered['query_species'].isna() &
                (otl_filtered['divergence_prct'] > th_genus)
            )

            cond_species = (
                otl_filtered['query_species'].notna() &
                (otl_filtered['divergence_prct'] > th_species)
            )

            threshold_mask = (
                cond_family |
                cond_genus |
                cond_species |
                otl_filtered['divergence_prct'].isin([np.inf, -np.inf]) |
                otl_filtered['divergence_prct'].isna()
            )

        otl_filtered = otl_filtered[threshold_mask].copy()

        otl_filtered = (
            otl_filtered
            .sort_values('divergence_prct', ascending=True)
            .groupby(['query_family', 'query_genus', 'query_species'], dropna=False)
            .head(1)
            .reset_index(drop=True)
        )

        return otl_filtered

    def get_taxonomic_resolution(self):
        """
        This method retrieves the divergence score for each primer set. The divergence score is
        the percentage of taxa with a percentage of divergence bellow a cutoff according to the
        total number of taxa in the OTL.

        Parameters:
        - otl_total_taxa_cont: int
            The total number of taxa in the OTL.

        Returns:
        """
        total_otl_taxa_count = self.otl_handler.total_taxa

        taxonomic_resolution_results = []

        for folder in os.listdir(self.catnip_dir):
            folder_path = Path(self.catnip_dir) / folder

            if not folder_path.is_dir():
                continue

            catnip_target_file = folder_path / f"catnip_target_{folder}_{self.country_name}.tsv"
            otl_target_file = folder_path / f"otl_target_{folder}_{self.country_name}.tsv"

            if not catnip_target_file.exists() or not otl_target_file.exists():
                print(
                    f"mozaiko WARNING: Missing target files for primer '{folder}'. "
                    f"Expected files: {catnip_target_file.name}, {otl_target_file.name}"
                )
                continue

            try:
                df_catnipt_all_on_target = pd.read_csv(catnip_target_file, sep="\t")
                df_otl_on_target = pd.read_csv(otl_target_file, sep="\t")

                # Find all non -NaN divergence_prct entries in both dataframes
                df_catnipt_all_on_target = df_catnipt_all_on_target[
                    df_catnipt_all_on_target['divergence_prct'].notna()
                ]

                df_otl_on_target = df_otl_on_target[
                    df_otl_on_target['divergence_prct'].notna()
                ]

                # Metric 1: Taxonomic resolution considering only OTL entries

                taxa_considered_otl = df_otl_on_target.shape[0]

                taxonomic_resolution_otl = taxa_considered_otl / total_otl_taxa_count

                # Metric 2: Ratio between Taxonomic resolution considering only OTL entries and Taxonomic resolution considering all taxa from catnipt

                taxa_considered_all_catnip = df_catnipt_all_on_target.shape[0]

                taxonomic_resolution_all_catnip = taxa_considered_all_catnip / total_otl_taxa_count

                ratio_taxonomic_resolution = (
                    taxonomic_resolution_all_catnip / taxonomic_resolution_otl
                )

                taxonomic_resolution_results.append(
                    {
                        "primer": folder,
                        "total_taxa": total_otl_taxa_count,
                        "n_taxa_above_cutoff": taxa_considered_otl,
                        "taxonomic_resolution": round(taxonomic_resolution_otl, 2),
                        "ratio_taxonomic_resolution": round(ratio_taxonomic_resolution, 2),
                    }
                )

            except Exception as e:
                print(f"mozaiko ERROR: Failed to process primer '{folder}': {e}")
                continue

        taxonomic_resolution_df = pd.DataFrame(taxonomic_resolution_results)

        return taxonomic_resolution_df

class MetricsSystemExecutor:
    """
    This class orchestrates the entire Metrics System, coordinating the execution of all evaluation
    categories by calling the necessary classes and methods.
    """

    def __init__(self, results_folder: str, otl: str, primer_table: str):
        # Initialize Reference Database Quality Category
        # self.ref_db = ReferenceDatabaseQuality(otl=otl)
        # Initialize OTL and related variables
        self.otl = otl
        self.country_name = Path(otl).stem.split('_')[0]
        self.otl_handler = OtlHandler(self.otl)
        self.otl_handler.import_otl()
        self.total_otl_taxa_count = self.otl_handler.total_taxa
        # Load files and folders into memory
        self.primer_table = primer_table
        self.results_folder = results_folder
        self.insert_folder_path = os.path.join(results_folder, "insert/filtered")
        self.amplicon_folder_path = os.path.join(results_folder, "amplicon")
        self.incomplete_pbs_path = os.path.join(
            results_folder, "incomplete_pbs/filtered"
        )
        self.complete_pbs_path = os.path.join(
            results_folder, "all_complete_pbs/filtered"
        )
        print(
            f"mozaiko INDO: Setting insert folder path to {self.insert_folder_path} and amplicon folder path to {self.amplicon_folder_path} and incomplete_pbs_path to {self.incomplete_pbs_path}."
        )

    def get_reference_database_quality(self):
        """
        Initializes the Reference Database Quality evaluation, calculating the Barcode Coverage
        Score (BCS) for each primer and stores the results.

        Output:
        - ref_bd_scores: DataFrame
            A DataFrane containing the percentage of taxa with more than five barcodes and the
            rounded Ratio of Barcoded Taxa for each primer.
        """
        print("mozaiko INFO: Retriving information on the Reference Database Quality.")


        red_bd_qual = ReferenceDatabaseQuality(
            otl=self.otl, all_inserts_path=self.complete_pbs_path
        )
        reference_db_quality = red_bd_qual.barcoded_taxa_ratio(
            total_taxa_count=self.total_otl_taxa_count
        )

        return reference_db_quality

    def get_primer_pbs_analysis(self):
        """
        Initializes the primer-pbs analysis and
        """
        binding = Binding(otl=self.otl)
        primer_pbs_dict, gc_df = binding.primer_pbs_analysis(
            insert_folder=self.insert_folder_path,
            amplicon_folder=self.amplicon_folder_path,
            primer_table=self.primer_table,
        )

        gc_df.index.name = "primer"

        return primer_pbs_dict, gc_df

    def comprehensive_primer_analysis(self, output_folder, save_otl_level_results=True):
        """
        Perform comprehensive primer analysis across multiple metrics.

        Parameters:
        - cls : MetricsSystemExecutor
            An instance of MetricsSystemExecutor with initialized data
        - output_folder : str
            Path to the output folder for storing results

        Returns:
        - binding_df: pd.DataFrame
            Aggregated results of primer analyses
        """
        binding = Binding(otl=self.otl)

        ref_qual = self.get_reference_database_quality()

        primer_pbs, gc_df = self.get_primer_pbs_analysis()

        primer_results = {}

        primers_to_analyze = list(primer_pbs.keys())

        self.otl_level_filenames = []

        for primer in primers_to_analyze:
            primer_metrics = {}

            # Mismatch Analyses
            tax_lev_max_ms_full_len = binding.process_analysis_per_taxon(
                primer_pbs[primer],
                operation="max",
                analysis_name="full_len_mismatch_sum",
            )

            mismatch_score = int(
                binding.process_analysis_across_taxon(
                    tax_lev_max_ms_full_len, operation="sum"
                )
            )

            total_taxa_considered =len(tax_lev_max_ms_full_len[~tax_lev_max_ms_full_len['full_len_mismatch_sum'].isnull()])

            normalized_mismatch_score = mismatch_score / total_taxa_considered if total_taxa_considered > 0 else 0

            primer_metrics["normalized_mismatch_score"] = float(round(normalized_mismatch_score, 2))

            # Priming Ratio
            tax_lev_max_ms_three_end = binding.process_analysis_per_taxon(
                primer_pbs[primer],
                operation="max",
                analysis_name="three_end_mismatch_sum",
            )

            priming_ratio_df = binding.get_priming_ratio(
                tax_lev_max_ms_full_len, tax_lev_max_ms_three_end
            )

            total_taxa_considered = len(
                priming_ratio_df[~priming_ratio_df['priming_ratio'].isnull()]
            )

            priming_ratio_sum = binding.process_analysis_across_taxon(priming_ratio_df, operation="sum")

            priming_ratio_normalized = (
                priming_ratio_sum / total_taxa_considered
                if total_taxa_considered > 0
                else 0
            )

            priming_ratio_normalized = float(round(priming_ratio_normalized, 2))

            primer_metrics["normalized_priming_ratio_sum"] = priming_ratio_normalized

            # GC Match Analysis
            binding.get_total_gc_matches(primer_pbs[primer])

            tax_lev_gc = binding.process_analysis_per_taxon(
                primer_pbs[primer], operation="min", analysis_name="gc_matches_score"
            )

            total_taxa_considered = len(
                tax_lev_gc[~tax_lev_gc['gc_matches_score'].isnull()]
            )

            gc_matches_sum = (
                binding.process_analysis_across_taxon(tax_lev_gc, operation="sum")
            )

            gc_matches_sum_normalized = (
                gc_matches_sum / total_taxa_considered
                if total_taxa_considered > 0
                else 0
            )
            gc_matches_sum_normalized = float(round(gc_matches_sum_normalized, 2))

            primer_metrics["normalized_gc_matches_across_taxon"] = gc_matches_sum_normalized

            # Temperature Melting (Tm) Analysis
            tax_lev_min_tm = binding.process_analysis_per_taxon(
                primer_pbs[primer], operation="min", analysis_name="min_tm"
            )
            primer_metrics["min_tm_cv"] = binding.process_analysis_across_taxon(
                tax_lev_min_tm, operation="coef_var"
            )

            # # Tm Score
            # primer_metrics["tm_score"] = binding.tm_score(primer_pbs[primer])

            # # Amplification Success (if applicable)
            # try:
            #     amp_succ = binding.calculate_amplification_success_score(output_folder)

            #     amp_succ_value = (
            #         amp_succ.loc[primer, "amplification_sucess_percent"]
            #         if primer in amp_succ.index
            #         else None
            #     )
            #     primer_metrics["amplification_success_percent"] = amp_succ_value

            # except Exception as e:
            #     print(
            #         f"mozaiko WARNING: Amplification success calculation failed for {primer}: {e}"
            #     )
            #     primer_metrics["amplification_success_percent"] = None

            # Store results for this primer
            primer_results[primer] = primer_metrics

            # Save OTL-level results
            if save_otl_level_results:
                merge_columns = ["family", "genus", "species"]
                otl_base = self.otl_handler.otl[['family', 'genus', 'species']]  # ensure this is a dataframe

                otl_lev_result = otl_base.merge(tax_lev_max_ms_full_len, on=merge_columns, how="left") \
                                        .merge(tax_lev_max_ms_three_end, on=merge_columns, how="left") \
                                        .merge(tax_lev_gc, on=merge_columns, how="left") \
                                        .merge(tax_lev_min_tm, on=merge_columns, how="left")

                os.makedirs(
                    os.path.join(output_folder, "otl_level_results/"), exist_ok=True
                )
                output_path_otl_results = os.path.join(
                    output_folder, "otl_level_results", f"otl_level_results_{primer}.tsv"
                )
                self.otl_level_filenames.append(output_path_otl_results)
                otl_lev_result.to_csv(output_path_otl_results, sep="\t", index=False)

        binding_df = pd.DataFrame.from_dict(primer_results, orient="index")

        if ref_qual is not None:
            binding_df = ref_qual.join(binding_df)

        binding_df.index.name = "primer"

        binding_df = binding_df

        return binding_df

    def get_traits_and_resolution(self):
        """
        Combine taxonomic resolution and genetic divergence analyses

        Returns:
        pd.DataFrame
            Combined analysis results with taxonomic resolution and divergence metrics
        """
        print("mozaiko INFO: Retriving information on the Taxonomic Resolution.")
        trait = TraitsAndResolution(otl=self.otl, results_folder=self.results_folder)

        # Get Divergence Score
        taxonomic_resolution = trait.get_taxonomic_resolution()

        if "primer" in taxonomic_resolution.columns:
            taxonomic_resolution = taxonomic_resolution.set_index("primer")

        combined_div_score_and_tax_res_results = pd.DataFrame(
            {"taxonomic_resolution": taxonomic_resolution["taxonomic_resolution"],
             "ratio_taxonomic_resolution": taxonomic_resolution["ratio_taxonomic_resolution"]}
        )

        traits_res_df = combined_div_score_and_tax_res_results

        return traits_res_df

    def join_analysis_results(self):
        """
        This method joins the results from the primer analysis and the traits and resolution analysis.
        """
        binding_dataframe = self.comprehensive_primer_analysis(
            self.results_folder, save_otl_level_results=True
        )
        traits_dataframe = self.get_traits_and_resolution()

        analysis_results = binding_dataframe.join(traits_dataframe, on="primer")

        return analysis_results

    def create_taxonomic_resolution_per_taxon(self, div_score_df_path):
        """
        Method to create the taxonomic_resolution_per_taxon file that aggregates catnipt output
        files from each primer, comprising values for taxonomic resolution and divergence scores.
        """
        tax_res_df = []

        self.catnip_dir = Path(self.results_folder) / "catnip"

        for folder in os.listdir(self.catnip_dir):
            folder_path = Path(self.catnip_dir) / folder
            if not folder_path.is_dir():
                continue
            output_file_path = folder_path / f"otl_target_{folder}_{self.country_name}.tsv"
            if not output_file_path.exists():
                print(f"mozaiko WARNING: No otl_target_{folder}_{self.country_name}.tsv found for {folder}")
                continue

            df = pd.read_csv(output_file_path, sep="\t")
            df = df.drop(columns=['target_family', 'target_genus', 'target_species'], errors='ignore')
            df = df.rename(columns={
                "query_family": "family",
                "query_genus": "genus",
                "query_species": "species",
                "divergence_prct": folder  # Rename divergence_prct to primer name
            })
            tax_res_df.append(df)

        # Merge all dataframes on taxonomy columns
        if tax_res_df:
            merged_df = tax_res_df[0]
            for df in tax_res_df[1:]:
                merged_df = merged_df.merge(
                    df,
                    on=['family', 'genus', 'species'],
                    how='outer'  # outer join keeps all taxonomies from all dataframes
                )

            # The result will have NaN for missing values, which represents empty cells

        merged_df.to_csv(div_score_df_path, sep='\t', index=False, header=True)

    def sort_otl_level_results(self, subdirectory_name: Optional[str] = None):
        """
        This method joins the divergence score results with the OTL-level binding results.
        Files are processed in the base directory, then moved to the specified subdirectory.

        Parameters:
            subdirectory_name (str, optional): Name of the subdirectory where processed files will be moved.
        """
        try:
            # Set up base directory and create subdirectory if needed
            base_directory = Path(self.results_folder) / "otl_level_results"
            base_directory.mkdir(parents=True, exist_ok=True)

            if subdirectory_name:
                target_directory = base_directory / subdirectory_name
                target_directory.mkdir(parents=True, exist_ok=True)

            div_score_df_path = base_directory / "taxonomic_resolution_per_taxon.tsv"

            self.create_taxonomic_resolution_per_taxon(div_score_df_path)

            if not div_score_df_path.exists():
                raise FileNotFoundError(
                    f"mozaiko ERROR: Divergence score file not found: {div_score_df_path}"
                )

            taxonomic_resolution_df = pd.read_csv(div_score_df_path, sep="\t")
            required_columns = ["family", "genus", "species"]

            if not all(
                col in taxonomic_resolution_df.columns for col in required_columns
            ):
                raise ValueError(
                    f"mozaiko WARNING: Missing required columns in divergence score file: {required_columns}"
                )

            for col in required_columns:
                taxonomic_resolution_df[col] = taxonomic_resolution_df[col].astype(str)
                taxonomic_resolution_df[col] = taxonomic_resolution_df[col].replace('nan', '')

            processed_files = []

            # Process files in base directory
            for file in os.listdir(base_directory):
                if file.startswith("otl_level_results_") and file.endswith(".tsv"):
                    file_path = base_directory / file
                    primer_name = (
                        file.replace("otl_level_results_", "")
                        .replace(".tsv", "")
                    )

                    if primer_name not in taxonomic_resolution_df.columns:
                        continue

                    try:
                        binding_df = pd.read_csv(file_path, sep="\t")

                        if not all(
                            col in binding_df.columns for col in required_columns
                        ):
                            continue

                        for col in required_columns:
                            binding_df[col] = binding_df[col].astype(str)
                            binding_df[col] = binding_df[col].replace('nan', '')

                        div_score_subset = taxonomic_resolution_df[
                            required_columns + [primer_name]
                        ].rename(columns={primer_name: "taxonomic_resolution"})

                        merged_df = pd.merge(
                            binding_df,
                            div_score_subset,
                            on=required_columns,
                            how="left",
                        )

                        # Save processed file back to base directory
                        merged_df.to_csv(file_path, sep="\t", index=False)
                        processed_files.append(file)

                        # Move processed file to subdirectory if specified
                        if subdirectory_name:
                            target_path = target_directory / file
                            shutil.move(str(file_path), str(target_path))

                    except Exception as e:
                        print(f"mozaiko ERROR: Error processing {file}: {str(e)}")

            if processed_files:
                # Move the divergence score file to subdirectory if specified
                if subdirectory_name:
                    target_div_score_path = (
                        target_directory / "taxonomic_resolution_per_taxon.tsv"
                    )
                    shutil.move(str(div_score_df_path), str(target_div_score_path))
                else:
                    os.remove(div_score_df_path)

        except Exception as e:
            raise Exception(f"mozaiko ERROR: Error in sort_otl_level_results: {str(e)}")

    def rank_primers_flat(self,
                     save_intermediate_ranks: bool = False,
                     output_path=None,
                     metrics_results_path=None
                     ):
        """
        This method ranks the primers performance based on the results of the Metric System.

        It assigns a ranking order for each relevant metric and sets the ranking order after
        joining all the results. The final rank is the sum of the ranks for each metric. This
        allows for a comprehensive ranking of the primers by setting the metric with the highest
        rank as the one with the lowest value (is first in most of the ranking metrics).

        Parameters:
        - save_intermediate_ranks: bool
            If True, the intermediate ranks will be saved to a TSV file. Default is False.
        - output_path: str
            Path to save the results. If None, the results will be saved to the results folder.
        - metrics_results_path: str
            Path to the dataframe containing the metrics values results. If None, all steps of the
        ranking process will be ran. If a path is provided, the primer ranking will be done using
        the provided dataframe as metrics values.
        """
        self.output_path = output_path

        if metrics_results_path:
            metrics_df = pd.read_csv(metrics_results_path, sep='\t')
            original_metrics = metrics_df.copy()
        else:
            metrics_df = self.join_analysis_results()

        metric_system = {
            "ref_db_qual": {
                "barcoded_taxa": "desc",
                "ratio_barcoded_taxa": "desc"
            },
            "binding_capacity": {
                "normalized_mismatch_score": "asc",
                "normalized_priming_ratio_sum": "asc",
                "normalized_gc_matches_across_taxon": "desc",
                "min_tm_cv": "asc"
            },
            "tax_res": {
                "taxonomic_resolution": "desc",
                "ratio_taxonomic_resolution": "desc"
            }
        }

        # Validation checks
        required_metrics = [
            metric
            for category_metrics in metric_system.values()
            for metric in category_metrics.keys()
        ]

        # 1. Check for missing metric columns
        missing_columns = [m for m in required_metrics if m not in metrics_df.columns]
        if missing_columns:
            raise ValueError(
                f"Primer ranking aborted: missing required metric columns: {missing_columns}"
            )

        # 2. Check for empty values in metric columns
        empty_metrics = {
            metric: metrics_df[metric].isna().sum()
            for metric in required_metrics
            if metrics_df[metric].isna().any()
        }

        if empty_metrics:
            details = ", ".join(
                f"{metric} ({count} NA values)"
                for metric, count in empty_metrics.items()
            )
            raise ValueError(
                f"Primer ranking aborted: empty values detected in metrics: {details}"
            )

        ranking_order = {}
        for category, metrics in metric_system.items():
            for metric, order in metrics.items():
                ranking_order[metric] = order

        for metric, order in ranking_order.items():
            if order == "desc":
                metrics_df[f"rank_{metric}"] = metrics_df[metric].rank(ascending=False)
            elif order == "asc":
                metrics_df[f"rank_{metric}"] = metrics_df[metric].rank(ascending=True)

        # get rank rum and convert to final rank
        rank_sum = metrics_df[[f"rank_{col}" for col in ranking_order]].sum(axis=1)
        metrics_df["final_rank"] = rank_sum.rank(ascending=True).astype(int)
        metrics_df_sorted = metrics_df.sort_values(
            by="final_rank", ascending=True
        ).reset_index(drop=False)
        metrics_df_final = metrics_df_sorted[
            ["primer"] + list(ranking_order.keys()) + ["final_rank"]
        ]

        if output_path is None:
            output_path = Path(self.results_folder) / "ranked_primers.tsv"
        metrics_df_final.to_csv(output_path, sep="\t", index=False)
        print(f"mozaiko INFO: Primer Ranking results saved to {output_path}.")

        # Save intermediate ranks if requested
        if save_intermediate_ranks:
            otl_name = Path(self.otl).stem
            intermediate_ranks_path = (
                Path(self.results_folder) / f"{otl_name}_intermediate_ranks.tsv"
            )
            metrics_df_intermediate = metrics_df_sorted[
                ["primer"] + [f"rank_{col}" for col in ranking_order]
            ]
            metrics_df_intermediate.to_csv(
                intermediate_ranks_path, sep="\t", index=False
            )
            print(
                f"mozaiko INFO: Intermediate ranks saved to {intermediate_ranks_path}."
            )

        return metrics_df_final

    def rank_primers_categorically_weighted(self,
                                            save_intermediate_ranks: bool = False,
                                            output_path=None,
                                            metrics_results_path=None
                                            ):
        """
        This method evaluates primers' performance by ranking them according to the Metric System.
        It calculates an overall rank by summing individual category ranks, prioritizing the
        metric with the lowest value as the highest rank across most metrics.

        Parameters:
        - save_intermediate_ranks: bool
            If True, the intermediate ranks will be saved to a TSV file. Default is False.
        - output_path: str
            Path to save the results. If None, the results will be saved to the results folder.
        - metrics_results_path: str
            Path to the dataframe containing the metrics values results. If None, all steps of the
        ranking process will be ran. If a path is provided, the primer ranking will be done using
        the provided dataframe as metrics values.
        """
        self.output_path = output_path

        if metrics_results_path:
            metrics_df = pd.read_csv(metrics_results_path, sep='\t')
            original_metrics = metrics_df.copy()
        else:
            metrics_df = self.join_analysis_results()

        metric_system = {
            "ref_db_qual": {
                "barcoded_taxa": "desc",
                "ratio_barcoded_taxa": "desc"
            },
            "binding_capacity": {
                "normalized_mismatch_score": "asc",
                "normalized_priming_ratio_sum": "asc",
                "normalized_gc_matches_across_taxon": "desc",
                "min_tm_cv": "asc"
            },
            "tax_res": {
                "taxonomic_resolution": "desc",
                "ratio_taxonomic_resolution": "desc"
            }
        }

        # Validation checks
        required_metrics = [
            metric
            for category_metrics in metric_system.values()
            for metric in category_metrics.keys()
        ]

        # 1. Check for missing metric columns
        missing_columns = [m for m in required_metrics if m not in metrics_df.columns]
        if missing_columns:
            raise ValueError(
                f"Primer ranking aborted: missing required metric columns: {missing_columns}"
            )

        # 2. Check for empty values in metric columns
        empty_metrics = {
            metric: metrics_df[metric].isna().sum()
            for metric in required_metrics
            if metrics_df[metric].isna().any()
        }

        if empty_metrics:
            details = ", ".join(
                f"{metric} ({count} NA values)"
                for metric, count in empty_metrics.items()
            )
            raise ValueError(
                f"Primer ranking aborted: empty values detected in metrics: {details}"
            )

        ranking_order = {}
        for category, metrics in metric_system.items():
            for metric, order in metrics.items():
                ranking_order[metric] = order

        for metric, order in ranking_order.items():
            if order == "desc":
                metrics_df[f"rank_{metric}"] = metrics_df[metric].rank(ascending=False)
            elif order == "asc":
                metrics_df[f"rank_{metric}"] = metrics_df[metric].rank(ascending=True)


        # Calculate rank sum for each category
        for category, metrics in metric_system.items():
            category_metrics = list(metrics.keys())
            category_rank_columns = [f"rank_{metric}" for metric in category_metrics]
            metrics_df[f"rank_sum_{category}"] = metrics_df[category_rank_columns].sum(axis=1)
            metrics_df[f"category_rank_{category}"] = metrics_df[f"rank_sum_{category}"].rank(ascending=True)


        # Calculate final rank as sum of category ranks
        category_rank_columns = [f"category_rank_{category}" for category in metric_system.keys()]
        metrics_df["final_score"] = metrics_df[category_rank_columns].sum(axis=1)
        metrics_df["final_rank"] = metrics_df["final_score"].rank(ascending=True).astype(int)

        metrics_df_sorted = metrics_df.sort_values(
            by="final_rank", ascending=True
        ).reset_index(drop=False)

        output_columns = ["primer"] + list(ranking_order.keys()) + ["final_rank"]
        metrics_df_final = metrics_df_sorted[output_columns]

        if output_path is None:
            output_path = Path(self.results_folder) / "ranked_primers.tsv"
        metrics_df_final.to_csv(output_path, sep="\t", index=False)
        print(f"mozaiko INFO: Primer Ranking results saved to {output_path}.")

        # Save intermediate ranks if requestedresults_folder
        if save_intermediate_ranks:
            otl_name = Path(self.otl).stem
            intermediate_ranks_path = (
                Path(self.results_folder) / f"{otl_name}_intermediate_ranks.tsv"
            )

            metric_rank_columns = [f"rank_{col}" for col in ranking_order]

            category_columns = []
            for category in metric_system.keys():
                category_columns.append(f"rank_sum_{category}")
                category_columns.append(f"category_rank_{category}")

            metrics_df_intermediate = metrics_df_sorted[
                ["primer"] + metric_rank_columns + category_columns + ["final_score", "final_rank"]
            ]

            metrics_df_intermediate.to_csv(
                intermediate_ranks_path, sep="\t", index=False
            )
            print(
                f"mozaiko INFO: Intermediate ranks saved to {intermediate_ranks_path}."
            )

        return metrics_df_final

    @staticmethod
    def evaluate_single_OTL(otl_path,
                            output_folder,
                            primer_table,
                            save_intermediate_ranks=True,
                            run_catnip=True,
                            thresholds: list | float = None,
                            ranking_mode: str = 'category'):
        """
        Evaluate a single OTL file and generate primer rankings.

        Parameters:
        - otl_path: str
            Path to the OTL file
        - output_folder: str
            Output directory for results
        - primer_table: str
            Path to primer table
        - save_intermediate_ranks: bool
            Whether to save intermediate ranking files
        - run_catnip: bool
            Whether to run catnip analysis
        - ranking_mode: str

        Ranking strategy to use. Options:
        - "category": category-weighted ranking
        - "flat": equally weighted per-metric ranking

        Returns:
        - ranked_df: pd.DataFrame
            Ranked primers results
        """
        if ranking_mode not in {"category", "flat"}:
            raise ValueError(
                "ranking_mode must be one of {'category', 'flat'}"
            )

        country_name = Path(otl_path).stem.split('_')[0]
        print("---------------------")
        print(f"mozaiko INFO: Starting evaluating process for {country_name} OTL.")

        incomplete_pbs_path = Path(output_folder) / "incomplete_pbs/filtered"
        catnip_dir = Path(output_folder) / "catnip"

        if run_catnip:
            for file in os.listdir(incomplete_pbs_path):
                if file.endswith(".fasta"):
                    primer_name = os.path.splitext(file)[0]

                    primer_output_dir = os.path.join(catnip_dir, primer_name)

                    # OTL Files
                    expected_output_all = os.path.join(primer_output_dir, f"catnip_target_{primer_name}_{country_name}.tsv")
                    expected_output_otl = os.path.join(primer_output_dir, f"otl_target_{primer_name}_{country_name}.tsv")

                    # Intermediate Results
                    output_interclust = os.path.join(primer_output_dir, f"final_output_interclst.tsv")
                    mapping = os.path.join(primer_output_dir, f"mapping_{primer_name}.tsv")

                    if os.path.exists(expected_output_all) and os.path.exists(expected_output_otl):
                        continue

                    trait = TraitsAndResolution(otl=otl_path, results_folder=output_folder)

                    if os.path.exists(output_interclust) and os.path.exists(mapping):
                        print(f"mozaiko INFO: catnip results already exist for {primer_name}, performing post-processing for {country_name}...")
                        trait.post_process_catnip_primer_results(thresholds)

                    else:
                        print("mozaiko INFO: Taxonomic Resolution files not found. Running catnip...")
                        trait.run_catnip(threshold=thresholds)
                        trait.post_process_catnip_primer_results(thresholds)

        output_path = os.path.join(output_folder, f'{country_name}_ranked_primers.tsv')

        executor = MetricsSystemExecutor(
            results_folder=output_folder,
            otl=otl_path,
            primer_table=primer_table
        )

        if ranking_mode == "category":
            ranked_df = executor.rank_primers_categorically_weighted(
                save_intermediate_ranks=save_intermediate_ranks,
                output_path=output_path
            )
        else:  # "flat"
            ranked_df = executor.rank_primers_flat(
                save_intermediate_ranks=save_intermediate_ranks,
                output_path=output_path
            )

        executor.sort_otl_level_results(subdirectory_name=country_name)

        # return ranked_df

    @staticmethod
    def evaluate_several_OTLs(otl_folder,
                            output_folder,
                            primer_table,
                            save_intermediate_ranks=True,
                            run_catnip=True,
                            thresholds: list | float = None,
                            ranking_mode: str = 'category'):
        """
        Evaluate multiple OTL files in a folder.

        Parameters:
        - otl_folder: str
            Path to folder containing OTL files
        - output_folder: str
            Output directory for results
        - primer_table: str
            Path to primer table
        - save_intermediate_ranks: bool
            Whether to save intermediate ranking files
        - run_catnip: bool
            Whether to run catnip analysis
        - ranking_mode: str
            Ranking strategy to use. Options:
            - "category": category-weighted ranking (default)
            - "flat": equally weighted per-metric ranking


        Returns:
        - results: dict
            Dictionary mapping country names to ranked DataFrames
        """
        results = {}

        # Get all TSV files in the folder
        otl_files = [f for f in os.listdir(otl_folder) if f.endswith('.tsv')]

        if not otl_files:
            print("No TSV files found in the specified folder.")
            return results

        print(f"mozaiko INFO: Found {len(otl_files)} OTL files to process.")

        for otl_file in otl_files:
            otl_path = os.path.join(otl_folder, otl_file)

            try:
                ranked_df = MetricsSystemExecutor.evaluate_single_OTL(
                    otl_path=otl_path,
                    output_folder=output_folder,
                    primer_table=primer_table,
                    save_intermediate_ranks=save_intermediate_ranks,
                    run_catnip=run_catnip,
                    thresholds=thresholds,
                    ranking_mode=ranking_mode
                )

                country_name = Path(otl_path).stem.split('_')[0]
                results[country_name] = ranked_df

            except Exception as e:
                print(f"mozaiko ERROR: Error processing {otl_file}: {str(e)}")
                continue

        print(f"Successfully processed {len(results)} out of {len(otl_files)} OTL files")
        # return results