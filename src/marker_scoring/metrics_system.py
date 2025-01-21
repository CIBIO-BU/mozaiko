import json
import os
import subprocess
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
        4) Create 'species' column populated from 'scientificName' column where 'rank' is 'species',
        'form', 'variety', 'subspecies'

        Returns:
        - otl: DataFrame
            The pre-processed OTL.
        """
        # 1) Transform entries with '-' to NA.
        self.otl = self.otl.replace("-", np.nan)

        # 2) Tranform entries to lower case.
        self.otl["rank"] = self.otl["rank"].str.lower()

        # 3) Remove entries with 'kingdom', 'phylum', 'class', 'order' in 'ranks'.
        ranks = ["kingdom", "phylum", "class", "order"]
        self.otl = self.otl[~self.otl["rank"].isin(ranks)]
        for rank in ranks:
            self.otl = self.otl.dropna(subset=[rank])

        # 4) Remove entries where 'scientificName' is NA.
        self.otl = self.otl.dropna(subset=["scientificName"])

        # 4) Clean ASCII characters from the 'scientificName' column.
        self.otl["scientificName"] = self.otl["scientificName"].apply(
        lambda x: self.fasta_handler.clean_header(x) if pd.notnull(x) else x
        )

        # 5) Create 'species' column populated from 'scientificName' column where 'rank' is
        # 'species',  'form', 'variety', 'subspecies'
        self.otl["rank"] = self.otl["rank"].replace({"form": "species", "variety": "species", "subspecies": "species"})
        self.otl["species"] = np.where(
            self.otl["rank"] == "species", self.otl["scientificName"], np.nan
)
        # 6) Extract the first two strings from the 'species' column
        self.otl["species"] = self.otl["species"].str.split().str[:2].str.join(" ")

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
        unique_otl_taxa = set()

        for entry in otl["scientificName"]:
            unique_otl_taxa.add(entry)

        total_taxa_count = len(unique_otl_taxa)

        self.total_taxa = total_taxa_count
        self.otl_taxa_set = unique_otl_taxa

        # Mapping between scientificName, family, genus, and species
        self.otl_taxa_mapping = {
            taxon: {"family": row["family"], "genus": row["genus"], "species": row["species"], "rank": row["rank"]}
            for _, row in otl.iterrows()
            for taxon in [row["scientificName"]]
        }

        return total_taxa_count, unique_otl_taxa

    def create_taxonomic_hierarchy(self):
        """
        Creates a nested dictionary representing the taxonomic hierarchy from the OTL.
        """
        self.taxa_hierarchy = {}

        otl = self.otl

        # ensure we have the minimum-required columns for taxonomy
        required_cols = ['family', 'genus', 'species']
        if not all(col in otl.columns for col in required_cols):
            raise ValueError("mozaico INFO: OTL must contain family, genus, and species columns.")

        # group by family and genus to create the hierarchy
        for _, row in otl.iterrows():
            family = row['family'] if pd.notna(row['family']) else None
            genus = row['genus'] if pd.notna(row['genus']) else None
            species = row['species'] if pd.notna(row['species']) else None

            if family not in  self.taxa_hierarchy:
                 self.taxa_hierarchy[family] = {'genera': {}, 'count': 0}

            if genus:
                if genus not in  self.taxa_hierarchy[family]['genera']:
                     self.taxa_hierarchy[family]['genera'][genus] = {'species': {}, 'count': 0}

                if species:
                     self.taxa_hierarchy[family]['genera'][genus]['species'][species] = {'count': 0}

        return  self.taxa_hierarchy

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
        This method calculates the total number of barcodes per unique taxonomy entry.

        Parameters:
        - all_inserts_file: File containing all inserts present in the reference database, either
        successfully amplified or not.

        Ouput:
        - barcodes_per_entry: Dictionary containing information on how many barcodes (values)
        exists per taxonomy entry (keys).
        """
        fasta_files = self.parse_fasta_files()
        barcodes_per_entry = {}

        for primer_pair, file_path in fasta_files.items():
            fasta_data = self.custom_fasta_import.read_fasta(
                file_path, check_taxid=False
            )
            # fasta_data.to_csv('data-fasta.csv', sep='\t')
            # print(fasta_data.head())

            barcodes_per_entry[primer_pair] = {}

            for taxa in fasta_data["taxa_info"].unique():
                sequences = fasta_data[fasta_data["taxa_info"] == taxa][
                    "sequence"
                ]

                barcodes_per_entry[primer_pair][taxa.strip()] = len(sequences)

        return barcodes_per_entry

    def calculate_number_of_barcodes_per_otl_taxonomy(
        self,
        barcodes_per_entry: dict,
        otl_hierarchical_taxonomy: dict):
        """
        Calculate barcode counts per taxonomic level, handling missing species/genus data.
        """
        # Replicate OTL taxonomic hierarchy for each primer
        # Allows us to add barcode counts based on OTL structure
        primer_taxa_data = {}
        for primer in barcodes_per_entry.keys():
            primer_taxa_data[primer] = json.loads(json.dumps(otl_hierarchical_taxonomy))

        # Create mapping between parsed taxonomy in headers and their occurences (counts)
        taxa_counts_mapping = {}
        for primer, unique_header_counts in barcodes_per_entry.items():
            for taxa_header, count in unique_header_counts.items():
                taxa_parts = taxa_header.split("|")
                family = taxa_parts[7]
                genus = taxa_parts[8] if len(taxa_parts) > 8 else None
                species = taxa_parts[9] if len(taxa_parts) > 9 else None

                mapping_key = (primer, family, genus, species)
                taxa_counts_mapping[mapping_key] = count

        def update_counts(primer_otl_hierarchy, taxa_counts_mapping, primer):
            """
            Updates counts for each taxonomic level in the OTL-based hierarchy according to the
            number of barcodes.
            """
            for family, family_data in primer_otl_hierarchy.items():
                family_total = 0

                # count entries that only have family level information
                family_only_count = sum(
                    count
                    for (p, f, g, s), count in taxa_counts_mapping.items()
                    if p == primer and f == family and (g == 'nan' or g is None) and (s == 'nan' or s is None)
                )

                if "genera" in family_data:
                    for genus, genus_data in family_data["genera"].items():
                        genus_total = 0

                        # count entries that only have genus level information
                        genus_only_count = sum(
                            count
                            for (p, f, g, s), count in taxa_counts_mapping.items()
                            if p == primer and f == family and g == genus and (s == 'nan' or s is None)
                        )

                        if "species" in genus_data:
                            for species_name, species_data in genus_data["species"].items():
                                # count entries with species level information
                                # excludes 'nan' values from count
                                species_total = sum(
                                    count
                                    for (p, f, g, s), count in taxa_counts_mapping.items()
                                    if p == primer and f == family and g == genus and s and s != 'nan' and species_name.startswith(s)
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

        return primer_taxa_data

    def calculate_percentage_of_taxa_w_x_barcodes(
        self,
        total_taxa_count: int,
        barcode_threshold: int = 1,
    ):
        """
        Calculates the percentage of taxa with more than X barcodes.
        """
        barcodes_per_entry = self.calculate_number_of_barcodes_per_fasta_entry()
        barcodes_per_species = self.calculate_number_of_barcodes_per_otl_taxonomy(
            barcodes_per_entry, self.taxa_hierarchy
        )

        def count_qualifying_taxa(taxa_data, threshold):
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
                            for species_name, species_data in genus_data["species"].items():
                                if isinstance(species_data, dict) and species_data.get("count", 0) >= threshold:
                                    qualifying_taxa += 1
                                    genus_has_qualifying_descendant = True
                                    family_has_qualifying_descendant = True

                        if not genus_has_qualifying_descendant and isinstance(genus_data, dict) and genus_data.get("count", 0) >= threshold:
                            qualifying_taxa += 1
                            family_has_qualifying_descendant = True

                if not family_has_qualifying_descendant and isinstance(family_data, dict) and family_data.get("count", 0) >= threshold:
                    qualifying_taxa += 1

            return qualifying_taxa

        results = {}
        for primer_pair, taxa_data in barcodes_per_species.items():
            taxa_meeting_threshold = count_qualifying_taxa(taxa_data, barcode_threshold)
            percentage = (taxa_meeting_threshold / total_taxa_count) * 100
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
        - barcoded_taxa_ratio: Dict
            A dictionary containing the raatio of barcoded taxa and the percentage of taxa with more than five barcodes per primer pair.
        """
        barcoded_taxa_five_plus = self.calculate_percentage_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=5
        )
        barcoded_taxa_one_plus = self.calculate_percentage_of_taxa_w_x_barcodes(
            total_taxa_count, barcode_threshold=1
        )

        barcoded_taxa_ratio = {}

        for primer_pair in barcoded_taxa_five_plus.keys():
            percent_5plus = barcoded_taxa_five_plus[primer_pair]
            percent_1plus = barcoded_taxa_one_plus[primer_pair]

            ratio_barcoded_taxa = (
                percent_5plus / percent_1plus if percent_5plus > 0 else 0
            )

            barcoded_taxa_ratio[primer_pair] = {
                "barcoded_taxa_one_plus": percent_1plus,
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
        self.otl_taxa_mapping = self.otl_handler.otl_taxa_mapping
        self.otl_unique_taxa = self.otl_handler.otl_taxa_set

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
                    pbs_table = self.get_pbs_table(amplicon_file, insert_file)

                    for _, pbs_row in pbs_table.iterrows():
                        seq_header = pbs_row["header"].replace(">", "")
                        seq_id = seq_header.split("|")[0].replace(" ", "")
                        taxon = seq_header.split("|")[2]
                        rank = seq_header.split("|")[3]
                        family = seq_header.split("|")[8]
                        genus = seq_header.split("|")[9]
                        species = seq_header.split("|")[10]
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
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
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
                                "family": family,
                                "genus": genus,
                                "species": species,
                                "rank": rank,
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
                    "family",
                    "genus",
                    "species",
                    "rank",
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
        primer_gc_df.set_index("primer_name", inplace=True)

        if save_results:
            primer_gc_df.to_csv("primer_gc_fractions.csv", index=False)

        return primer_pbs_df, primer_gc_df

    def add_missing_otl_taxa(self, primer_df):
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
        missing_taxon = self.otl_unique_taxa - taxon_on_df
        # intersection = taxon_on_df & missing_taxon
        # len_intersection = len(intersection)

        if missing_taxon:
            missing_data = {
                "taxon": list(missing_taxon),
                "seq_id": ["otl-import"] * len(missing_taxon),
            }

            # Add nan values for metrics columns
            metrics_columns = {
                col: np.nan
                for col in primer_df.columns
                if col not in ["taxon", "seq_id", "family", "genus", "species", "rank"]
            }

            # Add information for taxonomy columns based on the OTL
            taxonomy_columns = {
                "family": [self.otl_handler.otl_taxa_mapping[taxon]["family"] for taxon in missing_taxon],
                "genus": [self.otl_handler.otl_taxa_mapping[taxon]["genus"] for taxon in missing_taxon],
                "species": [self.otl_handler.otl_taxa_mapping[taxon]["species"] for taxon in missing_taxon],
                "rank": [self.otl_handler.otl_taxa_mapping[taxon]["rank"] for taxon in missing_taxon],
            }

            new_entries = pd.DataFrame({**missing_data, **metrics_columns, **taxonomy_columns})

            otl_populated_df = pd.concat([primer_df, new_entries], ignore_index=True)
        else:
            otl_populated_df = primer_df

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
                self.add_missing_otl_taxa(
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

        # Add missing taxa to the dataframe with values of Nan
        # primer_df = self.add_missing_otl_taxa(primer_df)

        grouped_taxa = primer_df.groupby(["family", "genus", "species"])

        if operation in {"min", "max", "sum", "mean"}:
            result = getattr(grouped_taxa[analysis_name], operation)().astype(float)
        elif operation == "coef_var":
            mean = grouped_taxa[analysis_name].mean()
            std = grouped_taxa[analysis_name].std()
            result = (std / mean.replace(0, pd.NA)) * 100

        else:
            raise ValueError(
                f"mozaiko ERROR: Unrecognized operation: '{operation}'. "
                "Please choose from 'min', 'max', 'sum', 'mean', or 'coef_var'."
            )

        result = pd.DataFrame(result)
        result.reset_index(inplace=True)
        result[['family', 'genus', 'species']] = result[['family', 'genus', 'species']].replace('nan', np.nan)
        # ref_otl = self.otl_handler.otl[['family', 'genus', 'species']]
        # ref_otl = ref_otl.drop_duplicates()

        # # Get entries only for the taxa that are in the OTL and add taxa that are missing (left join)
        # otl_based_result = pd.merge(
        #     result,
        #     ref_otl,
        #     on=["family", "genus", "species"],
        #     how="right"
        # )
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
        in_silico_amplified_inserts = os.path.join(results_folder, "insert/filtered")
        # Input B
        all_inserts_with_pbs = os.path.join(
            results_folder, "all_complete_pbs/filtered/filtered_intersection"
        )
        # Input C
        inserts_with_incomplete_pbs = os.path.join(
            results_folder, "incomplete_pbs/filtered/filtered_intersection"
        )

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


class TraitsAndResolution:
    def __init__(
        self,
        results_folder: Optional[str] = None,
        insert_folder_path: Optional[str] = None,
        amplicon_folder_path: Optional[str] = None,
        incomplete_pbs_folder_path: Optional[str] = None,
    ):
        print(
            f"results_folder: {results_folder}, insert_folder_path: {insert_folder_path}, amplicon_folder_path: {amplicon_folder_path}"
        )
        if results_folder is not None:
            self.results_folder = results_folder
            self.insert_folder_path = os.path.join(results_folder, "insert/filtered")
            self.amplicon_folder_path = os.path.join(
                results_folder, "amplicon/filtered"
            )
            self.incomplete_pbs_path = os.path.join(
                results_folder, "incomplete_pbs/filtered/filtered_intersection"
            )
            print(
                f"Set insert_folder_path to {self.insert_folder_path}, amplicon_folder_path to {self.amplicon_folder_path} and incomplete_pbs_path to {self.incomplete_pbs_path}."
            )
        elif (
            insert_folder_path is not None
            and amplicon_folder_path is not None
            and incomplete_pbs_folder_path is not None
        ):
            self.insert_folder_path = insert_folder_path
            self.amplicon_folder_path = amplicon_folder_path
            self.incomplete_pbs_path = incomplete_pbs_folder_path
            print(
                "Using provided insert_folder_path, amplicon_folder_path and incomplete_pbs_folder_path."
            )
        else:
            raise ValueError(
                "mozaiko ERROR: Either provide a path to the in-silico amplification results folder "
                "('results_folder') or paths to the insert ('insert_folder_path'), amplicon "
                "('amplicon_folder_path') and incomplete PBS ('incomplete_pbs_path') results folders."
            )

    def get_min_max_avg_seq_length_in_a_fasta(self, fasta_file):
        """
        This method calculates the  minimum, maximum, and average sequence lengths in a FASTA file.

        Parameters:
        - fasta_file: str
            Path to the FASTA file

        Returns:
            tuple: (min_length, max_length, avg_length)
        """
        try:
            seq_lengths = [
                len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")
            ]

            if not seq_lengths:
                return np.nan, np.nan, np.nan

            min_length = min(seq_lengths)
            max_length = max(seq_lengths)
            avg_length = sum(seq_lengths) / len(seq_lengths)

            return min_length, max_length, round(avg_length, 2)

        except Exception as e:
            print(f"mozaiko ERROR: Error processing {fasta_file}: {e}")
            return np.nan, np.nan, np.nan

    def get_length_stats_for_amplicon_and_insert(
        self, insert_folder_path=None, amplicon_folder_path=None
    ):
        """
        This method analyzes sequence lengths for insert and amplicon FASTA files.

        Parameters:
        - insert_folder_path (str, optional):
            Path to insert FASTA files
        - amplicon_folder_path (str, optional):
            Path to amplicon FASTA files

        Returns:
            pd.DataFrame: DataFrame with sequence length stats
        """
        if insert_folder_path is None:
            insert_folder_path = self.insert_folder_path
        if amplicon_folder_path is None:
            amplicon_folder_path = self.amplicon_folder_path

        if not (insert_folder_path and amplicon_folder_path):
            raise ValueError(
                "mozaiko ERROR: Insert or amplicon folder paths are not specified."
            )

        length_data = {"insert": {}, "amplicon": {}}

        for root, dirs, files in os.walk(insert_folder_path):
            for file in files:
                if file.endswith(".fasta"):
                    primer_name = os.path.splitext(file)[0]
                    file_path = os.path.join(root, file)

                    min, max, avg = self.get_min_max_avg_seq_length_in_a_fasta(
                        file_path
                    )

                    length_data["insert"][primer_name] = {"avg_length": avg}

        for root, dirs, files in os.walk(amplicon_folder_path):
            for file in files:
                if file.endswith(".fasta"):
                    primer_name = os.path.splitext(file)[0]
                    file_path = os.path.join(root, file)

                    min, max, avg = self.get_min_max_avg_seq_length_in_a_fasta(
                        file_path
                    )

                    length_data["amplicon"][primer_name] = {
                        "min_length": min,
                        "max_length": max,
                        "avg_length": avg,
                    }

        insert_df = pd.DataFrame.from_dict(length_data["insert"], orient="index")
        amplicon_df = pd.DataFrame.from_dict(length_data["amplicon"], orient="index")

        result_df = pd.concat(
            [insert_df.add_prefix("insert_"), amplicon_df.add_prefix("amplicon_")],
            axis=1,
        )

        return result_df.fillna(np.nan)

    def run_multibarcode_pipeline(self):
        """
        Runs MultiBarcodePipeline on the insert files.

        Zhu, T., & Iwasaki, W. (2023). MultiBarcodeTools: Easy selection of optimal primers for
        eDNA multi-metabarcoding. Environmental DNA, 5, 1793-1808. https://doi.org/10.1002/edn3.499
        """
        results_folder_base = os.path.dirname(self.insert_folder_path)
        multibarcode_output_folder = os.path.join(results_folder_base, "multibarcode")
        self.multibarcode_output_folder = multibarcode_output_folder
        os.makedirs(multibarcode_output_folder, exist_ok=True)

        multibarcode_outputdir = os.path.join(
            multibarcode_output_folder, "multibarcode_input.tsv"
        )

        multibarcode_file = create_MultiBarcodeTools_input(
            self.insert_folder_path, self.incomplete_pbs_path, multibarcode_outputdir
        )

        script_dir = os.path.dirname(os.path.abspath(__file__))
        multibarcode_script = os.path.join(script_dir, "multibarcodepipeline.sh")
        os.chmod(multibarcode_script, 0o755)

        try:
            result = subprocess.run(
                [
                    multibarcode_script,
                    self.multibarcode_output_folder,
                    multibarcode_file,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            # print(result.stdout)
            print(
                f"mozaiko INFO: MultiBarcodePipeline completed. Output in {self.multibarcode_output_folder}"
            )
            return result.stdout

        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: Error running MultiBarcodePipeline - {e}")
            print(f"Command output: {e.output}")
            raise

    def parse_multibarcode_output(self, result_stdout):
        """
        Method to parse MultiBarcode stdout into a DataFrame.

        Parameters:
            result_stdout: str
        Output string from MultiBarcode

        Returns:
            pandas.DataFrame: Parsed primer information
        """
        primers = []
        species_counts = []
        cumulative_counts = []

        # Regular expression to extract primer and species count
        pattern = r"(\d+):\s*([\w-]+),\s*resolve\s*(?:additional\s*)?(\d+)\s*species"

        cumulative = 0
        for match in re.finditer(pattern, result_stdout):
            primer = match.group(2)
            count = int(match.group(3))

            primers.append(primer)
            species_counts.append(count)

            cumulative += count
            cumulative_counts.append(cumulative)

        primer_resolv_species = pd.DataFrame(
            {
                "primer": primers,
                "additional_resolved_species": species_counts,
                "cumulative_resolved_species": cumulative_counts,
            }
        )

        self.primer_resolv_species = primer_resolv_species

        return self.primer_resolv_species

    def load_nucleotide_distance(self):
        """
        This method loads the nuclotide difference matrix into memory as outputed by MultiBarcode
        pipeline.
        """
        nuc_dist_matrix_path = os.path.join(
            self.multibarcode_output_folder, "matrix.xlsx"
        )
        nuc_dist_matrix = pd.read_excel(nuc_dist_matrix_path)

        return nuc_dist_matrix

    def compute_genetic_divergence_per_taxon(self):
        """
        Calculate the percentage of genetic divergence for each taxon.

        Returns:
        pd.DataFrame: Percentage of genetic divergence for each taxon
        """
        nuc_dist_matrix = self.load_nucleotide_distance()
        self.nuc_dist_matrix = nuc_dist_matrix
        insert_amplicon_len_stats = self.get_length_stats_for_amplicon_and_insert()
        insert_avg_len_matrix = insert_amplicon_len_stats["insert_avg_length"]

        # Create mapping to harmonize primer naming format between MultiBarcode and
        # Insert Length Matrix
        name_mapping = {
            col.replace("-", "_"): col for col in insert_avg_len_matrix.index
        }

        divergence_percentages = {}

        for taxon in nuc_dist_matrix.index:
            taxa_distances_per_primer = nuc_dist_matrix.loc[taxon]
            taxon_divergence = {}

            # Include nuc_dist_matrix species as index
            taxon_divergence["Species"] = taxa_distances_per_primer["Species"]

            for primer_set, nucleotide_distance in taxa_distances_per_primer.items():

                if primer_set == "Species":  # Skip processing for 'Species' column
                    continue

                try:
                    dist_numeric = float(nucleotide_distance)

                    # Match nuc_dist_matrix primer name with internal dataframe format
                    mapped_primer_set = name_mapping.get(primer_set, primer_set)

                    if mapped_primer_set in insert_avg_len_matrix.index:
                        # Convert insert length to float to handle potential division by zero
                        insert_length = float(insert_avg_len_matrix[mapped_primer_set])

                        if insert_length > 0:
                            divergence_percentage = round(
                                ((dist_numeric / insert_length) * 100), 1
                            )
                        else:
                            divergence_percentage = np.nan
                    else:
                        divergence_percentage = np.nan

                except (ValueError, TypeError):
                    # If conversion fails or any other error occurs set no Nan
                    divergence_percentage = np.nan

                taxon_divergence[primer_set] = divergence_percentage

            divergence_percentages[taxon] = taxon_divergence

        divergence_df = pd.DataFrame.from_dict(divergence_percentages, orient="index")

        # Reorder to have 'Species' as first col
        columns = ["Species"] + [
            col for col in divergence_df.columns if col != "Species"
        ]
        divergence_df = divergence_df[columns]

        return divergence_df

    def get_divergence_score(self, total_otl_taxa_count: int, cutoff: float = 2.0):
        """
        This method retrieves the divergence score for each primer set. The divergence score is
        the percentage of taxa with a percentage of divergence bellow a cutoff according to the
        total number of taxa in the OTL.

        Parameters:
        - otl_total_taxa_cont: int
            The total number of taxa in the OTL.
        - cutoff: float
            The cutoff for genetic divergence. Float number representing a percentage. Default is
            2.0, which is equivalent to 2.0%.

        Returns:
        """
        if total_otl_taxa_count <= 0:
            raise ValueError(
                "mozaiko ERROR: The total number of taxa in OTL must be above 0."
            )

        self.divergence_df = self.compute_genetic_divergence_per_taxon()

        primer_cols = self.divergence_df.loc[:, self.divergence_df.columns != "Species"]

        divergence_score_results = []

        for column in primer_cols:
            taxa_above_cutoff = self.divergence_df[self.divergence_df[column] > cutoff][
                column
            ].count()

            divergence_score = (taxa_above_cutoff / total_otl_taxa_count) * 100

            divergence_score_results.append(
                {
                    "primer": column,
                    "total_taxa": total_otl_taxa_count,
                    "n_taxa_above_cutoff": taxa_above_cutoff,
                    "divergence_score": round(divergence_score, 2),
                }
            )

        divergence_score_df = pd.DataFrame(divergence_score_results)

        return divergence_score_df


class MetricsSystemExecutor:
    """
    This class orchestrates the entire Metrics System, coordinating the execution of all evaluation
    categories by calling the necessary classes and methods.
    """

    def __init__(self, results_folder: str, otl: str, primer_table: str):
        # Initialize Reference Database Quality Category
        #self.ref_db = ReferenceDatabaseQuality(otl=otl)
        # Initialize OTL and related variables
        self.otl = otl
        self.otl_handler = OtlHandler(self.otl)
        self.otl_handler.import_otl()
        self.total_otl_taxa_count = self.otl_handler.total_taxa
        self.otl_unique_taxa_set = self.otl_handler.otl_taxa_set
        # Load files and folders into memory
        self.primer_table = primer_table
        self.results_folder = results_folder
        self.results_folder = results_folder
        self.insert_folder_path = os.path.join(results_folder, "insert/filtered")
        self.amplicon_folder_path = os.path.join(results_folder, "amplicon/filtered")
        self.incomplete_pbs_path = os.path.join(
            results_folder, "incomplete_pbs/filtered/filtered_intersection"
        )
        print(
            f"Set insert_folder_path to {self.insert_folder_path}, amplicon_folder_path to {self.amplicon_folder_path} and incomplete_pbs_path to {self.incomplete_pbs_path}."
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
        red_bd_qual = ReferenceDatabaseQuality(otl=self.otl, all_inserts_path=self.insert_folder_path)
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

    def comprehensive_primer_analysis(self, output_folder):
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

        primer_pbs, gc_df = self.get_primer_pbs_analysis()

        primer_results = {}

        ref_qual = self.get_reference_database_quality()

        primers_to_analyze = list(primer_pbs.keys())

        for primer in primers_to_analyze:
            primer_metrics = {}

            # Mismatch Analyses
            tax_lev_max_ms_full_len = binding.process_analysis_per_taxon(
                primer_pbs[primer],
                operation="max",
                analysis_name="full_len_mismatch_sum",
            )
            primer_metrics["max_mismatch_across_taxon"] = int(
                binding.process_analysis_across_taxon(
                    tax_lev_max_ms_full_len, operation="sum"
                )
            )

            # Priming Ratio
            tax_lev_max_ms_three_end = binding.process_analysis_per_taxon(
                primer_pbs[primer],
                operation="max",
                analysis_name="three_end_mismatch_sum",
            )
            primer_metrics["priming_ratio"] = binding.get_priming_ratio(
                tax_lev_max_ms_full_len, tax_lev_max_ms_three_end
            )

            # GC Match Analysis
            binding.get_total_gc_matches(primer_pbs[primer])
            tax_lev_gc = binding.process_analysis_per_taxon(
                primer_pbs[primer], operation="min", analysis_name="gc_matches_score"
            )
            primer_metrics["gc_matches_across_taxon"] = (
                binding.process_analysis_across_taxon(tax_lev_gc, operation="sum")
            )

            # Temperature Melting (Tm) Analysis
            tax_lev_min_tm = binding.process_analysis_per_taxon(
                primer_pbs[primer], operation="min", analysis_name="min_tm"
            )
            primer_metrics["tm_coefficient_var"] = (
                binding.process_analysis_across_taxon(
                    tax_lev_min_tm, operation="coef_var"
                )
            )

            # Tm Score
            primer_metrics["tm_score"] = binding.tm_score(primer_pbs[primer])

            # Amplification Success (if applicable)
            try:
                amp_succ = binding.calculate_amplification_success_score(output_folder)

                amp_succ_value = (
                    amp_succ.loc[primer, "amplification_sucess_percent"]
                    if primer in amp_succ.index
                    else None
                )
                primer_metrics["amplification_success_percent"] = amp_succ_value

            except Exception as e:
                print(
                    f"mozaico WARNING: Amplification success calculation failed for {primer}: {e}"
                )
                primer_metrics["amplification_success_percent"] = None

            # Store results for this primer
            primer_results[primer] = primer_metrics

        binding_df = pd.DataFrame.from_dict(primer_results, orient="index")

        if ref_qual is not None:
            binding_df = ref_qual.join(binding_df)

        binding_df.rename(index=lambda x: x.replace("-", "_"), inplace=True)
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
        trait = TraitsAndResolution(results_folder=self.results_folder)

        trait.multibarcode_output_folder = os.path.join(
            self.results_folder, "multibarcode"
        )

        # Run Multibarcode Pipeline
        output_str = trait.run_multibarcode_pipeline()
        trait.parse_multibarcode_output(output_str)

        # # Get Taxonomic Resolution Percentage
        # taxonomic_resolution = trait.get_taxonomic_resolution(
        #     total_otl_taxa_count=int(self.total_otl_taxa_count)
        # )

        # Get Divergence Score
        divergence_score = trait.get_divergence_score(
            total_otl_taxa_count=int(self.total_otl_taxa_count)
        )

        if "primer" in divergence_score.columns:
            divergence_score = divergence_score.set_index("primer")

        combined_div_score_and_tax_res_results = pd.DataFrame(
            {"divergence_score": divergence_score["divergence_score"]}
        )

        traits_res_df = combined_div_score_and_tax_res_results

        return traits_res_df

    def join_analysis_results(self):
        """
        This method joins the results from the primer analysis and the traits and resolution analysis.
        """
        binding_dataframe = self.comprehensive_primer_analysis(self.results_folder)
        traits_dataframe = self.get_traits_and_resolution()

        analysis_results = binding_dataframe.join(traits_dataframe, on="primer")

        return analysis_results

    def rank_primers(self, save_intermediate_ranks: bool = False, output_path=None):
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
        """
        ranking_order = {
            "barcoded_taxa_one_plus": "desc",
            "ratio_barcoded_taxa": "desc",
            "max_mismatch_across_taxon": "asc",
            "priming_ratio": "asc",
            "gc_matches_across_taxon": "desc",
            "tm_coefficient_var": "asc",
            "tm_score": "desc",
            "amplification_success_percent": "desc",
            "divergence_score": "asc",
        }

        metrics_df = self.join_analysis_results()

        for column, order in ranking_order.items():
            if order == "desc":
                metrics_df[f"rank_{column}"] = metrics_df[column].rank(ascending=False)
            elif order == "asc":
                metrics_df[f"rank_{column}"] = metrics_df[column].rank(ascending=True)

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
