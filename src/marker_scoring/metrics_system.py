import json
import os
import sys
from collections import defaultdict
from typing import Dict

import pandas as pd

from src.in_silico_analysis.amplification import InSilicoAmplification
from src.marker_scoring.scoring_utils import *
from src.reference_database.sequence_import import CustomFastaImport


class ReferenceDatabaseQuality:
    def __init__(self, all_inserts_file=None, otl=None):
        self.custom_fasta_import = CustomFastaImport()
        self.all_inserts_file = all_inserts_file
        self.otl = otl
        self.total_taxa = None

    def _validate_otl(self, otl=None):
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

        if self.all_inserts_file is None:
            print(
                "mozaiko INFO: To continue the evaluation, please provide a FASTA file containing \
                    all insert regions present in the original database, whether successfully \
                        amplified or not."
            )
            all_inserts_file = input("Please enter the path to a FASTA file: ")
            self.all_inserts_file = str(all_inserts_file)

        fasta_data = self.custom_fasta_import.read_fasta(
            self.all_inserts_file, check_taxid=False
        )

        barcodes_per_species = {}

        for taxa in fasta_data["taxa_info"].unique():
            unique_sequences = fasta_data[fasta_data["taxa_info"] == taxa][
                "sequence"
            ].unique()

            barcodes_per_species[taxa] = len(unique_sequences)

        return barcodes_per_species

    def calculate_percentage_of_taxa_w_x_barcodes(
        self, barcode_threshold=1, total_taxa=None
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

        threshold_valid_taxa = 0

        for key, value in barcodes_per_species.items():
            if value > barcode_threshold:
                threshold_valid_taxa += 1

        if total_taxa is not None:
            percentage_of_taxa_w_x_barcodes = (threshold_valid_taxa * 100) / total_taxa
        else:
            percentage_of_taxa_w_x_barcodes = (
                threshold_valid_taxa * 100
            ) / self.total_taxa
        percentage_of_taxa_w_x_barcodes = round(percentage_of_taxa_w_x_barcodes, 2)

        return percentage_of_taxa_w_x_barcodes

    def import_otl(self):
        """ """
        if self.otl is None:
            print(
                "mozaiko INFO: To continue the evaluation, a Operational Taxonomic List (OTL) is \
                    required. An OTL is a list contaning information on the taxonomic numenclature \
                        of all identifiable taxa in routine biomonitoring initiatives."
            )
            self.otl = input("Please enter the path to the OTL: ")

        self._validate_otl()

        otl = self.otl
        unique_taxa = set()

        for entry in otl["taxa"]:
            unique_taxa.add(entry)

        total_taxa = len(unique_taxa)

        self.total_taxa = total_taxa

    def ratio_barcoded_taxa(self):
        """
        This method calculates the Ratio of Barcoded Taxa (RBT).

        Parameters:
        - barcoded_taxa_one: float
            Percentage of taxa with more than five barcodes.
        - barcoded_taxa_two: float
            Percentage of taxa with more than one barcode.

        Output:
        - rbt: float
            Barcode Coverage Score
        """
        self.import_otl()

        barcoded_taxa_one = self.calculate_percentage_of_taxa_w_x_barcodes(
            barcode_threshold=5
        )
        barcoded_taxa_two = self.calculate_percentage_of_taxa_w_x_barcodes(
            barcode_threshold=1
        )

        rbt = barcoded_taxa_two / barcoded_taxa_one

        rbt_rounded = round(rbt, 2)

        return barcoded_taxa_one, rbt_rounded


class Binding:
    def __init__(self, number_of_mismatches=None):
        self.amplification_instance = InSilicoAmplification()
        if number_of_mismatches is None:
            number_of_mismatches = (
                self.amplification_instance.get_number_of_mismatches()
            )
        self.number_of_mismatches = number_of_mismatches

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

    def get_number_of_primer_pbs_mismatches(
        self,
        amplicon_folder,
        insert_folder,
        primer_table,
        only_at_three_end: bool = False,
        save_results: bool = False,
    ):
        """
        This method retrieves the number of mismatches between the foward and reverse primer-PBS
        sequence pair.

        It first retrieves the primer and PBS table. The last by iterating over the generated files
        in the analysis folder. It then processes both tables to retrieve the primer and primer-binding
        site sequences. Lastly it calls the method to calculate the number of mismatches and creates
        a dictionary with the results.

        Parameters:
        - amplicon_folder: Output folder from the in-silico amplification process. Its files should
        contain the amplicon sequences.
        - insert_folder: Output folder from the in-silico amplification process. Its files should
        contain the insert sequences.
        - primer_table: Path to the table containing the primer pair sequences and details.
        - only_at_three_end: Bool
            When set to True it calculates the number of mismatches only at the three end section
            of the sequences (last five nucleotides). If set to False, it calculates the number
            of mismatches at the overall length.
        - save_results: bool
            When set to True, the dictionary containing the results will be saved as a JSON file.

        Output:
        - max_mismatches_per_taxon: Dict
            A nested dictionary containing the primer pairs as keys and a dictionary as values,
            containing the sum of mismatches between the forward and reverse primer-PBS sequences.
        """
        self.get_primer_table(primer_table)
        matching_files = self.parse_files_with_same_extension_in_folders(
            amplicon_folder, insert_folder
        )

        if not matching_files:
            return None

        primer_pbs_mismatches = {}

        for primer_ind, primer_row in self.primer_table.iterrows():
            barcode_region = primer_row["barcode_region"]
            assay_name = primer_row["assay_name"]
            pbs_filename = barcode_region + "_" + assay_name
            primer_seq_fwd = primer_row["fwd_seq"]
            primer_seq_rev = primer_row["rev_seq"]

            primer_results = []

            for amplicon_file, insert_file in matching_files:
                amplicon_filename = os.path.splitext(os.path.basename(amplicon_file))[0]

                if pbs_filename == amplicon_filename:
                    pbs_table = self.get_pbs_table(amplicon_file, insert_file)

                    for pbs_ind, pbs_row in pbs_table.iterrows():
                        seq_header = pbs_row["header"]
                        seq_header = seq_header.replace(">", "")
                        seq_id = seq_header.split("|")[0]
                        seq_id = seq_id.replace(" ", "")
                        taxon = seq_header.split("|")[1]

                        pbs_fwd_seq = pbs_row["fwd_seq"]
                        pbs_rev_seq = pbs_row["rev_seq"]

                        rev_comp_pbs_rev_seq = str(
                            Seq(pbs_rev_seq).reverse_complement()
                        )

                        if only_at_three_end == True:
                            primer_seq_fwd = primer_seq_fwd[-5:]
                            primer_seq_rev = primer_seq_rev[-5:]
                            pbs_fwd_seq = pbs_fwd_seq[-5:]
                            rev_comp_pbs_rev_seq = rev_comp_pbs_rev_seq[-5:]

                        fwd_primer_pbs_mismatches = calculate_iupac_mismatches(
                            primer_seq_fwd, pbs_fwd_seq
                        )
                        rev_primer_pbs_mismatches = calculate_iupac_mismatches(
                            primer_seq_rev, rev_comp_pbs_rev_seq
                        )

                        mismatch_sum = (
                            fwd_primer_pbs_mismatches + rev_primer_pbs_mismatches
                        )

                        seq_result = {
                            "seq_id": seq_id,
                            "taxon": taxon,
                            "primer_seq_fwd": primer_seq_fwd,
                            "primer_seq_rev": primer_seq_rev,
                            "pbs_fwd_seq": pbs_fwd_seq,
                            "rev_comp_pbs_rev_seq": rev_comp_pbs_rev_seq,
                            "mismatch_sum": mismatch_sum,
                        }
                        primer_results.append(seq_result)

            if primer_results:
                primer_pbs_mismatches[pbs_filename] = primer_results

        if save_results == True:
            with open("primer_pbs_mismatches.json", "w") as fp:
                json.dump(primer_pbs_mismatches, fp)

        return primer_pbs_mismatches

    def get_max_mismatches_per_taxon(
        self, mismatches_dictionary: dict, save_results: bool = False
    ):
        """
        This method get the maximum number of mismatches per taxon for each primer pair.

        Parameters:
        - mismatches_dicitonary: Dict
            A dictionary containing the number of mismatches for every sequence for each primer pair.
        - save_results: bool
            When set to True, the dictionary containing the results will be saved as a JSON file.

        Ouput:
         - max_mismatches_per_taxon: Dict
            A nested dictionary containing the primer pairs as keys and a dictionary containing the
            max number of mismatches attained for each taxon.
        """
        grouped_taxons: Dict = {}

        for primer_pair, entries in mismatches_dictionary.items():
            grouped_taxons[primer_pair] = {}
            for entry in entries:
                taxon = entry["taxon"].strip()

                if taxon not in grouped_taxons[primer_pair]:
                    grouped_taxons[primer_pair][taxon] = []

                grouped_taxons[primer_pair][taxon].append(entry["mismatch_sum"])

        max_mismatches_per_taxon = {}

        for primer_pair, taxon_data in grouped_taxons.items():
            max_mismatches_per_taxon[primer_pair] = {
                taxon: max(mismatches) for taxon, mismatches in taxon_data.items()
            }

        if save_results == True:
            with open("max_mismatches_per_taxon.json", "w") as fp:
                json.dump(max_mismatches_per_taxon, fp)

        return max_mismatches_per_taxon

    def get_max_mismatches_count_across_taxon(
        self, max_mismatches_per_taxon_dict: dict
    ):
        """
        This method retrieves the sum of the maximum mismatches across taxon.

        Parameters:
        - max_mismatches_per_taxon_dict: Dict
            A nested dictionary containing the primer pairs as keys and a dictionary containing the
            max number of mismatches attained for each taxon.

        Output:
        - max_mismatches_count_per_primer: Dict
            A dictionary containing the primer pairs as keys and the total sum of the maximum number
            of mismatches from all taxons identified by the primer.
        """
        max_mismatches_count_per_primer = {}

        for primer_pair, max_taxon_mismatches in max_mismatches_per_taxon_dict.items():
            max_mismatches_count_per_primer[primer_pair] = sum(
                max_taxon_mismatches.values()
            )

        return max_mismatches_count_per_primer

    def get_priming_ration(
        self, mismatches_dictionary_all_len: dict, mismatches_dictionary_three_end: dict
    ):
        """
        This method computes the priming ration between the maximum number of mismatches per taxon
        and the maximum number of mismatches at the 3'end per taxon.

        Parameters:
        - mismatches_dictionary_all_len: Dict
            A nested dictionary containing the primer pairs as keys and a dictionary containing the
            max number of mismatches attained for each taxon, in across the total sequence length.
        - mismatches_dictionary_three_end: Dict
            A nested dictionary containing the primer pairs as keys and a dictionary containing the
            max number of mismatches attained for each taxon, in across the 3'end length.

        Output:
        - priming_ratio_dict: Dict
            A dictionary containing the sum of the ratios between the maximum number of mismatches
            per taxon and the maximum number of mismatches at the 3'end per taxon.
        """
        priming_ratio_dict = {}

        for primer_pair in mismatches_dictionary_three_end:
            if primer_pair in mismatches_dictionary_all_len:
                priming_ratio_dict[primer_pair] = 0

                three_end_taxa = set(
                    mismatches_dictionary_three_end[primer_pair].keys()
                )
                all_len_taxa = set(mismatches_dictionary_all_len[primer_pair].keys())
                common_taxa = three_end_taxa.intersection(all_len_taxa)

                for taxon in common_taxa:
                    three_end_value = mismatches_dictionary_three_end[primer_pair][
                        taxon
                    ]
                    all_len_value = mismatches_dictionary_all_len[primer_pair][taxon]

                    if all_len_value != 0:
                        ratio = three_end_value / all_len_value
                        priming_ratio_dict[primer_pair] += ratio

        for key, value in priming_ratio_dict.items():
            rounded_value = round(value, 2)
            priming_ratio_dict[key] = rounded_value

        return priming_ratio_dict


class MetricsSystemExecutor:
    """
    This class orchestrates the entire Metrics System, coordinating the execution of all evaluation
    categories by calling the necessary classes and methods.
    """

    def __init__(self, all_inserts_folder=None, otl=None):
        self.ref_db = ReferenceDatabaseQuality()
        self.otl = otl
        self.all_inserts_folder = all_inserts_folder

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
        ref_bd_scores = {}

        for filename in os.listdir(self.all_inserts_folder):
            if filename.endswith(".fasta"):
                fasta_file = os.path.join(self.all_inserts_folder, filename)
                primer_name = filename.split(".")[0]
                cls = ReferenceDatabaseQuality(fasta_file, self.otl)
                bt_1, rbt = cls.ratio_barcoded_taxa()
                ref_bd_scores[primer_name] = bt_1, rbt

        return ref_bd_scores


# if __name__ == '__main__':
#     amp_folder = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/DNAquaIMG/data/output_data/diat-barcode-test-amplicon/amplicon/filtered"
#     insert_folder = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/DNAquaIMG/data/output_data/diat-barcode-test-amplicon/insert/filtered"
#     primer_table = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/DNAquaIMG/data/input_data/diat-barcode-primers.tsv"


#     cls_instance = Binding()
#     mm_dict = cls_instance.get_number_of_primer_pbs_mismatches(amp_folder, insert_folder, primer_table)
#     mm_max = cls_instance.get_max_mismatches_per_taxon(mm_dict)
#     mm_max
