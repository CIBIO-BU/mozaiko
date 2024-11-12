import os
import sys

import pandas as pd

from src.reference_database.sequence_import import CustomFastaImport
from src.in_silico_analysis.amplification import InSilicoAmplification


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
    def __init__(self, number_of_mismatches: int = None):
        self.amplification_instance = InSilicoAmplification()
        if number_of_mismatches is None:
            number_of_mismatches = (
                self.amplification_instance.get_number_of_mismatches()
            )
        self.number_of_mismatches = number_of_mismatches


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
