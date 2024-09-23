"""
This module contains the methods needed to perform the in-silico amplification analysis.
"""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from src.reference_database.db_curation import CrabsScriptGenerator


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(
        self, data, extract_primers: bool = True, primer_table: pd.DataFrame = None
    ):
        self.data = data
        self.output_dir = "../data/output_data"
        self.extract_primers = extract_primers
        self.primers = (
            defaultdict(lambda: [None, None, None]) if extract_primers else None
        )
        self.primer_table = primer_table
        self.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fw_seq",
            "rev_seq",
            "average_size",
            "min_lenght",
            "max_lenght",
        ]  # todo: Define mandatory columns
        self.crabs_script_generator = CrabsScriptGenerator()

    def _check_if_cutadapt_installed(self):
        """
        Function to check if Cutadapt is installed.
        """
        print("mozaiko INFO: Checking if cutadapt is installed...")
        try:
            subprocess.run(
                ["cutadapt", "--version"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

            print("mozaiko INFO: cutadapt is installed.")

        except FileNotFoundError:
            print(
                "mozaiko INFO: Cutadapt is not installed. Please install Cutadapt before running \
                    this script."
            )
            print(
                "Cutadapt can be found at "
                + "https://cutadapt.readthedocs.io/en/stable/installation.html"
            )
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

    def _validate_primer_table(self, primer_table):
        """
        Function to validate the primer table.
        """

        if not os.path.exists(primer_table):
            print("mozaiko INFO: The primer table does not exist. Exiting...")
            sys.exit(1)

        _, file_extension = os.path.splitext(primer_table)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "tsv":
            print("mozaiko INFO: The primer table must be a TSV file. Exiting...")
            sys.exit(1)

    def read_primer_tables(self):
        """
        Method to read and extract the required properties from the primer table.
        """

        print(
            f"mozaiko INFO: To continue the analysis, a set of primers is needed. "
            f"This information should be uploaded as a TSV table and it should contain the "
            f"following fields: {self.primer_table_columns}"
        )

        primer_table = input("Please enter the path to the primer table: ")

        self._validate_primer_table(primer_table)

        primer_table = pd.read_csv(primer_table, sep="\t", header=0)

        primer_table_fields = primer_table.columns.tolist()

        if primer_table_fields != self.primer_table_columns:
            print(
                f"mozaiko INFO: The primer table must contain the following fields: "
                f"{self.primer_table_columns}"
            )
            sys.exit(1)

        for index, row in primer_table.iterrows():

            foward_primer = row["fw_seq"].replace("I", "N")
            reverse_primer = row["rev_seq"].replace("I", "N")

            correct_reverse_primer = str(Seq(reverse_primer).reverse_complement())

            forward_primer_lenght = len(foward_primer)
            correct_reverse_primer_lenght = len(correct_reverse_primer)

            overlap = min(forward_primer_lenght, correct_reverse_primer_lenght)
            adapter = foward_primer + "..." + correct_reverse_primer

            # function to set the maximum lenght of overlap allowed
            # This implementation was decided internally
            max_len_overlap_formula = (
                600 - correct_reverse_primer_lenght - forward_primer_lenght - 60
            )

            primer_table.at[index, "overlap"] = overlap
            primer_table.at[index, "adapter"] = adapter
            primer_table.at[index, "max_overlap"] = max_len_overlap_formula

        self.primer_table = primer_table

    def _detect_fwd_rev_primer_len(self, sequence):
        """
        Calculate the length of the forward and reverse primer.

        Parameters:
        sequence (str): DNA sequence.

        Returns:
        int: Length of the forward and reverse primer.
        """

        fwd_len = len(sequence) // 2
        rev_len = len(sequence) - fwd_len

        return fwd_len, rev_len

    def _calculate_ambiguous_percentage(self, sequence):
        """
        Calculate the percentage of ambiguous bases in a DNA sequence.

        Parameters:
        sequence (str): The DNA sequence.

        Returns:
        float: The percentage of ambiguous bases in the sequence.
        """
        ambiguous_bases = set("RYWSMKHBVDN")

        return sum(base in ambiguous_bases for base in sequence) / len(sequence)

    def write_filtered_sequence(self, output_handle, record):
        """
        Write a filtered sequence to the output file.
        """
        sequence = str(record.seq)
        output_handle.write(f">{record.description}\n{sequence}\n")

    def _filter_sequences_by_prcnt_ambiguous_bases(
        self, input_file, out_file, max_ambiguous_percentage=0.05
    ): # todo: check usage
        """
        Filter DNA sequences based on the maximum allowed percentage of ambiguous bases.

        Parameters:
        - input_file (str): Path to the input file containing DNA sequences in FASTA format.
        - out_file (str): Path to the output file to write the filtered sequences.

        Returns:
        - dict or None: Dictionary of extracted primers if `extract_primers` is True, otherwise
        None.
        """

        with open(out_file, "w", encoding="UTF-8") as output_handle:

            for record in SeqIO.parse(input_file, "fasta"):
                sequence = str(record.seq)
                ambiguous_percentage = self._calculate_ambiguous_percentage(sequence)

                if ambiguous_percentage <= max_ambiguous_percentage:

                    if self.extract_primers:
                        fwd_len, rev_len = self._detect_fwd_rev_primer_len(sequence)
                        accession_number = record.id
                        # stores forward primer
                        self.primers[accession_number][0] = sequence[:fwd_len]
                        # stores reverse primer
                        self.primers[accession_number][1] = sequence[-rev_len:]
                        # stores sequence lenght (total lenght minus lenght of both primers)
                        self.primers[accession_number][2] = (
                            len(sequence) - fwd_len - rev_len
                        )
                    self.write_filtered_sequence(output_handle, record)

        return self.primers if self.extract_primers else None

    def _validate_fasta(self):
        """
        This method validates the input fasta.
        """

        if not os.path.exists(self.data):
            print("mozaiko INFO: The input file does not exist. Exiting...")
            sys.exit(1)

        print("mozaiko INFO: Input FASTA exists. Validating file extension...")

        _, file_extension = os.path.splitext(self.data)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "fasta":
            print("mozaiko INFO: Input file must be a FASTA file. Exiting...")
            sys.exit(1)

    def run_in_silico_analysis(self):
        """
        This methods writes and runs the cutadapt command, responsible for the in-silico
        amplification.
        """

        self._check_if_cutadapt_installed()
        self.crabs_script_generator.check_if_crabs_installed()

        self._validate_fasta()

        self.read_primer_tables()
        print("mozaiko INFO: All set. Running in-silico amplification...")

        run_name = Path(input("Please enter a name for the folder where the analysis output will \
                              be stored: "))

        input_fasta = self.data

        for _, row in self.primer_table.iterrows():
            self.process_comands(row, run_name, input_fasta)

        print("mozaiko INFO: In-silico amplification analysis completed.")

    def process_comands(self, row, run_name, input_fasta):
        barcode_region = row["barcode_region"]
        assay_name = row["assay_name"]
        five_prime_adapter = row["adapter"]
        max_length = int(row["max_overlap"])
        overlap = int(row["overlap"])
        forward_primer = row["fw_seq"]
        reverse_primer = row["rev_seq"]

        output_dirs = {
            "successful_amplification": self.output_dir
            / run_name
            / "successful_amplification",
            "pbr_no_amplification": self.output_dir / run_name / "pbr_no_amplification",
            "inserts_pbr": self.output_dir / run_name / "inserts_pbr",
            "pga": self.output_dir / run_name / "pga",
        }

        for dir_path in output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        self._run_cutadapt_command(
            "successful_amplification",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["successful_amplification"],
        )

        self._run_cutadapt_command(
            "pbr_no_amplification",
            five_prime_adapter,
            input_fasta,  # todo: ref daxtabase
            overlap,
            None,
            barcode_region,
            assay_name,
            output_dirs["pbr_no_amplification"],
            error_rate=5,
        )

        self._run_cutadapt_command(
            "inserts_pbr",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["inserts_pbr"],
        )

        self._run_pga_command(
            input_fasta,
            forward_primer,
            reverse_primer,
            barcode_region,
            assay_name,
            output_dirs["pga"],
            output_dirs["inserts_pbr"],
        )

    def _run_cutadapt_command(
        self,
        command_type,
        adapter,
        input_file,
        overlap,
        max_length,
        barcode_region,
        assay_name,
        output_dir,
        error_rate=3,
    ):
        output_file = output_dir / f"{barcode_region}_{assay_name}.txt"

        base_command = [
            "cutadapt",
            "-g", adapter,
            "--output", str(output_file),
            str(input_file),
            "--no-indels",
            "-e", str(error_rate),
            "--overlap", str(overlap),
            "--revcomp",
            "--quiet"
        ]

        if command_type == "successful_amplification":
            additional_args = ["--action", "retain", "--discard-untrimmed", "--maximum-length", str(max_length)]
        elif command_type == "pbr_no_amplification":
            additional_args = ["--action", "trim", "--discard-untrimmed"]
        elif command_type == "inserts_pbr":
            additional_args = ["--action", "trim", "--discard-untrimmed", "--maximum-length", str(max_length)]

        full_command = base_command + additional_args

        #print(f"mozaiko INFO: Running cutadapt command as: {' '.join(full_command)}")
        # print(f"mozaiko INFO: Input file: {input_file}")
        # print(f"mozaiko INFO: Output file: {output_file}")

        try:
            result = subprocess.run(full_command, check=True, capture_output=True, text=True)
            #print(f"mozaiko INFO: Cutadapt stdout:\n{result.stdout}")
            #print(f"mozaiko INFO: Cutadapt stderr:\n{result.stderr}")

            # if output_file.stat().st_size == 0:
            #     print(f"mozaiko WARNING: Output file is empty: {output_file}")
            # else:
            #     print(f"mozaiko INFO: Output file size: {output_file.stat().st_size} bytes")
        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: cutadapt {command_type} command failed: {e}")
            print(f"mozaiko ERROR: Cutadapt stdout:\n{e.stdout}")
            print(f"mozaiko ERROR: Cutadapt stderr:\n{e.stderr}")
            raise  # Re-raise the exception instead of calling sys.exit()
        except FileNotFoundError:
            print("mozaiko ERROR: cutadapt command not found. Please ensure it's installed and in your PATH.")
            raise

    def _run_pga_command(
        self,
        input_file,
        forward_primer,
        reverse_primer,
        barcode_region,
        assay_name,
        output_dir,
        database_dir,
    ):
        output_file = output_dir / f"{barcode_region}_{assay_name}.txt"
        pga_database = database_dir / f"{barcode_region}_{assay_name}.txt"

        minimum_percentage_identity = 0.95
        minimum_alignment_coverage = 0.99

        pga_command = [
            "crabs",
            "pga",
            "--input", str(input_file),
            "--output", str(output_file),
            "--database", str(pga_database),
            "--fwd", forward_primer,
            "--rev", reverse_primer,
            "--speed", "slow",
            "--percid", str(minimum_percentage_identity),
            "--coverage", str(minimum_alignment_coverage),
            "--filter_method", "strict"
        ]

        try:
            # print(f"mozaiko INFO: Running cutadapt command as '{pga_command}'")
            subprocess.run(pga_command, check=True)

        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: CRABS PGA command failed: {e}")
            sys.exit(1)
