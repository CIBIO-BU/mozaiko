"""
This module contains the methods needed to perform the in-silico amplification analysis.
"""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from src.reference_database.db_curation import CrabsScriptGenerator


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(self, data, primer_table: pd.DataFrame = None):
        self.data = data
        self.output_dir = "../data/output_data"
        self.primer_table = primer_table
        self.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fw_seq",
            "rev_seq"
        ]
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

        primer_table = pd.read_csv(primer_table, sep="\t", header=0)

        primer_table_fields = primer_table.columns.tolist()

        required_fields = set(self.primer_table_columns)
        provided_fields = set(primer_table_fields)

        if not required_fields.issubset(provided_fields):
            missing_fields = required_fields - provided_fields
            print(
                f"mozaiko INFO: The primer table is missing the following required fields: "
                f"{', '.join(missing_fields)}"
            )
            print(f"Required fields are: {', '.join(self.primer_table_columns)}")
            sys.exit(1)

    def read_primer_tables(self, primer_table=None):
        """
        Method to read and extract the required properties from the primer table.
        """

        print(
            f"mozaiko INFO: To continue the analysis, a set of primers is needed. "
            f"This information should be uploaded as a TSV table and it should contain the "
            f"following fields: {self.primer_table_columns}"
        )

        if primer_table is None:
            primer_table = input("Please enter the path to the primer table: ")

        self._validate_primer_table(primer_table)

        primer_table = pd.read_csv(primer_table, sep="\t", header=0)

        for index, row in primer_table.iterrows():

            foward_primer = row["fw_seq"]
            reverse_primer = row["rev_seq"]

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

    def run_in_silico_analysis(self, primer_table=None):
        """
        This methods initiates the in-silico analysis. It does so by first veryfing if all required
        tools are installed in the machine. If installed, it requests the user to upload a table
        containing a list of primers to be evaluated and to provide a name for the folder where
        the results of the run will be stored. At last it processes the inputed primer table to
        iterate over each provided primer (row) and process the needed commands.
        """

        self._check_if_cutadapt_installed()
        self.crabs_script_generator.check_if_crabs_installed()

        self._validate_fasta()

        self.read_primer_tables(primer_table)
        print("mozaiko INFO: All set. Running in-silico amplification...")

        run_name = Path(
            input(
                "Please enter a name for the folder where the analysis output will \
                              be stored: "
            )
        )

        input_fasta = self.data

        for _, row in self.primer_table.iterrows():
            self.process_commands(row, run_name, input_fasta)

        print("mozaiko INFO: In-silico amplification analysis completed.")

    def process_commands(self, row, run_name, input_fasta):
        """
        This method creates variables from the user-inputted primer table to process commands for
        the in-silico ammplication.
        """
        barcode_region = row["barcode_region"]
        assay_name = row["assay_name"]
        five_prime_adapter = row["adapter"]
        max_length = int(row["max_overlap"])
        overlap = int(row["overlap"])
        forward_primer = row["fw_seq"]
        reverse_primer = row["rev_seq"]

        # Name output directories for each of the analysis steps
        output_dirs = {
            "amplicon": self.output_dir / run_name / "amplicon",
            "all_barcodes_w_pbr": self.output_dir / run_name / "all_barcodes_w_pbr",
            "insert": self.output_dir / run_name / "insert",
            "pga": self.output_dir / run_name / "pga",
        }

        for dir_path in output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        # "amplicon" comand makes use of --action=retain to trim the amplicon but not remove the
        # primer binding sites (sequences before and after the PBS are removed)
        self._run_cutadapt_command(
            "amplicon",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["amplicon"],
        )

        self._run_cutadapt_command(
            "all_barcodes_w_pbr",
            five_prime_adapter,
            input_fasta,
            overlap,
            None,
            barcode_region,
            assay_name,
            output_dirs["all_barcodes_w_pbr"],
            error_rate=5,
        )

        # "insert" makes use of --action=trim to remove the primer binding site (and the sequence
        # before or after it)
        self._run_cutadapt_command(
            "insert",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["insert"],
        )

        self._run_pga_command(
            input_fasta,
            forward_primer,
            reverse_primer,
            barcode_region,
            assay_name,
            output_dirs["pga"],
            output_dirs["all_barcodes_w_pbr"],
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
        """
        This method designs the commands to run with cutadapt (https://github.com/marcelm/cutadapt).
        It first defines the output file for each barcode regions and assay name. Then defines a
        base command, following a list of alternative commands for each of the analysis purposes.
        """
        output_file = output_dir / f"{barcode_region}_{assay_name}.txt"

        base_command = [
            "cutadapt",
            "-g",
            adapter,
            "--output",
            str(output_file),
            str(input_file),
            "--no-indels",
            "-e",
            str(error_rate),
            "--overlap",
            str(overlap),
            "--revcomp",
            "--quiet",
        ]

        if command_type == "amplicon":
            additional_args = [
                "--action",
                "retain",
                "--discard-untrimmed",
                "--maximum-length",
                str(max_length),
            ]
        elif command_type == "all_barcodes_w_pbr":
            additional_args = ["--action", "trim", "--discard-untrimmed"]
        elif command_type == "insert":
            additional_args = [
                "--action",
                "trim",
                "--discard-untrimmed",
                "--maximum-length",
                str(max_length),
            ]

        full_command = base_command + additional_args

        # print(f"mozaiko INFO: Running cutadapt command as: {' '.join(full_command)}")
        # print(f"mozaiko INFO: Input file: {input_file}")
        # print(f"mozaiko INFO: Output file: {output_file}")

        try:
            result = subprocess.run(
                full_command, check=True, capture_output=True, text=True
            )
            # print(f"mozaiko INFO: Cutadapt stdout:\n{result.stdout}")
            # print(f"mozaiko INFO: Cutadapt stderr:\n{result.stderr}")

            if output_file.stat().st_size == 0:
                print(f"mozaiko WARNING: Output file is empty: {output_file}")
            # else:
            # print(f"mozaiko INFO: Output file size: {output_file.stat().st_size} bytes")
        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: cutadapt {command_type} command failed: {e}")
            # print(f"mozaiko ERROR: Cutadapt stdout:\n{e.stdout}")
            # print(f"mozaiko ERROR: Cutadapt stderr:\n{e.stderr}")
            raise  # Re-raise the exception instead of calling sys.exit()
        except FileNotFoundError:
            print(
                "mozaiko ERROR: cutadapt command not found. Please ensure it's installed and in your PATH."
            )
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
        """
        This method designs the command to run Pairwise Global Alignment (PGA) with CRABS
        (https://github.com/gjeunen/reference_database_creator).
        It first defines the output file for each barcode regions and assay name; and the
        input file to be used as a database to be searched against the reference sequences.
        Finally, it states the base command to be ran.
        """
        output_file = output_dir / f"{barcode_region}_{assay_name}.txt"
        pga_database = database_dir / f"{barcode_region}_{assay_name}.txt"

        minimum_percentage_identity = 0.95
        minimum_alignment_coverage = 0.99

        pga_command = [
            "crabs",
            "pga",
            "--input",
            str(input_file),
            "--output",
            str(output_file),
            "--database",
            str(pga_database),
            "--fwd",
            forward_primer,
            "--rev",
            reverse_primer,
            "--speed",
            "slow",
            "--percid",
            str(minimum_percentage_identity),
            "--coverage",
            str(minimum_alignment_coverage),
            "--filter_method",
            "strict",
        ]

        try:
            # print(f"mozaiko INFO: Running cutadapt command as '{pga_command}'")
            subprocess.run(pga_command, check=True)

        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: CRABS PGA command failed: {e}")
            sys.exit(1)
