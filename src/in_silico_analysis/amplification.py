"""
This module contains the methods needed to perform the in-silico amplification analysis.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd
from Bio.Seq import Seq
from pandas import DataFrame

from ..marker_scoring.scoring_utils import filter_sequences_by_ambiguity
from ..reference_database.db_curation import CrabsScriptGenerator


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(self, data: Union[str, Path], primer_table: Optional[DataFrame] = None, run_name: Optional[str] = None):
        self.data = data
        self.base_output_dir = Path("../data/output_data")
        self.primer_table = primer_table
        self.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fw_seq",
            "rev_seq",
        ]
        self.crabs_script_generator = CrabsScriptGenerator()
        self.run_name: Optional[str] = run_name
        self.output_dirs: Optional[Dict[str, Path]] = None

    def _setup_output_directories(self, run_name: str) -> dict:
        """
        Sets up output directory structure based on run name.

        Parameters:
        - run_name: Name of the analysis run

        Returns:
        Dictionary containing Path objects for each output directory
        """
        if not run_name:
            raise ValueError("Run name must be provided")

        run_dir = self.base_output_dir / run_name

        output_dirs = {
            "amplicon": run_dir / "amplicon",
            "all_barcodes_w_pbr": run_dir / "all_barcodes_w_pbr",
            "insert": run_dir / "insert",
            "pga": run_dir / "pga",
        }

        for dir_path in output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        self.output_dirs = output_dirs

        return output_dirs

    def _ensure_output_dirs(self) -> Dict[str, Path]:
        if self.output_dirs is None:
            raise ValueError("Output directories not set up. Call run_in_silico_analysis first.")
        return self.output_dirs

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

        if not self.run_name:
            self.run_name = input(
                "Please enter a name for the folder where the analysis output will be stored: "
            ).strip()

            if not self.run_name:
                raise ValueError("Run name cannot be empty")

        self._setup_output_directories(self.run_name)

        if self.primer_table is None:
            raise ValueError("Primer table not initialized")

        for _, row in self.primer_table.iterrows():
            self.process_commands(row, self.data)

        directories_to_filter = ["amplicon", "all_barcodes_w_pbr", "pga"]
        for dir_name in directories_to_filter:
            try:
                # print(f"Filtering sequences with ambiguous bases in {dir_name} directory...")
                filter_sequences_by_ambiguity(
                    input_path=self.output_dirs,
                )
                # print(f"Completed filtering {dir_name} sequences")
            except Exception as e:
                print(f"Error filtering {dir_name} directory: {str(e)}")

        print("mozaiko INFO: In-silico amplification analysis completed.")

    def process_commands(self, row: dict, input_fasta: Path):
        """
        This method creates variables from the user-inputted primer table to process commands for
        the in-silico ammplication.
        """
        output_dirs = self._ensure_output_dirs()

        barcode_region = row["barcode_region"]
        assay_name = row["assay_name"]
        five_prime_adapter = row["adapter"]
        max_length = int(row["max_overlap"])
        overlap = int(row["overlap"])
        forward_primer = row["fw_seq"]
        reverse_primer = row["rev_seq"]

        # "amplicon" comand makes use of --action=retain to trim the amplicon but not remove the
        # primer binding sites (sequences before and after the PBS are removed)
        self.run_cutadapt_command(
            "amplicon",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["amplicon"],
        )

        # "insert" makes use of --action=trim to remove the primer binding site (and the sequence
        # before or after it)
        self.run_cutadapt_command(
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

        # Filter sequences per ambiguous bases prior to their use as reference database in PGA
        try:
            filter_sequences_by_ambiguity(
                input_path=output_dirs["all_barcodes_w_pbr"],
                output_dir=output_dirs["all_barcodes_w_pbr"] / "filtered",
                max_ambiguous_percentage=0.05,
            )
        except Exception as e:
            print(f"Error filtering insert directory: {str(e)}")

        # "insert" makes use of --action=trim
        self.run_cutadapt_command(
            "insert",
            five_prime_adapter,
            input_fasta,
            overlap,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["insert"],
        )

        self.run_pga_command(
            input_fasta,
            forward_primer,
            reverse_primer,
            barcode_region,
            assay_name,
            output_dirs["pga"],
            output_dirs["all_barcodes_w_pbr"] / "filtered",
        )

    def run_cutadapt_command(
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
        else:
            raise ValueError(f"Invalid command type: {command_type}")

        full_command = base_command + additional_args

        # print(f"mozaiko INFO: Running cutadapt command as: {' '.join(full_command)}")
        # print(f"mozaiko INFO: Input file: {input_file}")
        # print(f"mozaiko INFO: Output file: {output_file}")

        try:
            subprocess.run(
                full_command, check=True, capture_output=True, text=True
            )
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
                "mozaiko ERROR: cutadapt command not found. Please ensure it's installed and in \
                    your PATH."
            )
            raise

    def run_pga_command(
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
            print(f"mozaiko INFO: Running cutadapt command as '{pga_command}'")
            subprocess.run(pga_command, check=True)

        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: CRABS PGA command failed: {e}")
            sys.exit(1)
