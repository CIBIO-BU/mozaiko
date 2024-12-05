"""
This module contains the methods needed to perform the in-silico amplification analysis.

PBS: Primer Binding Site
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pandas import DataFrame

from src.marker_scoring.scoring_utils import filter_sequences_by_ambiguity
from src.reference_database.db_curation import CrabsScriptGenerator
from src.reference_database.sequence_import import CustomFastaImport


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(
        self,
        data: Optional[Path] = None,
        primer_table: Optional[DataFrame] = None,
        run_name: Optional[str] = None,
        number_of_mismatches: int = 3,
    ):
        self.data: Optional[Path] = data
        self.base_output_dir = Path("../data/output_data")
        self.primer_table = primer_table
        self.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fwd_seq",
            "rev_seq",
        ]
        self.crabs_script_generator = CrabsScriptGenerator()
        self.run_name: Optional[str] = run_name
        self.number_of_mismatches = number_of_mismatches

    def get_number_of_mismatches(self):
        """
        Gets the defined number of mismatches to use as variables in other external classes.
        """
        return self.number_of_mismatches

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

        self.run_dir = self.base_output_dir / run_name

        output_dirs = {
            "amplicon": self.run_dir / "amplicon",
            "insert": self.run_dir / "insert",
            "all_complete_pbs": self.run_dir / "all_complete_pbs",
            "incomplete_pbs": self.run_dir / "incomplete_pbs",
        }

        for dir_path in output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        self.output_dirs = output_dirs

        return output_dirs

    def _ensure_output_dirs(self) -> Dict[str, Path]:
        if self.output_dirs is None:
            raise ValueError(
                "Output directories not set up. Call run_in_silico_analysis first."
            )
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

    def validate_primer_table(self, primer_table):
        """
        Function to validate the primer table.
        """

        if not os.path.exists(primer_table):
            raise ValueError(
                "mozaiko INFO: The primer table does not exist. Exiting..."
            )

        _, file_extension = os.path.splitext(primer_table)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "tsv":
            raise ValueError(
                "mozaiko INFO: The primer table must be a TSV file. Exiting..."
            )

        primer_table = pd.read_csv(primer_table, sep="\t", header=0, dtype=str)

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

        listed_required_fields = list(required_fields)
        if primer_table[listed_required_fields].isnull().any().any():
            print(
                "mozaiko INFO: The primer table contains rows with missing values. Please fill in all required columns."
            )
            sys.exit(1)

        if not primer_table["assay_name"].is_unique:
            print(
                "mozaiko INFO: The 'assay_name' column contains duplicate values. Each assay name must be unique."
            )
            sys.exit(1)

    def read_primer_tables(
        self, primer_table=None, max_len_according_to_ilumina: bool = True
    ):
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

        self.validate_primer_table(primer_table)

        primer_table = pd.read_csv(primer_table, sep="\t", header=0)

        for index, row in primer_table.iterrows():

            foward_primer = row["fwd_seq"]
            reverse_primer = row["rev_seq"]

            correct_reverse_primer = str(Seq(reverse_primer).reverse_complement())

            forward_primer_length = len(foward_primer)
            correct_reverse_primer_length = len(correct_reverse_primer)

            min_fwd_overlap = str(forward_primer_length)
            min_rev_overlap = str(correct_reverse_primer_length)

            adapter = (
                foward_primer
                + ";"
                + "min_overlap="
                + min_fwd_overlap
                + "..."
                + correct_reverse_primer
                + ";"
                + "min_overlap="
                + min_rev_overlap
            )

            if max_len_according_to_ilumina is True:
                max_len_formula = (
                    600 - correct_reverse_primer_length - forward_primer_length
                )
                primer_table.at[index, "max_read_length"] = max_len_formula

            else:
                if (
                    "min_read_length" not in self.primer_table_columns
                    and "max_read_length" not in self.primer_table_columns
                ):
                    raise ValueError(
                        "mozaiko ERROR: When max_len_according_to_ilumina is False, the primer table must contain 'min_read_length' and 'max_read_length' columns."
                    )

            # primer_table.at[index, "min_fwd_overlap"] = min_fwd_overlap
            # primer_table.at[index, "min_rev_overlap"] = min_rev_overlap
            primer_table.at[index, "adapter"] = adapter
            primer_table.at[index, "correct_reverse_primer"] = correct_reverse_primer

        self.primer_table = primer_table

    def _validate_fasta(self):
        """
        This method validates the input fasta.
        """

        if not os.path.exists(self.data):
            print("mozaiko INFO: The input file does not exist. Exiting...")
            sys.exit(1)

        _, file_extension = os.path.splitext(self.data)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "fasta":
            print("mozaiko INFO: Input file must be a FASTA file. Exiting...")
            sys.exit(1)

    def count_sequences(self, fasta_file):
        """
        Count the number of sequences in a FASTA file.

        Parameter:
        - fasta_file: Path to the FASTA file

        Returns:
        - int: Number of sequences
        """
        return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    def remove_intersection_sequences(self, input_path, filter_path):
        """
        This method removes sequences that are present in both files. As such, it filters the
        sequences present in the filter file (filter_path) from the other (input_path).
        """
        input_path = Path(input_path)
        filter_path = Path(filter_path)

        is_input_dir = input_path.is_dir()

        # Create filtered output directory
        filtered_output_dir = (input_path.parent if not is_input_dir else input_path) / "filtered_intersection"
        filtered_output_dir.mkdir(parents=True, exist_ok=True)

        if input_path.is_file():
            input_files = [input_path]
        elif input_path.is_dir():
            input_files = list(input_path.glob("*.fasta"))
        else:
            raise ValueError(f"mozaiko ERROR: Input path {input_path} does not exist")

        if not input_files:
            raise ValueError(f"mozaiko ERROR: No FASTA files found in {input_path}")

        if filter_path.is_file():
            filter_files = [filter_path]
        elif filter_path.is_dir():
            filter_files = list(filter_path.glob("*.fasta"))
        else:
            raise ValueError(f"mozaiko ERROR: Input path {filter_path} does not exist")

        if not filter_files:
            raise ValueError(f"mozaiko ERROR: No FASTA files found in {filter_path}")

        filter_primer_mapping = {f.stem: f for f in filter_files}

        results = {}

        for input_file in input_files:
            if input_file.stem not in filter_primer_mapping:
                print(
                    f"mozaiko INFO: No incomplete PBSs found for {input_file}, all PBS are complete."
                )
                continue

            matching_primer = filter_primer_mapping[input_file.stem]

            # all_PBS_count = self.count_sequences(input_file)
            # incomplete_psb_count = self.count_sequences(matching_primer)

            output_file = filtered_output_dir / input_file.name

            ids_to_filter = set()
            for record in SeqIO.parse(matching_primer, "fasta"):
                # seq_string_filtered = str(record.seq).upper()
                id_filtered = record.id
                ids_to_filter.add(id_filtered)

            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
                retained_sequences = set()
                for record in SeqIO.parse(input_file, "fasta"):
                    seq_string_retained = str(record.seq).upper()
                    record_retained = record.description
                    record_id = record.id
                    if record_id not in ids_to_filter:
                        retained_sequences.add((seq_string_retained, record_retained))

                with open(output_file, "w") as output_handle:
                    for seq_string_retained, record_retained in retained_sequences:
                        output_handle.write(
                            f">{record_retained}\n{seq_string_retained}\n"
                        )

                results[input_file] = {
                    "retained_sequences": len(retained_sequences),
                    "filter_file_used": matching_primer,
                    "output_file": output_file,
                }

            print(
                f"    For {input_file.stem}: {results[input_file]['retained_sequences']} sequences were retained."
            )

    def add_taxonomy_to_pga_outputs(self, input_folders):
        """
        This method adds taxonomy information to FASTA file headers where it's missing (PGA outputs).

        Parameters:
        - input_folders (list):
        List of folder paths containing FASTA files to process
        """
        # Load mapping between taxonomy and seq-id
        self.custom_fasta_import = CustomFastaImport(self.data)
        self.custom_fasta_import.read_fasta(self.data)
        seq_id_taxonomy_dict = self.custom_fasta_import.get_mapping_between_seq_id_taxonomy()

        for folder_path in input_folders:
            folder_path = Path(folder_path)
            fasta_files = list(folder_path.glob("*.fasta"))

            for fasta_file in fasta_files:
                records = list(SeqIO.parse(fasta_file, "fasta"))

                modified = False
                for record in records:
                    if "|" not in record.description:
                        if record.id in seq_id_taxonomy_dict:
                            new_description = f"{record.description} | {seq_id_taxonomy_dict[record.id]}"
                            record.description = new_description
                            record.id = new_description.split()[0]
                            modified = True

                if modified:
                    with open(fasta_file, "w") as output_handle:
                        for record in records:
                            output_handle.write(
                                f">{record.description}\n{str(record.seq)}\n"
                            )

    def run_in_silico_analysis(
        self, primer_table=None, max_len_according_to_ilumina: bool = True
    ):
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

        if max_len_according_to_ilumina is True:
            self.read_primer_tables(primer_table, max_len_according_to_ilumina=True)
        if max_len_according_to_ilumina is False:
            self.read_primer_tables(primer_table, max_len_according_to_ilumina=False)

        if not self.run_name:
            self.run_name = input(
                "Please enter a name for the folder where the analysis output will be stored: "
            ).strip()

            if not self.run_name:
                raise ValueError("Run name cannot be empty")

        self._setup_output_directories(self.run_name)

        if self.primer_table is None:
            raise ValueError("Primer table not initialized")

        print("mozaiko INFO: All set. Running in-silico amplification...")

        for index, row in self.primer_table.iterrows():
            if self.data:
                self.process_commands(row, self.data)
            else:
                raise ValueError("mozaiko ERROR: No input data was found.")
            print(
                f"   --------   {index + 1}/{len(self.primer_table)} processed   --------   "
            )

        pga_directories = [
            self.output_dirs["all_complete_pbs"],
            self.output_dirs["incomplete_pbs"],
        ]
        self.add_taxonomy_to_pga_outputs(pga_directories)

        directories_to_filter = ["amplicon", "all_complete_pbs", "incomplete_pbs"]
        for dir_name in directories_to_filter:
            try:
                input_path = self.output_dirs[dir_name]
                filter_sequences_by_ambiguity(
                    input_path=input_path,
                )
            except Exception as e:
                print(f"mozaiko ERROR: Error filtering {dir_name} directory: {str(e)}")

        print("mozaiko INFO: Number of inserts that were amplified successfully...")
        filter_inserts_path = self.output_dirs["insert"] / "filtered"
        insert_files = list(filter_inserts_path.glob("*.fasta"))
        for file in insert_files:
            number_of_sequences = self.count_sequences(file)
            print(f"    For {file.stem}, {number_of_sequences} were retained.")

        # Make a copy of the inserts file to avoid alterations in downstream tasks (intersections)
        all_inserts_path = self.run_dir / "all_inserts"
        shutil.copytree(filter_inserts_path, all_inserts_path)

        print(
            "mozaiko INFO: Number of inserts with complete PBS that were not amplified..."
        )
        self.remove_intersection_sequences(
            self.output_dirs["all_complete_pbs"] / "filtered",
            self.output_dirs["incomplete_pbs"] / "filtered",
        )

        print("mozaiko INFO: Number of inserts with incomplete PBS...")
        self.remove_intersection_sequences(
            self.output_dirs["all_complete_pbs"] / "filtered",
            self.output_dirs["insert"] / "filtered",
        )

        print("mozaiko INFO: In-silico amplification analysis completed.")

    def process_commands(self, row: dict, input_fasta: Path):
        """
        This method creates variables from the user-inputted primer table to process commands for
        the in-silico ammplication.
        """
        output_dirs = self._ensure_output_dirs()

        barcode_region = row["barcode_region"]
        assay_name = row["assay_name"]
        adapter = row["adapter"]
        max_length = int(row["max_read_length"])
        if "min_read_length" in row.keys():
            min_length = int(row["min_read_length"])
        # min_fwd_overlap = int(row["min_fwd_overlap"])
        # min_rev_overlap = int(row["min_rec_overlap"])
        forward_primer = row["fwd_seq"]
        reverse_primer = row["correct_reverse_primer"]

        # "amplicon" comand makes use of --action=retain to trim the amplicon but not remove the
        # PBS (sequences before and after the PBS are removed)
        self.run_cutadapt_command(
            "amplicon",
            adapter,
            input_fasta,
            min_length,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["amplicon"],
        )

        # "insert" makes use of --action=trim
        self.run_cutadapt_command(
            "insert",
            adapter,
            input_fasta,
            min_length,
            max_length,
            barcode_region,
            assay_name,
            output_dirs["insert"],
        )

        # Filter sequences per ambiguous bases prior to their use as reference database in PGA
        try:
            filter_sequences_by_ambiguity(
                input_path=output_dirs["insert"],
                output_dir=output_dirs["insert"] / "filtered",
                max_ambiguous_percentage=0.05,
            )
        except Exception as e:
            print(f"Error filtering insert directory: {str(e)}")

        # Retrieve all inserts with a PBS
        self.run_pga_command(
            input_fasta,
            forward_primer,
            reverse_primer,
            barcode_region,
            assay_name,
            output_dirs["all_complete_pbs"],
            output_dirs["insert"] / "filtered",
            "relaxed",
        )

        # Retrieve all inserts with incomplete PBS
        self.run_pga_command(
            input_fasta,
            forward_primer,
            reverse_primer,
            barcode_region,
            assay_name,
            output_dirs["incomplete_pbs"],
            output_dirs["insert"] / "filtered",
            "strict",
        )

        print(f"mozaiko INFO: Completed analysis for {assay_name}.")

    def run_cutadapt_command(
        self,
        command_type,
        adapter,
        input_file,
        min_length,
        max_length,
        barcode_region,
        assay_name,
        output_dir,
        number_of_mismatches=3,
    ):
        """
        This method designs the commands to run with cutadapt (https://github.com/marcelm/cutadapt).
        It first defines the output file for each barcode regions and assay name. Then defines a
        base command, following a list of alternative commands for each of the analysis purposes.
        """
        if self.number_of_mismatches != number_of_mismatches:
            number_of_mismatches = self.number_of_mismatches

        output_dir = Path(output_dir)
        output_file = output_dir / f"{barcode_region}_{assay_name}.fasta"

        # debug_dir = None
        # if command_type == 'all_barcodes_w_pbr':
        #     debug_dir = output_dir / "debug_matrices" / f"{barcode_region}_{assay_name}.fasta"
        #     debug_dir.mkdir(parents=True, exist_ok=True)

        base_command = [
            "cutadapt",
            "-g",
            adapter,
            "--output",
            str(output_file),
            str(input_file),
            "--no-indels",
            "-e",
            str(number_of_mismatches),
            "--revcomp",
            "--quiet",
            "--minimum-length",
            str(min_length),
            "--maximum-length",
            str(max_length),
            "--discard-untrimmed",
        ]

        if command_type == "amplicon":

            additional_args = ["--action", "retain"]
        elif command_type == "insert":
            additional_args = ["--action", "trim"]
        else:
            raise ValueError(f"Invalid command type: {command_type}")

        full_command = base_command + additional_args

        # print(f"mozaiko INFO: Running cutadapt command as: {' '.join(full_command)}")
        # print(f"mozaiko INFO: Input file: {input_file}")
        # print(f"mozaiko INFO: Output file: {output_file}")

        try:
            result = subprocess.run(
                full_command,
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
            )
            if output_file.stat().st_size == 0:
                print(
                    f"mozaiko WARNING: No {command_type}s retrieved for {output_file.stem}."
                )

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
        filter,
    ):
        """
        This method designs the command to run Pairwise Global Alignment (PGA) with CRABS
        (https://github.com/gjeunen/reference_database_creator).
        It first defines the output file for each barcode regions and assay name; and the
        input file to be used as a database to be searched against the reference sequences.
        Finally, it states the base command to be ran.
        """
        output_file = output_dir / f"{barcode_region}_{assay_name}.fasta"
        pga_database = database_dir / f"{barcode_region}_{assay_name}.fasta"

        minimum_percentage_identity = 0.75  # Decimal Percentage [0 - 1.0]
        minimum_alignment_coverage = 99  # Whole Percentage [0 - 100]

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
            str(filter),
        ]

        try:
            # print(f"mozaiko INFO: Running CRABS command as '{' '.join(pga_command)}'")
            subprocess.run(
                pga_command,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )

        except subprocess.CalledProcessError as e:
            print(f"mozaiko ERROR: CRABS PGA command failed: {e}")
            sys.exit(1)


# if __name__ == "__main__":
#     data = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/DNAquaIMG/data/input_data/diat-barcode-taxa.fasta"
#     primer_table = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/DNAquaIMG/data/input_data/diat-barcode-primers.tsv"
#     run_name = "diat-barcode-test-checkalters"
#     cutadapt = InSilicoAmplification(data, run_name=run_name)
#     cutadapt.run_in_silico_analysis(primer_table=primer_table)
