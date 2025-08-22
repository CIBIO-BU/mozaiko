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

from src.marker_scoring.scoring_utils import *
from src.reference_database.db_curation import CrabsScriptGenerator
from src.reference_database.sequence_import import CustomFastaImport


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(
        self,
        database_fasta_file: Optional[Path] = None,
        primer_table: Optional[DataFrame] = None,
        run_name: Optional[str] = None,
        number_of_mismatches: int = 3,
    ):
        self.database_fasta_file: Optional[Path] = database_fasta_file
        self.base_output_dir = Path("./data/output_data")
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

        if not os.path.exists(self.database_fasta_file):
            print("mozaiko INFO: The input file does not exist. Exiting...")
            sys.exit(1)

        _, file_extension = os.path.splitext(self.database_fasta_file)

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

    def intersect_PGA_relaxed_and_strict(self, input_path, filter_path, output_dir_name="input_B"):
        """
        Currently solely used to create Input B. Where sequences from PGA Strict (Incomplete PBS)
        are filtered out of the sequences from PGA Relaxed (Complete PBS) to create "Inserts with
        complete PBS that were not in-silico amplified (> 3 mismatches)".

        Filter sequences present in the filter file out of the input file(s).
        Creates a directory with the specified name in self.run_dir and saves filtered sequences there.
        Filtering is based on sequence IDs found in filter files.

        Args:
            input_path: Path to input FASTA file or directory of FASTA files
            filter_path: Path to filter FASTA file or directory of FASTA files
            output_dir_name: Name of the output directory (default: "filtered_intersection")

        Returns:
            dict: Results summary with stats for each processed file
        """
        from pathlib import Path
        from Bio import SeqIO

        input_path = Path(input_path)
        filter_path = Path(filter_path)

        # Create output directory in self.run_dir
        filtered_output_dir = self.run_dir / output_dir_name
        filtered_output_dir.mkdir(parents=True, exist_ok=True)

        # Find input and filter files
        input_files = [input_path] if input_path.is_file() else list(input_path.glob("*.fasta"))
        filter_files = [filter_path] if filter_path.is_file() else list(filter_path.glob("*.fasta"))

        # Check if files exist
        if not input_files:
            raise ValueError(f"mozaiko ERROR: No FASTA files found in {input_path}")
        if not filter_files:
            raise ValueError(f"mozaiko ERROR: No FASTA files found in {filter_path}")

        # Create mapping from primer name to filter file
        filter_primer_mapping = {f.stem: f for f in filter_files}
        results = {}

        # Process each input file
        for input_file in input_files:
            if input_file.stem not in filter_primer_mapping:
                print(f"mozaiko INFO: No incomplete PBSs found for {input_file}, all PBS are complete.")
                continue

            # Get matching filter file and set up output file
            matching_filter = filter_primer_mapping[input_file.stem]
            output_file = filtered_output_dir / input_file.name

            # Get IDs to filter out
            ids_to_filter = {record.id for record in SeqIO.parse(matching_filter, "fasta")}

            # Read and filter sequences
            retained_sequences = []
            with open(output_file, "w") as output_handle:
                for record in SeqIO.parse(input_file, "fasta"):
                    if record.id not in ids_to_filter:
                        seq = str(record.seq).upper()
                        output_handle.write(f">{record.description}\n{seq}\n")
                        retained_sequences.append(record.id)

            # Record results
            results[input_file] = {
                "retained_sequences": len(retained_sequences),
                "filter_file_used": matching_filter,
                "output_file": output_file,
            }

            print(f" For {input_file.stem}: {len(retained_sequences)} sequences were retained.")

        return results

    def _load_taxonomy_mapping(self):
        """
        Load and cache taxonomy mapping from database file.
        """
        try:
            with open(self.database_fasta_file) as f:
                mapping = {}
                for line in f:
                    if line.startswith(">"):
                        parts = line[1:].strip().split("|")
                        if len(parts) >= 11:
                            accession = parts[0]
                            taxonomy = "|".join(parts[1:11])
                            mapping[accession] = taxonomy
                return mapping
        except Exception as e:
            print(f"mozaiko ERROR: Failed to load taxonomy mapping: {e}")
            raise

    def _process_fasta_file(self, fasta_file: Path, taxonomy_dict: Dict[str, str]):
        """
        Process a single FASTA file with error handling.
        """
        try:
            records = []
            modified = False

            # Read all records at once
            for record in SeqIO.parse(str(fasta_file), "fasta"):
                accession = record.id.split("|")[0]
                if accession in taxonomy_dict:
                    new_description = f"{accession}|{taxonomy_dict[accession]}"
                    record.description = new_description
                    record.id = new_description.split()[0]
                    modified = True
                records.append(record)

            # Only write if modifications were made
            if modified:
                with open(fasta_file, "w") as output_handle:
                    for record in records:
                        output_handle.write(
                            f">{record.description}\n{str(record.seq)}\n"
                        )

        except Exception as e:
            print(f"mozaiko ERROR: Error processing {fasta_file}: {e}")

    def add_taxonomy_to_pga_outputs(self, input_folders: list[str]) -> None:
        """
        Add taxonomy information to FASTA file headers where missing.

        Args:
            input_folders: List of folder paths containing FASTA files
            taxa_column_start: Starting column index for taxonomy information
            taxa_column_end: Ending column index for taxonomy information

        Raises:
            ValueError: If input parameters are invalid
            FileNotFoundError: If input folders don't exist
        """
        if not input_folders:
            raise ValueError("No input folders provided")

        taxonomy_dict = self._load_taxonomy_mapping()

        for folder_path in input_folders:
            folder = Path(folder_path)
            if not folder.exists():
                print(f"mozaiko ERROR: Folder not found: {folder}")
                continue

            fasta_files = list(folder.glob("*.fasta"))

            for fasta_file in fasta_files:
                self._process_fasta_file(fasta_file, taxonomy_dict)

    def sanity_check_on_mismatches(self, output_dir, testing = False):
        from src.marker_scoring.scoring_utils import (
            calculate_iupac_mismatches
        )
        if not output_dir:
            raise ValueError("mozaiko ERROR: Output directory not specified for mismatch sanity check.")

        if not isinstance(output_dir, Path):
            output_dir = Path(output_dir)

        if not output_dir.exists():
            raise FileNotFoundError(f"mozaiko ERROR: Output directory does not exist: {output_dir}")

        amplicon_folder = output_dir / "amplicon"
        insert_folder = output_dir / "insert"

        from src.marker_scoring.metrics_system import (
            Binding
        )
        matching_files = Binding.parse_files_with_same_extension_in_folders(
            amplicon_folder, insert_folder
        )

        seq_ids_to_remove = {}
        mismatch_per_seq_id = {}

        if not matching_files:
            print("mozaiko ERROR: No matching primer files found between the insert and amplicon folders.")
            return None, None

        for _primer_ind, primer_row in self.primer_table.iterrows():
            barcode_region = primer_row["barcode_region"]
            assay_name = primer_row["assay_name"]
            pbs_filename = f"{barcode_region}_{assay_name}"
            primer_seq_fwd = primer_row["fwd_seq"]
            primer_seq_rev = primer_row["rev_seq"]
            rev_comp_primer_seq_rev = str(Seq(primer_seq_rev).reverse_complement())

            for amplicon_file, insert_file in matching_files:
                amplicon_filename = os.path.splitext(os.path.basename(amplicon_file))[0]

                if pbs_filename == amplicon_filename:
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


                        full_fwd_mismatches = calculate_iupac_mismatches(
                                    primer_seq_fwd, pbs_fwd_seq
                                )
                        full_rev_mismatches = calculate_iupac_mismatches(
                            rev_comp_primer_seq_rev, pbs_rev_seq
                        )

                        full_len_mismatch_sum = full_fwd_mismatches + full_rev_mismatches

                        if seq_id not in mismatch_per_seq_id:
                            mismatch_per_seq_id[seq_id] = full_len_mismatch_sum

                        # mismatch_threshold = (self.number_of_mismatches * 2)
                        mismatch_threshold = 6

                        if full_len_mismatch_sum > mismatch_threshold:
                            #print(
                            # f"mozaiko INFO: For {pbs_filename}, found {full_len_mismatch_sum} /
                            # mismatches for sequence {seq_id}.")
                            seq_ids_to_remove.setdefault(
                                pbs_filename, set()
                            ).add(seq_id)


        # print(seq_ids_to_remove)
        if not seq_ids_to_remove:
            print("mozaiko INFO: No sequences with mismatches found.")
            return None, None

        for pbs_filename, seq_ids in seq_ids_to_remove.items():
            amplicon_file_path = amplicon_folder / f"{pbs_filename}.fasta"
            insert_file_path = insert_folder / "filtered" / f"{pbs_filename}.fasta"

            if not amplicon_file_path.exists():
                print(f"mozaiko ERROR: Amplicon file not found: {amplicon_file_path}")
                continue

            if not insert_file_path.exists():
                print(f"mozaiko ERROR: Insert file not found: {insert_file_path}")
                continue

            self._remove_sequences_with_mismatches(amplicon_file_path, seq_ids)

            self._remove_sequences_with_mismatches(insert_file_path, seq_ids)

        if testing:
            mismatch_per_seq_id = pd.DataFrame.from_dict(
                mismatch_per_seq_id, orient="index", columns=["mismatches"]
            ).reset_index()
            mismatch_per_seq_id.rename(columns={"index": "seq_id"}, inplace=True)
            return mismatch_per_seq_id

    def _remove_sequences_with_mismatches(self, fasta_file: Path, seq_ids: set):
        """
        Remove sequences from a FASTA file based on a set of sequence IDs.

        Parameters:
        - fasta_file: Path to the FASTA file
        - seq_ids: Set of sequence IDs to remove

        Returns:
        None
        """
        # print(f"mozaiko INFO: Removing sequences with more than {self.number_of_mismatches} mismatches from {fasta_file}...")
        try:
            records_to_keep = []
            number_sequences_with_mismatches = 0

            for record in SeqIO.parse(fasta_file, "fasta"):
                fasta_seq_id = record.id.split("|")[0]

                if fasta_seq_id not in seq_ids:
                    records_to_keep.append(record)
                else:
                    number_sequences_with_mismatches += 1
                    # print(f"mozaiko INFO: Found more than allowed mismatches "
                    #   f"for sequence {fasta_seq_id} in {fasta_file}. Removing it...")

            with open(fasta_file, "w") as output_handle:
                for record in records_to_keep:
                    output_handle.write(f">{record.description}\n{record.seq}\n")

            # print(f"mozaiko INFO: Removed {number_sequences_with_mismatches} sequences from {fasta_file}")

        except Exception as e:
            print(f"mozaiko ERROR: Sanity Check Error. Could not remove sequences with more than"
                  f" {(self.number_of_mismatches * 2)} combined mismatches from {fasta_file}: {e}")

    def run_in_silico_analysis(
        self, primer_table=None, max_len_according_to_ilumina: bool = True, minimum_percentage_identity=0.75, minimum_alignment_coverage=99
    ):
        """
        This methods initiates the in-silico analysis. It does so by first veryfing if all required
        tools are installed in the machine. If installed, it requests the user to upload a table
        containing a list of primers to be evaluated and to provide a name for the folder where
        the results of the run will be stored. At last it processes the inputed primer table to
        iterate over each provided primer (row) and process the needed commands.

        Parameters:
        - primer_table: Path to the primer table
        - max_len_according_to_ilumina: Boolean to determine if the max read length should be
        calculated according to Illumina's formula
        - taxa_column_number: Number of the column in the fasta file that contains the information
        for the taxa-level we want the analysis to be performed. The count starts from 0.
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
            if self.database_fasta_file:
                self.process_commands(row, self.database_fasta_file, minimum_percentage_identity, minimum_alignment_coverage)
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

        directories_to_filter = ["all_complete_pbs", "incomplete_pbs"]
        for dir_name in directories_to_filter:
            try:
                input_path = self.output_dirs[dir_name]
                filter_sequences_by_ambiguity(
                    input_path=input_path,
                )
            except Exception as e:
                print(f"mozaiko ERROR: Error filtering {dir_name} directory: {str(e)}")

        self.intersect_PGA_relaxed_and_strict(
            self.output_dirs["all_complete_pbs"] / "filtered",
            self.output_dirs["incomplete_pbs"] / "filtered",
            output_dir_name="input_B",
        )

        print("mozaiko INFO: Number of inserts that were amplified successfully...")
        filter_inserts_path = self.output_dirs["insert"] / "filtered"
        insert_files = list(filter_inserts_path.glob("*.fasta"))
        for file in insert_files:
            number_of_sequences = self.count_sequences(file)
            print(f"    For {file.stem}, {number_of_sequences} were retained.")

        print(
            "mozaiko INFO: Number of inserts with complete PBS that were not amplified..."
        )

        print("mozaiko INFO: In-silico amplification analysis completed.")

    def process_commands(self, row: dict, input_fasta: Path, minimum_percentage_identity, minimum_alignment_coverage):
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

        self.sanity_check_on_mismatches(output_dir=self.run_dir)

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
            minimum_percentage_identity,
            minimum_alignment_coverage,
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
            minimum_percentage_identity,
            minimum_alignment_coverage,
        )

        remove_rc_suffix_from_fasta_files(results_directory=self.run_dir)

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

        print(f"mozaiko INFO: Running cutadapt command as: {' '.join(full_command)}")
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
            print(f"mozaiko ERROR: Cutadapt stdout:\n{e.stdout}")
            print(f"mozaiko ERROR: Cutadapt stderr:\n{e.stderr}")
            # print(f"mozaiko ERROR: Command: {' '.join(e.cmd)}")
            # print(f"mozaiko ERROR: Return code: {e.returncode}")
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
        minimum_percentage_identity,
        minimum_alignment_coverage,
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

        if minimum_percentage_identity is None:
            minimum_percentage_identity = 0.75  # Decimal Percentage [0 - 1.0]
        if minimum_alignment_coverage is None:
            minimum_alignment_coverage = 99  # Whole Percentage [0 - 100]

        print(
            f"mozaiko INFO: Running PGA with {minimum_alignment_coverage} coverage and {minimum_percentage_identity} identity."
        )

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
#     data = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/mozaico/data/input_data/diat-barcode-taxa_harmonized.fasta"
#     primer_table = "/home/camilababo/Documents/coding-projects/DNAquaIMG-tool/mozaico/data/input_data/diat-barcode-primers.tsv"
#     run_name = "diat-barcode-test-check-headers"
#     cutadapt = InSilicoAmplification(data, run_name=run_name)
#     cutadapt.run_in_silico_analysis(primer_table=primer_table)
