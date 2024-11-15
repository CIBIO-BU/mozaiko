"""
This modules generates scripts to borrow methods from CRBAS. Its use is intended to assing Tax IDs,
generate lineage files and dereplicate sequences.

The CrabsScriptGenerator class contains the following methods:
    - _check_if_crabs_installed: Function to check if CRBAS is installed.
    - _load_parameters: Function to load the parameters from the JSON file.
    - _download_taxonomy_files: Functions to download the taxonomy files from the NCBI database.
    - run_assign_tax_command: Function to run the assign_tax command from CRBAS.
    - run_dereplicate_command: Function to write the script to dereplicate sequences.
"""

import json
import os
import subprocess
import sys

from src.reference_database.sequence_import import CustomFastaImport


class CrabsScriptGenerator:
    """
    This class generates commands to be executed in the shell to borrow methods from CRBAS.
    """

    def __init__(self):
        self.params = {}
        self.fasta_import = CustomFastaImport()

    def check_if_crabs_installed(self):
        """
        Function to check if CRBAS is installed.
        """
        print("mozaiko INFO: Checking if CRBAS is installed...")
        try:
            subprocess.run(
                ["crabs", "--version"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

            print("mozaiko INFO: CRBAS is installed.")

        except FileNotFoundError:
            print(
                "mozaiko INFO: CRBAS is not installed. Please install CRBAS before running this \
                    script."
            )
            print(
                "CRABS can be found at "
                + "https://github.com/gjeunen/reference_database_creator/tree/main"
            )
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

    # def _update_assign_tax_parameters(self):
    #     """
    #     Function to update needed parameters for the assign_tax by user request."
    #     """
    #     # Required Parameters
    #     print(
    #         "To assign the taxonomic IDs and generate the lineage files, "
    #         + "the following parameters are required:"
    #     )

    #     # Retrieve processed fasta file as input
    #     if self.fasta_import.fasta_file is None:
    #         print("1. Path to the input FASTA file:")
    #         fasta_file = input("Enter the path to the input FASTA file: ")
    #         self.fasta_import.read_fasta(fasta_file)

    #     self.params["input"] = self.fasta_import.fasta_file

    #     print("2. Path to the output file")
    #     self.params["output"] = input("Enter the path to the output file: ")

    #     if not self.params["output"]:
    #         print("No output file specified. Exiting...")
    #         sys.exit(1)

    #     extension = self.params["output"].split(".")[1]

    #     if extension != "tsv":
    #         print("Output file must be a TSV file. Exiting...")
    #         sys.exit(1)

    #     # Optional Parameters
    #     print("The following requirements are optional:")

    #     # print("Proceed with default taxonomic ranks? " +
    #     #       "(superkingdom+phylum+class+order+family+genus+species)" +
    #     #       "Ente yes to proceed with default ranks, or no to specify custom ranks.")

    #     # self.params["ranks"] = input("Enter yes or no: ")

    #     # if self.params["ranks"] != "yes":
    #     #     print("Disabling default taxonomic ranks")
    #     #     self.params["ranks"] = input("Enter the taxonomic ranks to use: (separated by '+')")

    #     print(
    #         "Retrieve missing taxonomic information through a web search? (yes/no, default: no)"
    #     )
    #     self.params["web"] = input("Enter yes or no: ")

    #     if self.params["web"] != "yes":
    #         print(" Web-based retrival of missing taxonomic information: Disabled")

    #     if self.params["web"] == "yes":
    #         missing_file = input(
    #             "Enter the path to the file with missing taxonomic information: "
    #         )
    #         self.params["web"] = missing_file

    #     print(
    #         "Write sequences for which no taxonomic lineage was found to a file?"
    #         + "(yes/no, default: no)"
    #     )
    #     self.params["missing"] = input("Enter yes or no: ")

    #     if self.params["missing"] != "yes":
    #         print("Writing sequences with missing taxonomic lineage: Disabled")
    #         self.params["missing"] = "no"

    #     if self.params["missing"] == "yes":
    #         self.params["missing"] = "yes"

    #     print("mozaiko INFO: Parameters updated.")

    def _download_taxonomy_files(self):
        """
        Functions to download the taxonomy files from the NCBI database.
        """
        print("mozaiko INFO: Checking if taxonomy files already exist...")

        if os.path.exists("taxonomy_files"):
            print("mozaiko INFO: Taxonomy files found. Proceeding with the analysis.")
            return

        else:
            print(
                "mozaiko INFO: Taxonomy files not found. Creating taxonomy_files folder..."
            )
            os.makedirs("taxonomy_files")

        os.chdir("taxonomy_files")

        print("mozaiko INFO: Downloading files...")
        command = "crabs db_download --source taxonomy"

        subprocess.run(command, shell=True, check=True)

        print("mozaiko INFO: Taxonomy files downloaded.")

    def _validate_json_file(self, json_file):
        """
        Function to validate the JSON file before loading.
        """
        print("mozaiko INFO: Validating JSON file...")

        if not os.path.exists(json_file):
            print("mozaiko INFO: The JSON file does not exist. Exiting...")
            raise FileNotFoundError(f"mozaiko INFO: {json_file} not found.")

        _, file_extension = os.path.splitext(json_file)

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "json":
            raise ValueError(f"mozaiko INFO: Invalid file type: {file_extension}. Expected a .json file.")

    def _load_parameters(self, json_file):
        """
        Function to load the parameters from JSON files.
        """
        self._validate_json_file(json_file)

        print(
            "mosaiko INFO: Validation Successful. Loading parameters from JSON file..."
        )
        with open(json_file, encoding="UTF-8") as file:
            self.params = json.load(file)

        print("mosaiko INFO: Parameters loaded.")

    def run_assign_tax_command(self, json_file):
        """
        Function to run the assign_tax command from CRBAS.
        """
        self.check_if_crabs_installed()

        self._download_taxonomy_files()

        self._load_parameters(json_file)

        print(" mozaiko INFO: All set. Running taxonomy assignment task..")

        # print("Generating script to assign taxonomic IDs...")

        command = (
            f"crabs assign_tax --input {self.params['input']} --output {self.params['output']}"
            f" --acc2tax {self.params['acc2tax']} --taxid {self.params['taxid']}"
            f" --name {self.params['name']} --web {self.params['web']}"
            f" --rank {self.params['ranks']} --missing {self.params['missing']}"
        )

        # print(f"Script generated: {command}")

        # print("Running script...")
        subprocess.run(command, shell=True, check=True)

    def run_dereplicate_command(self, json_file):
        """
        Function to write the script to dereplicate sequences.
        """
        self.check_if_crabs_installed()

        self._load_parameters(json_file)

        _, file_extension = os.path.splitext(self.params["input"])

        file_extension = file_extension.lstrip(".")

        if file_extension.lower() != "tsv":
            print("mozaiko INFO: Input file must be a TSV file. Exiting...")
            sys.exit(1)

        print("mozaiko INFO: All set. Running dereplication...")

        # print("Generating script to dereplicate sequences...")

        command = (
            f"crabs dereplicate --input {self.params['input']}"
            f" --output {self.params['output']}"
            f" --method {self.params['method']}"
            f" --ranks {self.params['ranks']}"
        )

        # print(f"Script generated: {command}")

        # print("Running script...")
        subprocess.run(command, shell=True, check=True)
