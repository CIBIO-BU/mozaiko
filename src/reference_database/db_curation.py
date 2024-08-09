"""
This modules generates scripts to borrow methods from CRBAS. Its use is intended to assing Tax IDs,
generate lineage files and dereplicate sequences.
"""

import json
import subprocess
import sys
import os

from reference_database.sequence_import import CustomFastaImport


class CrabsScriptGenerator:
    """
    This class generates commands to be executed in the shell to borrow methods from CRBAS.
    """

    def __init__(self):
        _base_path = os.path.dirname(os.path.abspath(__file__))
        self.json_file = os.path.join(_base_path, "crabs_parameters.json")
        self.params = {}
        self.fasta_import = CustomFastaImport()

    def _check_if_crabs_installed(self):
        """
        Function to check if CRBAS is installed.
        """
        print("Checking if CRBAS is installed...")
        try:
            subprocess.run(["crabs", "-h"],
                           check=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except  FileNotFoundError:
            print(
                "CRBAS is not installed. Please install CRBAS before running this script."
            )
            print(
                "CRABS can be found at "
                + "https://github.com/gjeunen/reference_database_creator/tree/main"
            )
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

    def _load_parameters(self):
        """
        Function to load the parameters from the JSON file.
        """
        print("Loading parameters from JSON file...")
        with open(self.json_file, encoding="UTF-8") as file:
            self.params = json.load(file)

    def _update_parameters(self):
        """
        Function to request user to input the needed parameters in the JSON file.
        """
        # Required Parameters
        print(
            "To assign the taxonomic IDs and generate the lineage files, "
            + "the following parameters are required:"
        )

        # Retrieve processed fasta file as input
        self.params["input"] = self.fasta_import.fasta_file

        print("2. Path to the output file")
        self.params["output"] = input("Enter the path to the output file: ")

        # Optional Parameters
        print("The following requirements are optional:")

        # print("Proceed with default taxonomic ranks? " +
        #       "(superkingdom+phylum+class+order+family+genus+species)" +
        #       "Ente yes to proceed with default ranks, or no to specify custom ranks.")

        # self.params["ranks"] = input("Enter yes or no: ")

        # if self.params["ranks"] != "yes":
        #     print("Disabling default taxonomic ranks")
        #     self.params["ranks"] = input("Enter the taxonomic ranks to use: (separated by '+')")

        print(
            "Retrieve missing taxonomic information through a web search? (yes/no, default: no)"
        )
        self.params["missing"] = input("Enter yes or no: ")

        if self.params["missing"] != "yes":
            print("Disabling missing taxonomic information retrieval")
            self.params["missing"] = "no"

        if self.params["missing"] == "yes":
            self.params["missing"] = "yes"

        print(
            "Write sequences for which no taxonomic lineage was found to a file?"
            + "(yes/no, default: no)"
        )
        self.params["missing"] = input("Enter yes or no: ")

        if self.params["missing"] != "yes":
            print("Disabling writing of sequences with missing taxonomic lineage")
            self.params["missing"] = "no"

        if self.params["missing"] == "yes":
            self.params["missing"] = "yes"

        print("Parameters updated.")

    def _download_taxonomy_files(self):
        """
        Functions to download the taxonomy files from the NCBI database.
        """
        print("Downloading taxonomy files...")

        print("Checking if taxonomy files already exist...")

        if os.dir.exists("taxonomy_file"):
            print("Taxonomy files found. Proceeding with the analysis.")
            return

        print("Taxonomy files not found. Downloading files...")
        command = "crabs db_download --source taxonomy"

        subprocess.run(command, shell=True, check=True)

        print("Taxonomy files downloaded.")

    def run_assign_tax_command(self):
        """
        Function to run the assign_tax command from CRBAS.
        """
        self._check_if_crabs_installed()

        self._load_parameters()

        self._update_parameters()

        print("Generating script to assign taxonomic IDs...")
        command = (
            f"crabs assign_tax --input {self.params['input']} --output {self.params['output']}"
            f" --acc2tax {self.params['acc2tax']} --taxid {self.params['taxid']}"
            f" --name {self.params['name']} --missing {self.params['missing']}"
        )

        print(f"Script generated: {command}")

        print("Running script...")
        subprocess.run(command, shell=True, check=True)

    def generate_dereplicate_script(self):
        """
        Function to write the script to dereplicate sequences.
        """
