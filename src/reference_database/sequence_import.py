"""
This module contains the CustomFastaImport class which is responsible for handling and
transforming custom data into desired inputs/outputs.

The CustomFastaImport class contains the following methods:
- _validate_input: Validates the input provided by the user.
- _add_taxids: Joins the TaxIDs of sequences in the data.
- check_for_taxids: Checks if fasta file contains TaxIDs.
- _request_lineage_file: Requests users to upload lineage file if no taxids are found in fasta file.
- read_fasta: Reads a fasta file.
- get_taxids: Returns the TaxIDs of the sequences in the data.
- get_number_of_sequences: Returns the number of sequences in the data.
- get_sequence_lengths: Returns the lengths of the sequences in the data.
- get_sequence_ids: Returns the sequence IDs in the data.
- get_sequences: Returns the sequences in the data.
- df2csv: Write the data frame to a csv file.
"""

import argparse
import os
import re

import pandas as pd
from Bio import SeqIO


class CustomFastaImport:
    """
    Handles and transforms fasta into desired outputs.
    """

    def __init__(self, data=None):
        """
        Initializes the CustomFastaImport class.
        """
        self.data = data
        self.lineage_file_loader = LineageFileLoader()
        self.lineage_file = None
        self.fasta_file = None

    def _validate_input(self, input_file):
        """
        Validades the input provided by the user.

        Parameters
        input_file (str): Path to fasta file.
        """

        if not isinstance(input_file, str):
            raise ValueError("Directory must be a string.")

        if not os.path.exists(input_file):
            raise FileNotFoundError("Input file does not exist in the directory.")

        if not input_file.endswith(".fasta"):
            raise ValueError("Input file must be a fasta file.")

        if os.path.getsize(input_file) == 0:
            raise ValueError("Input file is empty.")

    def add_taxids(self, input_file):
        """
        Joins the TaxIDs of sequences in the data.

        Returns
        pd.DataFrame
        """

        with open(input_file, "r", encoding="UTF-8") as fasta_file:
            records = SeqIO.parse(fasta_file, "fasta")
            taxids = []

            for record in records:
                match = re.search(r"(?<=taxid=)([0-9]+)", record.description)
                if match:
                    taxids.append(match.group(1))

            self.data["taxid"] = taxids

        return self.data

    def check_for_taxids(self, input_file):
        """
        Checks if fasta file contains TaxIDs.
        """
        with open(input_file, "r", encoding="UTF-8") as fasta_file:
            taxid_found = False
            for record in SeqIO.parse(fasta_file, "fasta"):

                if "taxid" in record.description.lower():
                    self.add_taxids(input_file)
                    taxid_found = True
                    break

            if not taxid_found:
                print(
                    "mozaiko INFO: No TaxIDs found in the FASTA file. Starting lineage file upload process."
                )
                self.lineage_file = self.lineage_file_loader.load_lineage_file()

    def read_fasta(self, input_file, sep="|", check_taxid=False, harmonized: bool = False):
        """
        Reads a fasta file.

        Parameters
        input_file (str): Path to the fasta file.

        Returns
        pd.DataFrame
        """

        self._validate_input(input_file)

        with open(input_file, "r", encoding="UTF-8") as fasta_file:
            records = SeqIO.parse(fasta_file, "fasta")
            data_dict = {"seq_id": [], "sequence": [], "length": []}
            if not check_taxid:
                data_dict["taxa_info"] = []

            for seq in records:
                name, sequence, description = seq.id, str(seq.seq), seq.description
                seq_len = len(sequence)

                name = name.split("|")[0]
                data_dict["seq_id"].append(name)
                data_dict["sequence"].append(sequence)
                data_dict["length"].append(seq_len)

                if not check_taxid:
                    description_parts = description.split(sep)
                    if harmonized == True:
                        taxa_info = (
                            description_parts[2].strip()
                            if len(description_parts) > 1
                            else ""
                        )
                    else:
                        taxa_info = (
                            description_parts[1].strip()
                            if len(description_parts) > 1
                            else ""
                        )
                    data_dict["taxa_info"].append(taxa_info)

            self.data = pd.DataFrame(data_dict)

        if check_taxid:
            self.check_for_taxids(input_file)

        return self.data

    def get_taxids(self):
        """
        Returns the TaxIDs of the sequences in the data.

        Returns
        pd.Series
        """

        taxids = self.data["taxid"]
        taxids = taxids.to_list()

        return taxids

    def get_number_of_sequences(self):
        """
        Returns the number of sequences in the data.

        Returns
        int
        """

        return len(self.data)

    def get_sequence_lengths(self):
        """
        Returns the lengths of the sequences in the data.

        Returns
        pd.Series
        """

        seq_lengths = self.data["length"]
        seq_lengths = seq_lengths.to_list()

        return seq_lengths

    def get_sequence_ids(self):
        """
        Returns the sequence IDs in the data.

        Returns
        pd.Series
        """

        seq_ids = self.data["seq_id"]
        seq_ids = seq_ids.to_list()

        return seq_ids

    def get_sequences(self):
        """
        Returns the sequences in the data.
        """

        seq_list = self.data["sequence"]
        seq_list = seq_list.to_list()

        return seq_list

    def df2csv(self, output_name: str = "data/output_data/processed_input_fasta.csv"):
        """
        Write the data frame to a csv file.
        """

        self.data.to_csv(output_name, index=False)

    def df2fasta(
        self, output_name: str = "data/output_data/processed_input_fasta.fasta"
    ):
        """
        Write the data frame to a fasta file.
        """

        if not output_name.lower().endswith((".fasta", ".fa")):
            raise ValueError(
                f"mozaiko ERROR: Invalid output file name. File must have a '.fasta' extension."
            )

        output_dir = os.path.dirname(output_name)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(output_name, "w") as file:
            for index, row in self.data.iterrows():
                file.write(f">{row['seq_id']}\n{row['sequence']}\n")

        self.fasta_file = output_name

        return self.fasta_file

    def get_mapping_between_seq_id_taxonomy(self):
        """
        This method creates a mapping between the seq-id and the taxonomy present in the inputed
        FASTA file.

        Returns:
        - seq_id_taxonomy_dict: Dict
            A dictionary with seq-ids as keys and the taxonomy information as values.
        """
        seq_id_taxonomy_dict = dict(zip(self.data["seq_id"], self.data["taxa_info"]))

        return seq_id_taxonomy_dict


class LineageFileLoader:
    """
    This class is responsible for loading the lineage file.
    """

    def __init__(self):
        """
        Initializes the LineageFileLoader class.
        """
        self.header_requirements = [
            "seq_id",
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "subkingdom",
            "kingdom",
            "empire",
        ]
        self.str_requirements = ", ".join(self.header_requirements)
        self.lineage_file = None
        self.help_message_template = """
        --------------------------------- help message ---------------------------------------
        Taxonomic information is needed to continue the in-silico analysis.
        The inputed FASTA file does not contain Taxonomic IDs.
        Please upload a TSV file containing the sequence IDs and taxonomic lineage.

        The TSV file must contain the columns:
        [{columns}].
        Please make sure columns follow the same order before uploading.

        For the taxonomic levels where no assignment is available, please leave the cells blank.

        If your FASTA file does contain Taxonomic IDs for all sequences, please make sure these
        are present in each sequence header.
        For correct reading, these must be identified with 'taxid=' beforehand.
        For example: 'CM074756.1|taxid=8481'.
        ---------------------------------------------------------------------------------------
        """

    def _print_help_message(self):
        """
        Prints help message to guide users on how to upload the lineage
        """
        print(self.help_message_template.format(columns=self.str_requirements))

    def _validate_file(self, input_file):
        """
        Validates the input file provided by the user.

        Parameters
        input_file (str): Path to the lineage file.

        Returns
        str: Error message if validation
        """

        if not input_file.endswith(".tsv"):
            return "File must be a TSV file. Please try again."

        if not os.path.isfile(input_file):
            return "Invalid input. Please try again."

        if os.path.getsize(input_file) == 0:
            return "File is empty. Please try again."

        if not os.path.exists(input_file):
            return "File does not exist in the directory. Please try again."

        return None

    def read_lineage_file(self, input_file):
        """
        Reads the lineage file.

        Parameters
        input_file (str): Path to the lineage file.

        Returns
        pd.DataFrame
        """
        try:
            lineage_file = pd.read_csv(input_file, header=0, sep="\t")

        except pd.errors.ParserError as e:
            raise ValueError(f"Error reading the file: {e}. Please try again.") from e

        except Exception as e:
            raise RuntimeError(
                f"An unexpected error occurred while reading the file: {e}."
                + "Please try again."
            ) from e

        lineage_header = [column.lower() for column in lineage_file.columns]

        if lineage_header != self.header_requirements:
            raise ValueError(
                f"Columns in TSV file do not match the requirements: [{self.str_requirements}]."
                " Please try again."
            )

        self.lineage_file = lineage_file
        return self.lineage_file

    def load_lineage_file(self):
        """
        Requests users to upload lineage file if no taxids are found in fasta file.
        Prints help message to guide users on how to upload the lineage file.
        Allows users to exit the operation.

        Returns
        pd.DataFrame
        """

        self._print_help_message()

        while True:
            input_file = input(
                "Please type the directory of the TSV file to upload "
                + "(or 'exit' to quit this operation): "
            )

            if input_file.strip().lower() == "exit":
                print("Operation canceled. Data currently in memory: ")
                return None

            if input_file == "":
                print("Error: No input_file provided. Please try again.")
                continue

            error_message = self._validate_file(input_file)

            if error_message:
                print(f"Error: {error_message}")
                user_input = (
                    input(
                        "Type 'help' for more information, or press Enter to try again: "
                    )
                    .strip()
                    .lower()
                )

                if user_input == "help":
                    self._print_help_message()

                continue

            self.read_lineage_file(input_file)
            print("File uploaded successfully.")

            break

        return self.lineage_file


def create_parser():
    """
    Creates an ArgumentParser object to parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Tool to import custom fasta files for downstream analysis."
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input fasta file.",
    )

    return parser
