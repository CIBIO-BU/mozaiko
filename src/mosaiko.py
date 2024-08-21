#!/usr/bin/env python3

"""
This module contains the command line interface for mosaiko.
"""

import argparse
import logging

from src.reference_database.db_curation import CrabsScriptGenerator
from src.reference_database.sequence_import import CustomFastaImport


def create_parser():
    """
    Create an argument parser for the mosaiko's CLI.
    """
    parser = argparse.ArgumentParser(
        description="mosaiko: Piecing Together Complete Genetic Coverage for Biomonitoring"
    )
    parser.add_argument("-i", "--input", help="Path to the input FASTA file")
    parser.add_argument(
        "-o", "--output", help="Path to processed FASTA file (optional)"
    )
    parser.add_argument(
        "--load_custom_fasta", action="store_true", help="Load custom FASTA file."
    )
    parser.add_argument("--json_file", help="Path to the JSON file with parameters.")
    parser.add_argument(
        "--assign_tax", action="store_true", help="Assign taxonomic information."
    )
    parser.add_argument(
        "--dereplicate", action="store_true", help="Dereplicate sequences."
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose logging."
    )
    return parser


def handle_custom_fasta_import(args):
    """
    Handle the import of a custom FASTA file.
    """
    logging.info("mosaiko INFO: Initiating custom sequence import...")
    fasta_import = CustomFastaImport()
    fasta_import.read_fasta(args.input)

    if args.output:
        fasta_import.df2fasta(args.output)

    logging.info(
        f"mosaiko INFO: Processed {fasta_import.get_number_of_sequences()} sequences."
    )


def handle_taxonomic_assignment(args):
    """
    Handle the taxonomic assignment of sequences.
    """
    logging.info("mosaiko INFO: Initiating taxonomic assignment...")
    crabs_generator = CrabsScriptGenerator()
    crabs_generator.run_assign_tax_command(args.json_file)


def main():
    """
    Main function for mosaiko's CLI.
    """
    parser = create_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.load_custom_fasta and not args.input:
        parser.error(
            "mosaiko INFO: --input is required when --load_custom_fasta is used."
        )

    if args.load_custom_fasta:
        handle_custom_fasta_import(args)

    if args.assign_tax and args.json_file:
        handle_taxonomic_assignment(args)

    else:
        logging.error(
            "No JSON file specified. Please specify a JSON file with parameters."
        )
        logging.error("Exiting...")
        return


if __name__ == "__main__":
    main()
