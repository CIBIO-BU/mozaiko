#!/usr/bin/env python3

"""
This module contains the command line interface for mozaiko.
"""

import argparse
import logging

from src.in_silico_analysis.amplification import InSilicoAmplification
from src.reference_database.db_curation import CrabsScriptGenerator
from src.reference_database.sequence_import import CustomFastaImport


def create_parser():
    """
    Create an argument parser for the mozaiko's CLI.
    """
    parser = argparse.ArgumentParser(
        description="mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring"
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
    parser.add_argument(
        "--in_silico_analysis", action="store_true", help="Run in-silico analysis."
    )
    return parser


def handle_custom_fasta_import(args):
    """
    Handle the import of a custom FASTA file.
    """
    print("mozaiko INFO: Initiating custom sequence import...")
    fasta_import = CustomFastaImport()
    fasta_import.read_fasta(args.input)

    print(
        f"mozaiko INFO: Processed {fasta_import.get_number_of_sequences()} sequences."
    )

    if args.output:
        fasta_import.df2fasta(args.output)


def handle_taxonomic_assignment(args):
    """
    Handle the taxonomic assignment of sequences.
    """
    print("mozaiko INFO: Initiating taxonomic assignment...")
    crabs_generator = CrabsScriptGenerator()
    crabs_generator.run_assign_tax_command(args.json_file)


def handle_dereplication(args):
    """
    Handle the dereplication of sequences.
    """
    print("mozaiko INFO: Initiating sequence dereplication...")
    crabs_generator = CrabsScriptGenerator()
    crabs_generator.run_dereplicate_command(args.json_file)


def handle_in_silico_analysis(args):
    """
    Handle the in-silico analysis process.
    """
    print("mozaiko INFO: Initiating in-silico analysis...")
    in_silico_generator = InSilicoAmplification(args.input)
    in_silico_generator.run_in_silico_analysis()


def main():
    """
    Main function for mozaiko's CLI.
    """
    parser = create_parser()
    args = parser.parse_args()
    print(f"Parsed args: {args}")

    # Verbose logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load custom FASTA file
    if args.load_custom_fasta and not args.input:
        parser.error(
            "mozaiko INFO: No FASTA file specified. Please specify a FASTA file with parameter --input."
        )

    if args.load_custom_fasta and args.input:
        print("Loading custom FASTA...")
        handle_custom_fasta_import(args)

    else:
        print("Skipping custom FASTA load...")

    # Assign taxonomic information
    if args.assign_tax and args.json_file:
        handle_taxonomic_assignment(args)

    elif args.assign_tax:
        logging.error(
            "mozaiko INFO: No JSON file specified. Please specify a JSON file with parameter --json_file."
        )
        logging.error("Exiting...")
        return

    # Dereplicate sequences
    if args.dereplicate and args.json_file:
        handle_dereplication(args)

    elif args.dereplicate:
        logging.error(
            "mozaiko INFO: No JSON file specified. Please specify a JSON file with parameter --json_file."
        )
        logging.error("Exiting...")
        return

    # In-silico Analysis
    if args.in_silico_analysis and args.input:
        handle_in_silico_analysis(args)

    elif args.in_silico_analysis:
        logging.error(
            "mozaiko INFO: No FASTA file specified. Please specify a FASTA file with parameter --input."
        )
        return


if __name__ == "__main__":
    main()
