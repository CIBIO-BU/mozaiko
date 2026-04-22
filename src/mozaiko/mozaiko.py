#!/usr/bin/env python3

"""
This module contains the command line interface for mozaiko.
"""

import argparse
import logging

from mozaiko.in_silico_analysis.amplification import InSilicoAmplification
from mozaiko.reference_database.db_curation import CrabsScriptGenerator
from mozaiko.reference_database.sequence_import import CustomFastaImport

__version__ = "0.1.0"


def create_parser():
    """
    Create an argument parser for mozaiko's CLI.
    """
    parser = argparse.ArgumentParser(
        description="mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring",
        epilog="""
Examples:
  # Pre-process a FASTA database
  %(prog)s --pre_process_db -i input.fasta -o output.tsv

  # Pre-process with taxonomic harmonization
  %(prog)s --pre_process_db --harmonized -i input.fasta -o output.tsv

  # Assign taxonomy to sequences
  %(prog)s --assign_tax --json_file config.json

  # Dereplicate sequences
  %(prog)s --dereplicate --json_file derep_config.json

  # Run in-silico PCR analysis
  %(prog)s --in_silico_analysis -i database.fasta --run_name my_analysis --primer_table primers.tsv

For more information, visit: https://github.com/CIBIO-BU/mozaiko
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False  # We'll add a custom help group below
    )

    # Help & Version options (shown first in --help output)
    info_group = parser.add_argument_group('Information Options')
    info_group.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    info_group.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program version number and exit"
    )

    # Input/Output options
    io_group = parser.add_argument_group('Input/Output Options')
    io_group.add_argument(
        "-i", "--input",
        help="Path to the input FASTA file (required for --pre_process_db and --in_silico_analysis)",
        metavar="FILE"
    )
    io_group.add_argument(
        "-o", "--output",
        help="Path where the processed FASTA file will be saved (required for --pre_process_db)",
        metavar="FILE"
    )

    # Workflow Processing options
    process_group = parser.add_argument_group('Workflow Processing Options')
    process_group.add_argument(
        "--pre_process_db",
        action="store_true",
        help="Pre-process input database (cleans and standardizes FASTA file)"
    )
    process_group.add_argument(
        "--harmonized",
        action="store_true",
        help="Pre-processed databases as if taxonomic harmonization was previously applied (use with --pre_process_db)"
    )
    process_group.add_argument(
        "--assign_tax",
        action="store_true",
        help="Assign taxonomic information using CRABS (requires --json_file)"
    )
    process_group.add_argument(
        "--dereplicate",
        action="store_true",
        help="Remove duplicate sequences using CRABS method (requires --json_file)"
    )

    # Analysis options
    analysis_group = parser.add_argument_group('Analysis Options')
    analysis_group.add_argument(
        "--in_silico_analysis",
        action="store_true",
        help="Perform in-silico PCR amplification analysis (requires -i, --run_name, and --primer_table)"
    )
    analysis_group.add_argument(
        "--run_name",
        help="Name for the analysis results folder (required for --in_silico_analysis)",
        metavar="NAME"
    )
    analysis_group.add_argument(
        "--primer_table",
        help="Path to TSV file containing primer sequences to evaluate (required for --in_silico_analysis)",
        metavar="FILE"
    )

    # Configuration options
    config_group = parser.add_argument_group('Configuration Options')
    config_group.add_argument(
        "--json_file",
        help="Path to JSON configuration file with parameters for CRABS operations (--assign_tax, --dereplicate)",
        metavar="FILE"
    )
    config_group.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed logging output for debugging"
    )

    return parser


def database_pre_process(args):
    """
    Handle the import of a custom FASTA file.
    """
    print("mozaiko INFO: Initiating database pre-process...")

    try:
        fasta_import = CustomFastaImport()
        fasta_import.read_fasta(args.input)

        print(
            f"mozaiko INFO: Processed {fasta_import.get_number_of_sequences()} sequences."
        )

        if args.harmonized:
            print(f"mozaiko INFO: Processing taxonomic harmonization.")
            fasta_import.pre_process_harmonized_fasta_database()

        print(f"mozaiko INFO: Saving output to {args.output}...")
        fasta_import.df2csv(args.output)

    except Exception as e:
        logging.error(f"mozaiko ERROR: Failed to process the FASTA file: {e}")
        raise


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
    in_silico_generator = InSilicoAmplification(database_fasta_file=args.input, run_name=args.run_name)
    in_silico_generator.run_in_silico_analysis(primer_table=args.primer_table)


def main():
    """
    Main function for mozaiko's CLI.
    """
    parser = create_parser()
    args = parser.parse_args()

    # Print help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        return

    print(f"Parsed args: {args}")

    if args.pre_process_db and args.dereplicate:
        parser.error(
            "mozaiko INFO: The options --pre_process_db and --dereplicate cannot be used together."
        )

    # Verbose logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load custom FASTA file
    if args.pre_process_db and not args.input and not args.output:
        parser.error(
            "mozaiko INFO: Please specify a FASTA file with parameter --input and an output file name with --output."
        )

    if args.pre_process_db:
        print("Pre-processing database...")
        database_pre_process(args)

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
            "mozaiko INFO: No JSON file specified. Please specify a JSON file with parameters --json_file for dereplication. Refer to https://github.com/CIBIO-BU/mozaiko/blob/main/data/test_data/test_dereplication.json for an example."
        )
        logging.error("Exiting...")
        return

    # In-silico Analysis
    if args.in_silico_analysis and args.input and args.run_name and args.primer_table:
        handle_in_silico_analysis(args)

    elif args.in_silico_analysis:
        logging.error(
            "mozaiko INFO: Please specify a FASTA file with parameter --input, a name for the results folder with --run_name, and a table containing primer information with --primer_table."
        )
        return


if __name__ == "__main__":
    main()