#!/usr/bin/env python3

import argparse

from src.reference_database.db_curation import CrabsScriptGenerator
from src.reference_database.sequence_import import CustomFastaImport


def create_parser():
    parser = argparse.ArgumentParser(
        description="mosaiko: Piecing Together Complete Genetic Coverage for Biomonitoring"
    )
    parser.add_argument("-i", "--input", help="Path to the input FASTA file")
    parser.add_argument(
        "-o", "--output", help="Path to output files (for FASTA and CSV)"
    )
    parser.add_argument("--csv", action="store_true", help="Generate CSV output")
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
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    # Ensure --input is required only when --load_custom_fasta is specified
    if args.load_custom_fasta and not args.input:
        parser.error("--input is required when --load_custom_fasta is used.")

    if args.load_custom_fasta:
        print("Initiating custom sequence import...")
        fasta_import = CustomFastaImport()
        fasta_import.read_fasta(args.input)

        if args.output:
            fasta_import.df2fasta(args.output)

        if args.csv:
            csv_output = (
                args.output.rsplit(".", 1)[0] + ".csv"
                if args.output
                else "processed_input_fasta.csv"
            )
            fasta_import.df2csv(csv_output)

        print(f"Processed {fasta_import.get_number_of_sequences()} sequences.")

    crabs_generator = CrabsScriptGenerator()

    if args.assign_tax and not args.json_file:
        print(
            "mosaiko INFO: No JSON file specified. Please specify a JSON file with parameters."
        )
        print("Exiting...")
        return

    if args.assign_tax and args.json_file:
        print("mosaiko INFO: Initiating taxonomic assignment ...")
        crabs_generator.run_assign_tax_command(args.json_file)


if __name__ == "__main__":
    main()
