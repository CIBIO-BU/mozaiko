"""
This is a stand-alone script to convert a TSV file to a FASTA file; and remove alignment dashes from FASTA files.
"""

import pandas as pd


def tsv2fasta(tsv_file, output_fasta_file):
    """
    This function converts a TSV file to a FASTA file.

    Parameters
    tsv_file (str): The input TSV file.
    output_fasta_file (str): The output FASTA file.
    """
    df = pd.read_csv(tsv_file, sep="\t")
    sequences = df["sequence"]
    seq_ids = df["seq_ID"]
    species = df["species"]
    with open(output_fasta_file, "w", encoding="UTF-8") as f:
        for i, sequence in enumerate(sequences):
            f.write(f">{seq_ids[i]} {species[i]}\n")
            f.write(f"{sequence}\n")

    print("TSV file converted to FASTA file successfully.")


def remove_alignment_dashes(fasta_file, output_fasta_file):
    """
    This function removes alignment dahses from FASTA files.

    These alignment dashes '-' are added to the sequences to align their lengths.

    Parameters
    fasta_file (str): The input FASTA file.
    output_fasta_file (str): The output FASTA file.
    """
    with open(fasta_file, "r", encoding="UTF-8") as f:
        lines = f.readlines()
    with open(output_fasta_file, "w", encoding="UTF-8") as f:
        for i, line in enumerate(lines):
            if i % 2 == 0:
                f.write(line)
            else:
                f.write(line.replace("-", ""))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert a TSV file into a FASTA file")
    parser.add_argument("--tsv_file", "-t", help="The input TSV file")
    parser.add_argument("--fasta_file", "-f", help="The input FASTA file.")
    parser.add_argument("--output_fasta_file", "-o", help="The output FASTA file.")
    parser.add_argument(
        "--task",
        "-a",
        choices=["tsv2fasta", "remove_alignment_dashes"],
        help="The task to perform: tsv2fasta or remove_alignment_dashes",
    )
    args = parser.parse_args()

    if args.task == "tsv2fasta":
        tsv2fasta(args.tsv_file, args.output_fasta_file)
    elif args.task == "remove_alignment_dashes":
        remove_alignment_dashes(args.fasta_file, args.output_fasta_file)
    else:
        print(
            "Invalid task."
            + "Please provide a valid task ('tsv2fasta', 'remove_alignment_dashes') "
            + "with the use of -a or --task argument."
        )
