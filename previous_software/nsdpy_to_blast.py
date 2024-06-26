
import argparse
import os
import pandas as pd
import sys


#".blast_files/NSDPY_results/2023-03-13_11-19-35/tsv/sequences.tsv"


def nspdy_to_blast(sequence_table):
    """
    Convert NSDPY results to BLAST format.

    Parameters:
    sequence_table (str): Path to the sequence table file (TSV format)

    Returns:
    str: Path to the output fasta file

    """
    print("Converting NSDPY results to BLAST format")
    
    # read the sequence table
    try:
        read_file = pd.read_csv(sequence_table, sep="\t")
        print("Read sequence table successfully.")
    except Exception as e:
        print(f"Failed to read sequence table: {e}")
        return

    # write the sequence table to a fasta file
    directory = os.path.dirname(sequence_table)
    base_name = os.path.basename(directory)

    output_fasta_path = os.path.join(directory, f"reference_{base_name}.fasta")
    
    with open(output_fasta_path, "w") as output_fasta:
        for i in range(len(read_file.index)):
            output_fasta.write(">{};taxid={};\n{}\n".format(read_file.loc[i, "SeqID"],read_file.loc[i, "TaxID"], read_file.loc[i, "sequence"]))

    print(f"NSDPY results have been converted to BLAST format. Output file: {output_fasta_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert NSDPY results to BLAST format")
    parser.add_argument("sequence_table", help="Path to the sequence table file (TSV format)")

    args = parser.parse_args()
    nspdy_to_blast(args.sequence_table)