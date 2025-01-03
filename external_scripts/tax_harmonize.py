import sys
import pandas as pd
from pygbif import species
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm
import argparse


def parse_database(database: object) -> list:
    """Function to read the database in form of a .fasta file.

    Args:
        database (object): Path to the fasta file to be read.

    Returns:
        list: Generator of dereplicated species names found in the input file.
    """
    # extract taxa names from the headers of fasta file
    taxa_names = [
        header.split(" | ")[-1]
        for (header, seq) in SimpleFastaParser(open(Path(database)))
    ]

    # remove duplicate entries
    taxa_names = [name for name in set(taxa_names)]

    # return the taxa names
    return taxa_names


def query_api(taxa_names: list) -> dict:
    """Function to query the GBIF API.

    Args:
        taxa_names (list): List of taxa names to query.

    Returns:
        dict: Dict in the form of {taxa_name: (scientificName, rank)}
    """

    # gather the results here
    scientific_names_dict = {}

    for name in tqdm(taxa_names):
        r = species.name_backbone(name)

        if r["matchType"] == "NONE":
            scientific_names_dict[name] = ("-", "-")
        else:
            scientific_names_dict[name] = (r["scientificName"], r["rank"])

    return scientific_names_dict


def update_fasta(
    scientific_names: dict, database: str, database_name: str, database_dir: str
):
    """Function to save results to an updated .fasta

    Args:
        scientific_names (dict): Dict with scientific names and taxon names.
        database (str): Path to the database
        database_name (str): Name ofthe database
        database_dir (str): Directory of the database
    """
    # generate the output name
    output_name = database_dir.joinpath("{}_harmonized.fasta".format(database_name))

    with open(Path(database), "r", encoding="utf-8") as in_handle:
        with open(output_name, "a", encoding="utf-8") as out_handle:
            for header, seq in SimpleFastaParser(in_handle):
                header = header.split(" | ")
                header_start = header[0]
                header_species = header[1]
                scientificName = scientific_names[header_species][0]
                rank = scientific_names[header_species][1]
                out_handle.write(
                    ">{}|{}|{}|{}\n{}\n".format(
                        header_start, header_species, scientificName, rank, seq
                    )
                )


def parse_otl(otl: object) -> list:
    """Function to parse the otl file.

    Args:
        otl (object): Path to the otl file.

    Returns:
        list: List of taxa names
    """
    otl = pd.read_csv(Path(otl), sep="\t")
    otl = list(set([taxa_name for taxa_name in otl["taxa"]]))

    return otl


def update_otl(
    otl_scientific_names: dict, otl: object, otl_name: object, otl_dir: object
):
    """Function to save a new otl with scientific names and ranks

    Args:
        otl_scientific_names (dict): Dict that holds the original names and theAPI responses.
        otl (object): Path to the otl
        otl_name (object): Name of the otl
        otl_dir (object): Directory of the otl
    """
    # generate the output name
    output_name = otl_dir.joinpath("{}_harmonized.tsv".format(otl_name))

    # read the input
    otl = pd.read_csv(Path(otl), sep="\t")

    # update the otl
    scientific_names = {key: value[0] for key, value in otl_scientific_names.items()}
    ranks = {key: value[1] for key, value in otl_scientific_names.items()}

    # update the otl
    otl["scientific_name"] = otl["taxa"].map(scientific_names)
    otl["rank"] = otl["taxa"].map(ranks)

    # save the updated otl
    otl.to_csv(output_name, sep="\t", index=False)

def harmonize_database(database):
    fasta_path = Path(database)
    taxa_names = parse_database(database)
    scientific_names = query_api(taxa_names)
    update_fasta(scientific_names, database, fasta_path.stem, fasta_path.parent)

def harmonize_otl(otl):
    otl_path = Path(otl)
    otl_taxa_names = parse_otl(otl)
    otl_scientific_names = query_api(otl_taxa_names)
    update_otl(otl_scientific_names, otl, otl_path.stem, otl_path.parent)

def main():
    parser = argparse.ArgumentParser(description="Harmonize taxonomy of database and/or OTL")
    parser.add_argument("--database", "-d", help="Path to database file (.fasta)")
    parser.add_argument("--otl", "-o", help="Path to OTL file (.tsv)")

    args = parser.parse_args()

    if not args.database and not args.otl:
        parser.error("At least one of --database or --otl must be provided")

    if args.database:
        harmonize_database(args.database)

    if args.otl:
        harmonize_otl(args.otl)

if __name__ == "__main__":
    main()