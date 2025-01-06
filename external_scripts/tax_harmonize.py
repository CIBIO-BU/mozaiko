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
    """
    Query GBIF API for unique taxa only.
    """
    scientific_names_dict = {}
    taxonomic_fields = ['scientificName', 'rank', 'kingdom', 'phylum',
                       'class', 'order', 'family', 'genus']

    for name in tqdm(taxa_names, desc="Querying GBIF API"):
        if name not in scientific_names_dict:
            r = species.name_backbone(name, strict=False)
            if r["matchType"] == "NONE":
                scientific_names_dict[name] = tuple(["-"] * len(taxonomic_fields))
            else:
                scientific_names_dict[name] = tuple(r.get(field, "-")
                                                  for field in taxonomic_fields)

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
    output_name = database_dir.joinpath(f"{database_name}_harmonized.fasta")
    with open(Path(database), "r", encoding="utf-8") as in_handle:
        with open(output_name, "a", encoding="utf-8") as out_handle:
            for header, seq in SimpleFastaParser(in_handle):
                header = header.split(" | ")
                header_start, header_species = header[0], header[1]
                tax_info = scientific_names[header_species]
                header_str = f">{header_start}|{header_species}|" + "|".join(tax_info)
                out_handle.write(f"{header_str}\n{seq}\n")


def parse_otl(otl: object) -> tuple:
    """
    Return unique taxa and original dataframe.
    """
    otl_df = pd.read_csv(Path(otl), sep="\t")
    unique_taxa = list(otl_df['taxa'].unique())
    return unique_taxa, otl_df

def update_otl(otl_scientific_names: dict, otl_df: pd.DataFrame, otl_name: str, otl_dir: object):
    """
    Update OTL using vectorized operations.
    """
    output_name = otl_dir.joinpath(f"{otl_name}_harmonized.tsv")
    fields = ['scientificName', 'rank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']

    # Create a mapping DataFrame
    mapping_data = pd.DataFrame.from_dict(otl_scientific_names, orient='index', columns=fields)

    # Merge with original dataframe
    result = otl_df.merge(mapping_data, left_on='taxa', right_index=True, how='left')
    result.to_csv(output_name, sep="\t", index=False)

def harmonize_database(database):
    fasta_path = Path(database)
    taxa_names = parse_database(database)
    scientific_names = query_api(taxa_names)
    update_fasta(scientific_names, database, fasta_path.stem, fasta_path.parent)

def harmonize_otl(otl):
    otl_path = Path(otl)
    unique_taxa, otl_df = parse_otl(otl)
    print(f"Found {len(unique_taxa)} unique taxa out of {len(otl_df)} total entries")

    otl_scientific_names = query_api(unique_taxa)
    update_otl(otl_scientific_names, otl_df, otl_path.stem, otl_path.parent)

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