#!/usr/bin/env python3

"""
This module contains the command line interface for mozaiko.
"""

import argparse
import logging

from mozaiko.in_silico_analysis.amplification import InSilicoAmplification
from mozaiko.reference_database.db_curation import CrabsScriptGenerator
from mozaiko.reference_database.sequence_import import CustomFastaImport
from mozaiko.marker_scoring.metrics_system import MetricsSystemExecutor

__version__ ="0.1.3"

def create_parser():
    parser = argparse.ArgumentParser(
        description="mozaiko CLI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--verbose", action="store_true")

    subparsers = parser.add_subparsers(dest="command")

    # -------------------------
    # PREPROCESS
    # -------------------------
    preprocess = subparsers.add_parser("preprocess", help="Pre-process FASTA database")
    preprocess.add_argument("-i", "--input", required=True)
    preprocess.add_argument("-o", "--output", required=True)
    preprocess.add_argument("--harmonized", action="store_true")
    preprocess.set_defaults(func=database_pre_process)

    # -------------------------
    # ASSIGN TAX
    # -------------------------
    assign = subparsers.add_parser("assign-tax", help="Assign taxonomy using CRABS")
    assign.add_argument("--json_file", required=True)
    assign.set_defaults(func=handle_taxonomic_assignment)

    # -------------------------
    # DEREPLICATE
    # -------------------------
    derep = subparsers.add_parser("dereplicate", help="Dereplicate sequences using CRABS")
    derep.add_argument("--json_file", required=True)
    derep.set_defaults(func=handle_dereplication)

    # -------------------------
    # IN SILICO
    # -------------------------
    insilico = subparsers.add_parser("insilico", help="Run in-silico PCR")
    insilico.add_argument("-i", "--input", required=True)
    insilico.add_argument("--run_name", required=True)
    insilico.add_argument("--primer_table", required=True)
    insilico.add_argument("--minimum_percentage_identity", type=float, default=0.5)
    insilico.set_defaults(func=handle_in_silico_analysis)

    # -------------------------
    # EVALUATE MULTIPLE
    # -------------------------
    eval_multi = subparsers.add_parser("evaluate-multi", help="Evaluate multiple OTLs")
    eval_multi.add_argument("--otl_folder", required=True)
    eval_multi.add_argument("--output_folder", required=True)
    eval_multi.add_argument("--primer_table", required=True)
    eval_multi.add_argument("--thresholds", nargs="+", type=float, default=[10.0, 5.0, 2.0])
    eval_multi.add_argument("--ranking_mode", default="flat")
    eval_multi.add_argument("--run_catnip", action="store_true")
    eval_multi.add_argument("--save_intermediate_ranks", action="store_true")
    eval_multi.set_defaults(func=handle_evaluate_multiple_otls)

    # -------------------------
    # EVALUATE SINGLE
    # -------------------------
    eval_single = subparsers.add_parser("evaluate-single", help="Evaluate single OTL")
    eval_single.add_argument("--otl_path", required=True)
    eval_single.add_argument("--output_folder", required=True)
    eval_single.add_argument("--primer_table", required=True)
    eval_single.add_argument("--thresholds", nargs="+", type=float, default=[10.0, 5.0, 2.0])
    eval_single.add_argument("--ranking_mode", default="flat")
    eval_single.add_argument("--run_catnip", action="store_true")
    eval_single.add_argument("--save_intermediate_ranks", action="store_true")
    eval_single.set_defaults(func=handle_evaluate_single_otl)

    # -------------------------
    # PIPELINE RUN
    # -------------------------

    run = subparsers.add_parser("run-pipeline", help="Run full pipeline from JSON config file.")
    run.add_argument("--config", required=True, help="Path to JSON configuration file.")
    run.set_defaults(func=handle_pipeline_run)

    return parser

def database_pre_process(args):
    print("mozaiko INFO: Pre-processing database...")

    try:
        fasta = CustomFastaImport()
        fasta.read_fasta(args.input)

        if args.harmonized:
            fasta.pre_process_harmonized_fasta_database()
            fasta.df2csv(args.output)
            print("mozaiko INFO: Database pre-processing completed successfully. File saved as .tsv and _processed.fasta. The latter can be used for in-silico analysis.")
        else:
            print("mozaiko INFO: Database pre-processing is not available for non-harmonized databases.")
            print("Please harmonize your database against the GBIF taxonomic backbone:")
            print("https://github.com/CIBIO-BU/DNAquaIMG/blob/main/tax_harmonize.py")
            print("Alternatively, include the following taxonomic hierarchy columns:")
            print("['scientificName', 'rank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']")
            print("mozaiko INFO: If your database is already harmonized, please use the --harmonized flag to pre-process it.")
            # fasta.df2csv(args.output)
            # print("mozaiko INFO: Database pre-processing completed successfully. File saved as .tsv and original fasta file overwritten with clean headers.")

    except Exception as e:
        logging.error(f"mozaiko ERROR: Failed to process the FASTA file: {e}")
        raise

def handle_taxonomic_assignment(args):
    print("mozaiko INFO: Assigning taxonomy (CRABS)...")
    CrabsScriptGenerator().run_assign_tax_command(args.json_file)


def handle_dereplication(args):
    print("mozaiko INFO: Dereplicating sequences (CRABS)...")
    CrabsScriptGenerator().run_dereplicate_command(args.json_file)


def handle_in_silico_analysis(args):
    print("mozaiko INFO: Running in-silico PCR...")
    insil = InSilicoAmplification(
        database_fasta_file=args.input,
        run_name=args.run_name
    )
    insil.run_in_silico_analysis(primer_table=args.primer_table,
                                 minimum_percentage_identity=args.minimum_percentage_identity)

def handle_evaluate_multiple_otls(args):
    print("mozaiko INFO: Evaluating multiple OTLs...")

    MetricsSystemExecutor.evaluate_several_OTLs(
        otl_folder=args.otl_folder,
        output_folder=args.output_folder,
        primer_table=args.primer_table,
        save_intermediate_ranks=args.save_intermediate_ranks,
        run_catnip=args.run_catnip,
        thresholds=args.thresholds,
        ranking_mode=args.ranking_mode
    )


def handle_evaluate_single_otl(args):
    print("mozaiko INFO: Evaluating single OTL...")

    MetricsSystemExecutor.evaluate_single_OTL(
        otl_path=args.otl_path,
        output_folder=args.output_folder,
        primer_table=args.primer_table,
        save_intermediate_ranks=args.save_intermediate_ranks,
        run_catnip=args.run_catnip,
        thresholds=args.thresholds,
        ranking_mode=args.ranking_mode
    )

import json
import os

def handle_pipeline_run(args):
    with open(args.config) as f:
        config = json.load(f)

    paths = config["paths"]
    steps = config["steps"]
    run_name = config["run_name"]

    output_root = paths["output_root"]
    os.makedirs(output_root, exist_ok=True)
    run_folder=os.path.join(output_root, run_name)
    os.makedirs(run_folder, exist_ok=True)

    # -------------------------
    # PREPROCESS
    # -------------------------
    print("mozaiko INFO: Pre-processing database...")
    if steps["preprocess"]["enabled"]:
        fasta = CustomFastaImport()
        fasta.read_fasta(paths["input_fasta"])

        if steps["preprocess"]["harmonized"]:
            processed_fasta = fasta.pre_process_harmonized_fasta_database()
        else:
            print("mozaiko INFO: Database pre-processing is not available for non-harmonized databases.")
            print("Please harmonize your database against the GBIF taxonomic backbone:")
            print("https://github.com/CIBIO-BU/DNAquaIMG/blob/main/tax_harmonize.py")
            print("Alternatively, include the following taxonomic hierarchy columns:")
            print("['scientificName', 'rank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']")
            print("mozaiko INFO: If your database is already harmonized, please use the --harmonized flag to pre-process it.")

        fasta.df2csv()

    print("mozaiko INFO: Database pre-processing completed successfully.")

    # -------------------------
    # IN SILICO
    # -------------------------
    print("mozaiko INFO: Running in-silico PCR...")
    if steps["insilico"]["enabled"]:
        insil = InSilicoAmplification(
            database_fasta_file=processed_fasta,
            run_name=run_name
        )

        insil.run_in_silico_analysis(
            primer_table=paths["primer_table"],
            minimum_percentage_identity=steps["insilico"].get(
                "minimum_percentage_identity", 0.5
            )
        )

    # -------------------------
    # EVALUATION
    # -------------------------
    if steps["evaluate_multiple_otl"]["enabled"]:
        MetricsSystemExecutor.evaluate_several_OTLs(
            otl_folder=paths["otl_folder"],
            output_folder=run_folder,
            primer_table=paths["primer_table"],
            save_intermediate_ranks=steps["evaluate_multiple_otl"]["save_intermediate_ranks"],
            run_catnip=steps["evaluate_multiple_otl"]["run_catnip"],
            thresholds=steps["evaluate_multiple_otl"]["thresholds"],
            ranking_mode=steps["evaluate_multiple_otl"]["ranking_mode"]
        )

def main():
    parser = create_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not hasattr(args, "func"):
        parser.print_help()
        return

    args.func(args)


if __name__ == "__main__":
    main()