"""
Unit tests for the mozaiko.py module.
"""

import argparse
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch, call

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from src.mozaiko.mozaiko import (
    create_parser,
    database_pre_process,
    handle_dereplication,
    handle_taxonomic_assignment,
    main,
    handle_evaluate_multiple_otls,
    handle_evaluate_single_otl,
    handle_in_silico_analysis,
    handle_pipeline_run,
)


class Testmozaiko(unittest.TestCase):
    """
    Class to test the mozaiko.py module.
    """

    def test_create_parser(self):
        """
        Test that the create_parser function returns an ArgumentParser object with the correct
        description.
        """
        parser = create_parser()
        self.assertIsInstance(parser, argparse.ArgumentParser)
        self.assertEqual(
            parser.description,
            "mozaiko CLI",
        )

    @patch("src.mozaiko.mozaiko.CustomFastaImport")
    def test_database_pre_process(self, mock_fasta_import):
        """
        Test database preprocessing for harmonized databases.
        """
        mock_fasta = MagicMock()
        mock_fasta_import.return_value = mock_fasta

        args = argparse.Namespace(
            input="data/test_data/fasta_example_file.fasta",
            output="output.tsv",
            harmonized=True,
        )

        database_pre_process(args)

        mock_fasta.read_fasta.assert_called_once_with(
            "data/test_data/fasta_example_file.fasta"
        )
        mock_fasta.pre_process_harmonized_fasta_database.assert_called_once()
        mock_fasta.df2csv.assert_called_once_with("output.tsv")

        if Path("output.tsv").exists():
            Path("output.tsv").unlink()

    @patch("src.mozaiko.mozaiko.CrabsScriptGenerator.run_assign_tax_command")
    def test_handle_taxonomic_assignment(self, mock_run_assign_tax_command):
        """
        Test to check if the handle_taxonomic_assignment function runs the assign_tax command with
        the correct parameters.
        """
        args = argparse.Namespace(json_file="dummy.json")
        handle_taxonomic_assignment(args)
        mock_run_assign_tax_command.assert_called_with("dummy.json")

        # remove 'taxonomy_files' folder if it was created
        taxonomy_folder = Path("taxonomy_files")
        if taxonomy_folder.exists() and taxonomy_folder.is_dir():
            for file in taxonomy_folder.iterdir():
                file.unlink()
            taxonomy_folder.rmdir()

    @patch("src.mozaiko.mozaiko.CrabsScriptGenerator.run_dereplicate_command")
    def test_handle_dereplication(self, mock_run_dereplicate_command):
        """
        Test to check if the handle_dereplication function runs the dereplication command with the
        correct parameters.
        """
        args = argparse.Namespace(json_file="dummy.json")
        handle_dereplication(args)
        mock_run_dereplicate_command.assert_called_with("dummy.json")

    @patch("src.mozaiko.mozaiko.create_parser")
    def test_main_without_subcommand(self, mock_create_parser):
        """
        Test main() when no subcommand is provided.
        """
        mock_parser = MagicMock()
        mock_args = MagicMock()

        del mock_args.func

        mock_parser.parse_args.return_value = mock_args
        mock_create_parser.return_value = mock_parser

        main()

        mock_parser.print_help.assert_called_once()

    @patch("src.mozaiko.mozaiko.create_parser")
    def test_main_dispatches_function(self, mock_create_parser):
        """
        Test that main dispatches to the selected subcommand function.
        """
        mock_parser = MagicMock()

        mock_func = MagicMock()

        mock_args = argparse.Namespace(
            verbose=False,
            func=mock_func,
        )

        mock_parser.parse_args.return_value = mock_args
        mock_create_parser.return_value = mock_parser

        main()

        mock_func.assert_called_once_with(mock_args)

    @patch("src.mozaiko.mozaiko.create_parser")
    def test_main_sets_verbose_logging(self, mock_create_parser):
        """
        Test verbose logging activation.
        """
        mock_parser = MagicMock()

        mock_func = MagicMock()

        mock_args = argparse.Namespace(
            verbose=True,
            func=mock_func,
        )

        mock_parser.parse_args.return_value = mock_args
        mock_create_parser.return_value = mock_parser

        with patch("logging.getLogger") as mock_get_logger:
            main()

            mock_get_logger.return_value.setLevel.assert_called_once()

    @patch("builtins.print")
    @patch("src.mozaiko.mozaiko.CustomFastaImport")
    def test_database_pre_process_non_harmonized(
        self,
        mock_fasta_import,
        mock_print,
    ):
        """
        Test database preprocessing for non-harmonized databases.
        """
        mock_fasta = MagicMock()
        mock_fasta_import.return_value = mock_fasta

        args = argparse.Namespace(
            input="data/test_data/fasta_example_file.fasta",
            output="output.tsv",
            harmonized=False,
        )

        database_pre_process(args)

        mock_fasta.read_fasta.assert_called_once_with(
            "data/test_data/fasta_example_file.fasta"
        )

        mock_fasta.pre_process_harmonized_fasta_database.assert_not_called()
        mock_fasta.df2csv.assert_not_called()

        mock_print.assert_any_call(
            "mozaiko INFO: Database pre-processing is not available for non-harmonized databases."
        )

        mock_print.assert_any_call(
            "Please harmonize your database against the GBIF taxonomic backbone:"
        )

    @patch("src.mozaiko.mozaiko.logging.error")
    @patch("src.mozaiko.mozaiko.CustomFastaImport")
    def test_database_pre_process_exception(
        self,
        mock_fasta_import,
        mock_logging_error,
    ):
        """
        Test exception handling during FASTA preprocessing.
        """
        mock_fasta = MagicMock()
        mock_fasta.read_fasta.side_effect = Exception("FASTA read failed")

        mock_fasta_import.return_value = mock_fasta

        args = argparse.Namespace(
            input="bad.fasta",
            output="output.tsv",
            harmonized=True,
        )

        with self.assertRaises(Exception):
            database_pre_process(args)

        mock_logging_error.assert_called_once_with(
            "mozaiko ERROR: Failed to process the FASTA file: FASTA read failed"
        )

    @patch("src.mozaiko.mozaiko.InSilicoAmplification")
    def test_handle_in_silico_analysis(self, mock_insilico):
        """
        Test in-silico PCR execution.
        """
        mock_instance = MagicMock()
        mock_insilico.return_value = mock_instance

        args = argparse.Namespace(
            input="database.fasta",
            run_name="test_run",
            primer_table="primers.tsv",
            minimum_percentage_identity=0.75,
        )

        handle_in_silico_analysis(args)

        mock_insilico.assert_called_once_with(
            database_fasta_file="database.fasta",
            run_name="test_run",
        )

        mock_instance.run_in_silico_analysis.assert_called_once_with(
            primer_table="primers.tsv",
            minimum_percentage_identity=0.75,
        )

    @patch(
        "src.mozaiko.mozaiko.MetricsSystemExecutor.evaluate_several_OTLs"
    )
    def test_handle_evaluate_multiple_otls(
        self,
        mock_evaluate,
    ):
        """
        Test evaluation of multiple OTLs.
        """
        args = argparse.Namespace(
            otl_folder="otls/",
            output_folder="output/",
            primer_table="primers.tsv",
            save_intermediate_ranks=True,
            run_catnip=False,
            thresholds=[10.0, 5.0],
            ranking_mode="flat",
        )

        handle_evaluate_multiple_otls(args)

        mock_evaluate.assert_called_once_with(
            otl_folder="otls/",
            output_folder="output/",
            primer_table="primers.tsv",
            save_intermediate_ranks=True,
            run_catnip=False,
            thresholds=[10.0, 5.0],
            ranking_mode="flat",
        )

    @patch(
    "src.mozaiko.mozaiko.MetricsSystemExecutor.evaluate_single_OTL"
    )
    def test_handle_evaluate_single_otl(
        self,
        mock_evaluate,
    ):
        """
        Test evaluation of a single OTL.
        """
        args = argparse.Namespace(
            otl_path="otl.tsv",
            output_folder="output/",
            primer_table="primers.tsv",
            save_intermediate_ranks=False,
            run_catnip=True,
            thresholds=[10.0, 5.0],
            ranking_mode="weighted",
        )

        handle_evaluate_single_otl(args)

        mock_evaluate.assert_called_once_with(
            otl_path="otl.tsv",
            output_folder="output/",
            primer_table="primers.tsv",
            save_intermediate_ranks=False,
            run_catnip=True,
            thresholds=[10.0, 5.0],
            ranking_mode="weighted",
        )


    @patch("src.mozaiko.mozaiko.os.makedirs")
    @patch("src.mozaiko.mozaiko.MetricsSystemExecutor.evaluate_several_OTLs")
    @patch("src.mozaiko.mozaiko.InSilicoAmplification")
    @patch("src.mozaiko.mozaiko.CustomFastaImport")
    @patch("builtins.open", new_callable=mock_open)
    @patch("src.mozaiko.mozaiko.json.load")
    def test_handle_pipeline_run(
        self,
        mock_json_load,
        mock_open_file,
        mock_fasta_import,
        mock_insilico,
        mock_evaluate,
        mock_makedirs,
    ):
        """
        Test full pipeline execution from config.
        """

        mock_json_load.return_value = {
            "run_name": "test_run",
            "paths": {
                "output_root": "output",
                "input_fasta": "database.fasta",
                "primer_table": "primers.tsv",
                "otl_folder": "otls/",
            },
            "steps": {
                "preprocess": {
                    "enabled": True,
                    "harmonized": True,
                },
                "insilico": {
                    "enabled": True,
                    "minimum_percentage_identity": 0.8,
                },
                "evaluate_multiple_otl": {
                    "enabled": True,
                    "save_intermediate_ranks": True,
                    "run_catnip": False,
                    "thresholds": [10.0, 5.0],
                    "ranking_mode": "flat",
                },
            },
        }

        mock_fasta = MagicMock()
        mock_fasta.pre_process_harmonized_fasta_database.return_value = (
            "processed.fasta"
        )

        mock_fasta_import.return_value = mock_fasta

        mock_insilico_instance = MagicMock()
        mock_insilico.return_value = mock_insilico_instance

        args = argparse.Namespace(config="config.json")

        handle_pipeline_run(args)

        mock_open_file.assert_called_once_with("config.json")

        mock_fasta.read_fasta.assert_called_once_with(
            "database.fasta"
        )

        mock_fasta.pre_process_harmonized_fasta_database.assert_called_once()

        mock_insilico.assert_called_once_with(
            database_fasta_file="processed.fasta",
            run_name="test_run",
        )

        mock_insilico_instance.run_in_silico_analysis.assert_called_once_with(
            primer_table="primers.tsv",
            minimum_percentage_identity=0.8,
        )

        mock_evaluate.assert_called_once_with(
            otl_folder="otls/",
            output_folder="output/test_run",
            primer_table="primers.tsv",
            save_intermediate_ranks=True,
            run_catnip=False,
            thresholds=[10.0, 5.0],
            ranking_mode="flat",
        )

        mock_makedirs.assert_has_calls(
            [
                call("output", exist_ok=True),
                call("output/test_run", exist_ok=True),
            ]
        )

    @patch("builtins.print")
    @patch("src.mozaiko.mozaiko.CustomFastaImport")
    @patch("builtins.open", new_callable=mock_open)
    @patch("src.mozaiko.mozaiko.json.load")
    def test_handle_pipeline_run_non_harmonized(
        self,
        mock_json_load,
        mock_open_file,
        mock_fasta_import,
        mock_print,
    ):
        """
        Test pipeline with non-harmonized preprocessing.
        """

        mock_json_load.return_value = {
            "run_name": "test_run",
            "paths": {
                "output_root": "output",
                "input_fasta": "database.fasta",
                "primer_table": "primers.tsv",
                "otl_folder": "otls/",
            },
            "steps": {
                "preprocess": {
                    "enabled": True,
                    "harmonized": False,
                },
                "insilico": {
                    "enabled": False,
                },
                "evaluate_multiple_otl": {
                    "enabled": False,
                },
            },
        }

        mock_fasta = MagicMock()
        mock_fasta_import.return_value = mock_fasta

        args = argparse.Namespace(config="config.json")

        handle_pipeline_run(args)

        mock_print.assert_any_call(
            "mozaiko INFO: Database pre-processing is not available for non-harmonized databases."
        )