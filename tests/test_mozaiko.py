"""
Unit tests for the mozaiko.py module.
"""

import argparse
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from src.mozaiko.mozaiko import (
    create_parser,
    database_pre_process,
    handle_dereplication,
    handle_taxonomic_assignment,
    main,
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