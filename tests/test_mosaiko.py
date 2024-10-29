"""
Unit tests for the mosaiko.py module.
"""

import argparse
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from src.mozaiko import (
    create_parser,
    handle_custom_fasta_import,
    handle_dereplication,
    handle_taxonomic_assignment,
    main,
)


class TestMosaiko(unittest.TestCase):
    """
    Class to test the mosaiko.py module.
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
            "mosaiko: Piecing Together Complete Genetic Coverage for Biomonitoring",
        )

    @patch("builtins.open", mock_open(read_data=""))
    @patch("os.path.exists", return_value=True)
    @patch("os.path.getsize", return_value=100)
    @patch("builtins.input", return_value="exit")
    def test_handle_custom_fasta_import(self, _mock_exists, _mock_getsize, _mock_input):
        """
        Test that the handle_custom_fasta_import function reads a FASTA file and asks for a lineage
        file.
        """
        args = argparse.Namespace(input="input.fasta", output="output.fasta")
        handle_custom_fasta_import(args)

    @patch(
        "src.reference_database.db_curation.CrabsScriptGenerator.run_assign_tax_command"
    )
    def test_handle_taxonomic_assignment(self, mock_run_assign_tax_command):
        """
        Test to check if the handle_taxonomic_assignment function runs the assign_tax command with
        the correct parameters.
        """
        args = argparse.Namespace(json_file="dummy.json")
        handle_taxonomic_assignment(args)
        mock_run_assign_tax_command.assert_called_with("dummy.json")

    @patch(
        "src.reference_database.db_curation.CrabsScriptGenerator.run_dereplicate_command"
    )
    def test_handle_dereplication(self, mock_run_dereplicate_command):
        """
        Test to check if the handle_dereplication function runs the dereplication command with the
        correct parameters.
        """
        args = argparse.Namespace(json_file="dummy.json")
        handle_dereplication(args)
        mock_run_dereplicate_command.assert_called_with("dummy.json")

    @patch("argparse.ArgumentParser.parse_args")
    @patch("src.mosaiko.CustomFastaImport")
    def test_main_load_custom_fasta(self, mock_custom_fasta_import, mock_parse_args):
        mock_fasta_instance = MagicMock()
        mock_custom_fasta_import.return_value = mock_fasta_instance

        mock_args = MagicMock(
            load_custom_fasta=True,
            input="data/test_data/fasta_example_file_taxid.fasta",
            assign_tax=False,
            dereplicate=False,
            json_file=None,
            verbose=False,
            in_silico_analysis=False,
            output=None,
        )
        mock_parse_args.return_value = mock_args

        main()

        mock_custom_fasta_import.assert_called_once()
        mock_fasta_instance.read_fasta.assert_called_once_with(mock_args.input)
        mock_fasta_instance.get_number_of_sequences.assert_called_once()

    @patch("argparse.ArgumentParser.parse_args")
    @patch("mosaiko.handle_taxonomic_assignment")
    def test_main_assign_tax_without_json(
        self, mock_taxonomic_assignment, mock_parse_args
    ):
        mock_parse_args.return_value = MagicMock(
            load_custom_fasta=False,
            input=None,
            assign_tax=True,
            json_file=None,
            verbose=False,
        )

        with self.assertLogs(level="ERROR") as log:
            main()

        mock_taxonomic_assignment.assert_not_called()
        self.assertIn(
            "mosaiko INFO: No JSON file specified. Please specify a JSON file with parameter --json_file.",
            log.output[0],
        )
        self.assertIn("Exiting...", log.output[1])

    @patch("argparse.ArgumentParser.parse_args")
    @patch("mosaiko.handle_dereplication")
    def test_main_dereplicate_without_json(self, mock_dereplication, mock_parse_args):
        mock_parse_args.return_value = MagicMock(
            load_custom_fasta=False,
            input=None,
            assign_tax=False,
            dereplicate=True,
            json_file=None,
            verbose=False,
        )

        with self.assertLogs(level="ERROR") as log:
            main()

        mock_dereplication.assert_not_called()
        self.assertIn(
            "mosaiko INFO: No JSON file specified. Please specify a JSON file with parameter --json_file.",
            log.output[0],
        )
        self.assertIn("Exiting...", log.output[1])
