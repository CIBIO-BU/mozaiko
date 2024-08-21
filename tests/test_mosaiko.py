import argparse
import unittest
from unittest.mock import mock_open, patch

from src.mosaiko import (
    create_parser,
    handle_custom_fasta_import,
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
    def test_handle_custom_fasta_import(self, mock_exists, mock_getsize, mock_input):
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
        Test that the handle_taxonomic_assignment function runs the assign_tax command with the
        correct parameters.
        """
        args = argparse.Namespace(json_file="dummy.json")
        handle_taxonomic_assignment(args)
        mock_run_assign_tax_command.assert_called_with("dummy.json")

    @patch("src.mosaiko.handle_custom_fasta_import")
    @patch("src.mosaiko.handle_taxonomic_assignment")
    def test_main(
        self, mock_handle_taxonomic_assignment, mock_handle_custom_fasta_import
    ):
        """
        Test that the main function calls the handle_custom_fasta_import and
        handle_taxonomic_assignment functions with the correct arguments.
        """
        args = argparse.Namespace(
            input="input.fasta",
            output="output.fasta",
            load_custom_fasta=True,
            json_file="dummy.json",
            assign_tax=True,
            verbose=True,
        )
        with patch("argparse.ArgumentParser.parse_args", return_value=args):
            main()
            mock_handle_custom_fasta_import.assert_called_with(args)
            mock_handle_taxonomic_assignment.assert_called_with(args)
