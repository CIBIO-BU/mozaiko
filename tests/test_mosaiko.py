"""
Unit tests for the mosaiko.py module.
"""
import argparse
import unittest
from unittest.mock import mock_open, patch

from src.mosaiko import (
    create_parser,
    handle_custom_fasta_import,
    handle_dereplication,
    handle_taxonomic_assignment,
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
