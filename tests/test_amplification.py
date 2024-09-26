"""
Unit tests for the InSilicoAmplification class.
"""

import subprocess
import sys
import unittest
from io import StringIO
from unittest.mock import patch

import pandas as pd

from src.in_silico_analysis.amplification import InSilicoAmplification


class TestInSilicoAmplification(unittest.TestCase):
    """
    Class to test the InSilicoAmplification class.
    """

    def setUp(self):
        """
        Set up the test class and data.
        """
        self.data_dir = "data/test_data"
        self.primer_list = self.data_dir + "/test_primer_table.tsv"
        self.input_data = self.data_dir + "/fasta_example_file_taxid.fasta"
        self.amplification = InSilicoAmplification(self.input_data)

    @patch("subprocess.run")
    def test_check_if_cutadapt_installed(self, mock_subprocess):
        """
        Test that the _check_if_cutadapt_installed method runs.
        """
        self.amplification._check_if_cutadapt_installed()
        mock_subprocess.assert_called_with(
            ["cutadapt", "--version"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    def test_validate_primer_table_not_exist(self):
        """
        Test that the _validate_primer_table method raises a SystemExit when the primer table does
        not exist.
        """
        with self.assertRaises(SystemExit):
            self.amplification._validate_primer_table("nonexistent_file.tsv")

    def test_validate_primer_table_valid(self):
        """
        Test that the _validate_primer_table method validates the primer table correctly.
        """
        self.amplification._validate_primer_table(self.primer_list)

    def test_validate_primer_table_wrong_extension(self):
        """
        Test that the _validate_primer_table method raises a SystemExit when the primer table has
        the wrong extension.
        """
        with self.assertRaises(SystemExit):
            self.amplification._validate_primer_table(self.input_data)

    def test_read_primer_tables_valid(self):
        """
        Test that the _validate_primer_table method reads the primer table correctly.
        """
        sys.stdin = StringIO(f"{self.primer_list}\n")

        self.amplification.read_primer_tables()

        self.assertIsNotNone(self.amplification.primer_table)
        self.assertGreater(len(self.amplification.primer_table), 0)
        self.assertEqual(
            self.amplification.primer_table.loc[1]["target_group"], "Vertebrate"
        )
        self.assertEqual(
            self.amplification.primer_table.loc[0]["adapter"],
            "ACACCGCCCGTCACTCTC...GAGGAGGCAAGTCGTAACATG",
        )

    def test_validate_primer_table_wrong_columns(self):
        """
        Test that the _validate_primer_table method raises a SystemExit when the primer table does
        not contain the required columns.
        """
        self.amplification.primer_table_columns = ["x", "y"]

        user_input = StringIO(f"{self.primer_list}\n")
        captured_output = StringIO()
        sys.stdin, sys.stdout = user_input, captured_output

        with self.assertRaises(SystemExit):
            self.amplification.read_primer_tables()

        output = captured_output.getvalue()

        expected_error = "mozaiko INFO: The primer table must contain the following fields: ['x', 'y']"
        self.assertIn((expected_error), output)

        sys.stdin, sys.stdout = sys.__stdin__, sys.__stdout__

    def test_validate_fasta_not_exist(self):
        """
        Test that the _validate_fasta method raises a SystemExit when the input data does not exist.
        """
        test_class = InSilicoAmplification(data="nonexistent_fasta.fasta")
        with self.assertRaises(SystemExit):
            test_class._validate_fasta()

    def test_validate_fasta_valid(self):
        """
        Test that the _validate_fasta method validates the input data correctly.
        """
        test_class = InSilicoAmplification(self.input_data)
        test_class._validate_fasta()

    def test_validate_fasta_wrong_extension(self):
        """
        Test that the _validate_fasta method raises a SystemExit when the input data has the wrong
        extension.
        """
        test_class = InSilicoAmplification(self.primer_list)
        with self.assertRaises(SystemExit):
            test_class._validate_fasta()