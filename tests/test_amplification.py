"""
Unit tests for the InSilicoAmplification class.
"""

import os
import subprocess
import sys
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

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

    @patch("subprocess.run")
    def test_check_if_cutadapt_installed_not_found(self, mock_run):
        """
        Test if _check_if_cutadapt_installed method raises FileNotFoundError() when no cutadapt
        installation is found.
        """
        mock_run.side_effect = FileNotFoundError()
        with self.assertRaises(
            SystemExit
        ) as context:  # captures the exception into context
            # redirect sys.stdout to a StringIO object
            with patch("sys.stdout", new=StringIO()) as fake_out:
                self.amplification._check_if_cutadapt_installed()
            self.assertEqual(context.exception.code, 1)  # 1 -> SystemExist
            message = "mozaiko INFO: Cutadapt is not installed. Please install Cutadapt before running this script. \n Cutadapt can be found at https://cutadapt.readthedocs.io/en/stable/installation.html"
            self.assertEqual(message, fake_out.get_value())

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

    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification._check_if_cutadapt_installed"
    )
    @patch(
        "src.reference_database.db_curation.CrabsScriptGenerator.check_if_crabs_installed"
    )
    @patch("src.in_silico_analysis.amplification.InSilicoAmplification._validate_fasta")
    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.read_primer_tables"
    )
    @patch("src.in_silico_analysis.amplification.Path")
    @patch("builtins.input", side_effect=["path_to_table.tsv"])
    def test_run_in_silico_analysis_calls(
        self,
        mock_check_cutadapt,
        mock_check_crabs,
        mock_validate_fasta,
        mock_read_tables,
        mock_path,
        _mock_input,
    ):
        """
        Test that run_in_silico_analysis calls the all required methods.
        """

        mock_path.return_value = "test_output_folder"
        self.amplification.primer_table = MagicMock()

        self.amplification.run_in_silico_analysis()

        mock_check_cutadapt.assert_called_once()
        mock_check_crabs.assert_called_once()
        mock_validate_fasta.assert_called_once()
        mock_read_tables.assert_called_once()
        mock_path.assert_called_once()

    @patch("builtins.input", side_effect=["test_folder"])
    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.process_commands"
    )
    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.read_primer_tables"
    )
    def test_run_in_silico_analysis_calls_process_commands(
        self, mock_read_tables, mock_process_commands, _mock_input
    ):
        """
        Test that run_in_silico_analysis calls the process_commands the correct number of times and
        with the correct arguments.
        """
        self.amplification.primer_table = MagicMock()
        self.amplification.primer_table.iterrows.return_value = [
            (
                0,
                {
                    "target_group": "Chondrichthyes",
                    "primer_name": "Chon01",
                    "forward_sequence": "ACACCGCCCGTCACTCTC",
                    "reverse_sequence": "CATGTTACGACTTGCCTCCTC",
                    "amplicon_length": "43",
                    "expected_size": "23",
                },
            ),
            (
                1,
                {
                    "target_group": "Vertebrate",
                    "primer_name": "12S-V5-c",
                    "forward_sequence": "AGGGATAACAGCGCAATC",
                    "reverse_sequence": "TCGTTGAACAAACGAACC",
                    "amplicon_length": "74",
                    "expected_size": "24",
                },
            ),
        ]
        mock_read_tables.return_value = None

        self.amplification.run_in_silico_analysis()

        self.assertEqual(mock_process_commands.call_count, 2)

        mock_process_commands.assert_any_call(
            {
                "target_group": "Chondrichthyes",
                "primer_name": "Chon01",
                "forward_sequence": "ACACCGCCCGTCACTCTC",
                "reverse_sequence": "CATGTTACGACTTGCCTCCTC",
                "amplicon_length": "43",
                "expected_size": "23",
            },
            Path("test_folder"),
            self.amplification.data,
        )
        mock_process_commands.assert_any_call(
            {
                "target_group": "Vertebrate",
                "primer_name": "12S-V5-c",
                "forward_sequence": "AGGGATAACAGCGCAATC",
                "reverse_sequence": "TCGTTGAACAAACGAACC",
                "amplicon_length": "74",
                "expected_size": "24",
            },
            Path("test_folder"),
            self.amplification.data,
        )
