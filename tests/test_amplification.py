"""
Unit tests for the InSilicoAmplification class.
"""

import os
import subprocess
import sys
import unittest
from io import StringIO
from pathlib import Path, PosixPath
from unittest.mock import MagicMock, mock_open, patch

import pandas as pd
from Bio.Seq import Seq

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
            with patch("sys.stdout", new=StringIO()) as mock_stdout:
                self.amplification._check_if_cutadapt_installed()
            self.assertEqual(context.exception.code, 1)  # 1 -> SystemExist
            message = "mozaiko INFO: Cutadapt is not installed. Please install Cutadapt before running this script. \n Cutadapt can be found at https://cutadapt.readthedocs.io/en/stable/installation.html"
            self.assertEqual(message, mock_stdout.get_value())

    def test_validate_primer_table_not_exist(self):
        """
        Test that the _validate_primer_table method raises a SystemExit when the primer table does
        not exist.
        """
        with self.assertRaises(SystemExit):
            self.amplification.validate_primer_table("nonexistent_file.tsv")

    def test_validate_primer_table_valid(self):
        """
        Test that the _validate_primer_table method validates the primer table correctly.
        """
        self.amplification.validate_primer_table(self.primer_list)

    def test_validate_primer_table_wrong_extension(self):
        """
        Test that the _validate_primer_table method raises a SystemExit when the primer table has
        the wrong extension.
        """
        with self.assertRaises(SystemExit):
            self.amplification.validate_primer_table(self.input_data)

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

        with patch("sys.stdin", user_input), patch("sys.stdout", captured_output):
            with self.assertRaises(SystemExit):
                self.amplification.read_primer_tables()

            output = captured_output.getvalue()

            self.assertIn(
                "mozaiko INFO: The primer table is missing the following required fields:",
                output,
            )
            self.assertIn("x", output)
            self.assertIn("y", output)
            self.assertIn("Required fields are: x, y", output)

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

    @patch("shutil.copytree")
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
    @patch("builtins.input", side_effect=["test_output_folder"])
    @patch("pathlib.Path.is_dir", return_value=True)
    @patch("pathlib.Path.glob", return_value=[Path("dummy.fasta")])
    @patch(
        "builtins.open", new_callable=mock_open, read_data=">seq1\nATGC\n>seq2\nATGC"
    )
    @patch("Bio.SeqIO.parse")
    def test_run_in_silico_analysis_calls(
        self,
        _mock_parse,
        _mock_open,
        _mock_glob,
        _mock_dir,
        _mock_input,
        mock_read_tables,
        mock_validate_fasta,
        mock_check_crabs,
        mock_check_cutadapt,
        _mock_copytree,
    ):
        """
        Test that run_in_silico_analysis calls all required methods.
        """
        self.amplification.primer_table = MagicMock()

        self.amplification.run_in_silico_analysis()

        mock_check_cutadapt.assert_called_once()
        mock_check_crabs.assert_called_once()
        mock_validate_fasta.assert_called_once()
        mock_read_tables.assert_called_once()
        self.assertEqual(self.amplification.run_name, "test_output_folder")
        self.assertIsNotNone(self.amplification.output_dirs)

    @patch("shutil.copytree")
    @patch("builtins.input", side_effect=["test_folder"])
    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.process_commands"
    )
    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.read_primer_tables"
    )
    @patch("pathlib.Path.is_dir", return_value=True)
    @patch("pathlib.Path.glob", return_value=[Path("dummy.fasta")])
    @patch("os.path.exists", return_value=True)
    @patch(
        "builtins.open", new_callable=mock_open, read_data=">seq1\nATGC\n>seq2\nATGC"
    )
    @patch("Bio.SeqIO.parse")
    def test_run_in_silico_analysis_calls_process_commands(
        self,
        _mock_parse,
        _mock_open,
        _mock_exists,
        _mock_glob,
        _mock_is_dir,
        mock_read_tables,
        mock_process_commands,
        _mock_input,
        _mock_copytree,
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
                    "assay_name": "Chon01",
                    "fwd_seq": "ACACCGCCCGTCACTCTC",
                    "rev_seq": "CATGTTACGACTTGCCTCCTC",
                },
            ),
            (
                1,
                {
                    "target_group": "Vertebrate",
                    "assay_name": "12S-V5-c",
                    "fwd_seq": "AGGGATAACAGCGCAATC",
                    "rev_seq": "TCGTTGAACAAACGAACC",
                },
            ),
        ]
        mock_read_tables.return_value = None

        self.amplification.run_in_silico_analysis()

        self.assertEqual(mock_process_commands.call_count, 2)

        mock_process_commands.assert_any_call(
            {
                "target_group": "Chondrichthyes",
                "assay_name": "Chon01",
                "fwd_seq": "ACACCGCCCGTCACTCTC",
                "rev_seq": "CATGTTACGACTTGCCTCCTC",
            },
            self.amplification.data,
        )
        mock_process_commands.assert_any_call(
            {
                "target_group": "Vertebrate",
                "assay_name": "12S-V5-c",
                "fwd_seq": "AGGGATAACAGCGCAATC",
                "rev_seq": "TCGTTGAACAAACGAACC",
            },
            self.amplification.data,
        )

    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.run_cutadapt_command"
    )
    @patch("src.in_silico_analysis.amplification.InSilicoAmplification.run_pga_command")
    def test_process_commands(self, mock_run_pga_command, mock_run_cutadapt_command):
        row = {
            "target_group": "Chondrichthyes",
            "barcode_region": "12S",
            "assay_name": "Chon01",
            "fwd_seq": "ACACCGCCCGTCACTCTC",
            "correct_reverse_primer": "CATGTTACGACTTGCCTCCTC",
            "adapter": "ACACCGCCCGTCACTCTC...GAGGAGGCAAGTCGTAACATG",
            "min_read_length": 60,
            "max_read_length": 600,
            "overlap": 120,
        }

        self.amplification.run_name = "test_run"
        self.amplification.output_dirs = self.amplification._setup_output_directories(
            "test_run"
        )

        self.amplification.process_commands(row, self.input_data)

        self.assertEqual(mock_run_cutadapt_command.call_count, 2)
        self.assertEqual(mock_run_pga_command.call_count, 2)

    @patch("subprocess.run")
    @patch("pathlib.Path.stat")
    def test_run_cutadapt_command(self, mock_stat, mock_subprocess_run):
        output_dir = Path("output")

        # mock stat to return non-zero sized files
        mock_stat.return_value = os.stat_result((0, 0, 0, 0, 0, 0, 100, 0, 0, 0))

        common_args = [
            "cutadapt",
            "-g",
            "ADAPTER",
            "--output",
            str(output_dir / "12S_Chon01.fasta"),
            str(self.input_data),
            "--no-indels",
            "-e",
            "3",
            "--overlap",
            "20",
            "--revcomp",
            "--quiet",
            "--minimum-length",
            "30",
            "--maximum-length",
            "100",
            "--discard-untrimmed",
        ]

        # Test case 1: 'amplicon' command type
        mock_subprocess_run.reset_mock()  # reset state to complete more tests
        self.amplification.run_cutadapt_command(
            "amplicon",
            "ADAPTER",
            self.input_data,
            20,
            30,
            100,
            "12S",
            "Chon01",
            output_dir,
        )
        expected_args = common_args + ["--action", "retain"]
        mock_subprocess_run.assert_called_with(
            expected_args, check=True, capture_output=True, text=True, encoding="utf-8"
        )

        # Test case 3: 'insert' command type
        mock_subprocess_run.reset_mock()
        self.amplification.run_cutadapt_command(
            "insert",
            "ADAPTER",
            self.input_data,
            20,
            30,
            100,
            "12S",
            "Chon01",
            output_dir,
        )
        expected_args = common_args + ["--action", "trim"]
        mock_subprocess_run.assert_called_with(
            expected_args, check=True, capture_output=True, text=True, encoding="utf-8"
        )

        # Test case 4: Invalid command type
        with self.assertRaises(ValueError):
            self.amplification.run_cutadapt_command(
                "invalid_type",
                "ADAPTER",
                self.input_data,
                20,
                30,
                100,
                "12S",
                "Chon01",
                output_dir,
            )

        # Test case 5: Error handling (subprocess.CalledProcessError)
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, "cutadapt")
        with self.assertRaises(subprocess.CalledProcessError):
            self.amplification.run_cutadapt_command(
                "amplicon",
                "ADAPTER",
                self.input_data,
                20,
                30,
                100,
                "12S",
                "Chon01",
                output_dir,
            )

        # Test case 6: Error handling (FileNotFoundError)
        mock_subprocess_run.side_effect = FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            self.amplification.run_cutadapt_command(
                "amplicon",
                "ADAPTER",
                self.input_data,
                20,
                30,
                100,
                "12S",
                "Chon01",
                output_dir,
            )

        # Test case 7: Empty output file
        mock_stat.return_value = os.stat_result((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        mock_subprocess_run.side_effect = None
        self.amplification.run_cutadapt_command(
            "amplicon",
            "ADAPTER",
            self.input_data,
            20,
            30,
            100,
            "12S",
            "Chon01",
            output_dir,
        )

    @patch("subprocess.run")
    @patch("pathlib.Path.stat")
    def test_run_pga_command(self, mock_stat, mock_subprocess_run):
        output_dir = Path("output")
        database_dir = Path("database")

        mock_stat.return_value = os.stat_result((0, 0, 0, 0, 0, 0, 100, 0, 0, 0))

        # Test case 1: Successful run
        mock_subprocess_run.return_value = MagicMock(returncode=0)

        self.amplification.run_pga_command(
            self.input_data,
            "FORWARD",
            "REVERSE",
            "12S",
            "Assay1",
            output_dir,
            database_dir,
            "relaxed",
        )
        mock_subprocess_run.assert_called_with(
            [
                "crabs",
                "pga",
                "--input",
                str(self.input_data),
                "--output",
                str(output_dir / "12S_Assay1.fasta"),
                "--database",
                str(database_dir / "12S_Assay1.fasta"),
                "--fwd",
                "FORWARD",
                "--rev",
                "REVERSE",
                "--speed",
                "slow",
                "--percid",
                "0.75",
                "--coverage",
                "99",
                "--filter_method",
                "relaxed",
            ],
            check=True,
            stdout=-3,
            stderr=-2,
        )

        # Test case 2: Failed run
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(1, "crabs")
        with self.assertRaises(SystemExit):
            self.amplification.run_pga_command(
                self.input_data,
                "FORWARD",
                "REVERSE",
                "12S",
                "Assay1",
                output_dir,
                database_dir,
                "strict",
            )

    def test_read_primer_tables_max_len_true(self):
        self.amplification.read_primer_tables(
            self.primer_list, max_len_according_to_ilumina=True
        )
        self.assertIsNotNone(self.amplification.primer_table)

        for index, row in self.amplification.primer_table.iterrows():
            foward_primer = row["fwd_seq"]
            reverse_primer = row["rev_seq"]
            correct_reverse_primer = str(Seq(reverse_primer).reverse_complement())
            forward_primer_length = len(foward_primer)
            correct_reverse_primer_length = len(correct_reverse_primer)
            expected_max_length = (
                600 - correct_reverse_primer_length - forward_primer_length
            )
            self.assertEqual(row["max_read_length"], expected_max_length)

    def test_read_primer_tables_max_len_false(self):
        self.amplification.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fwd_seq",
            "rev_seq",
            "min_read_length",
            "max_read_length",
        ]
        self.amplification.read_primer_tables(
            self.primer_list, max_len_according_to_ilumina=False
        )
        self.assertIsNotNone(self.amplification.primer_table)
        self.assertIn("min_read_length", self.amplification.primer_table.columns)
        self.assertIn("max_read_length", self.amplification.primer_table.columns)

        test_primer_table = pd.read_csv(self.primer_list, sep="\t")

        for index, row in self.amplification.primer_table.iterrows():
            expected_min_read_length = test_primer_table.loc[index, "min_read_length"]
            expected_max_read_length = test_primer_table.loc[index, "max_read_length"]
            self.assertEqual(row["min_read_length"], expected_min_read_length)
            self.assertEqual(row["max_read_length"], expected_max_read_length)

    def test_read_primer_tables_max_len_false_missing_cols(self):
        self.amplification.primer_table_columns = [
            "target_group",
            "barcode_region",
            "assay_name",
            "fwd_seq",
            "rev_seq",
        ]
        with self.assertRaises(ValueError):
            self.amplification.read_primer_tables(
                self.primer_list, max_len_according_to_ilumina=False
            )

    def tearDown(self):
        if os.path.exists("dummy.fasta"):
            os.remove("dummy.fasta")
