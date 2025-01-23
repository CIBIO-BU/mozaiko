"""
Unit tests for scoring_utils.py
"""

import shutil
import tempfile
import unittest
from io import StringIO
from pathlib import Path

import pandas as pd

from src.marker_scoring.scoring_utils import *


class TestScoringUtils(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path("data/test_data")

    def test_calculate_iupac_mismatches(self):
        """
        Function to test the calculate_iupac_mismatches function.
        """
        mismatches_cases = self.test_directory / "iupac_mismatches_cases.txt"

        with open(mismatches_cases, "r") as file:
            for line in file:
                seq1, seq2, expected_mismatches, gc_matches = line.strip().split(" ")
                mismatches = calculate_iupac_mismatches(seq1, seq2)
                self.assertEqual(mismatches, int(expected_mismatches))

    def test_calculate_gc_matches(self):
        """
        Function to test the search_gc_clamp functonality, within the calculate_iupac_mismatches function.
        """
        mismatches_cases = self.test_directory / "iupac_mismatches_cases.txt"

        with open(mismatches_cases, "r") as file:
            for line in file:
                print(line)
                seq1, seq2, expected_mismatches, expected_gc_matches = (
                    line.strip().split(" ")
                )
                mismatches, real_gc_matches = calculate_iupac_mismatches(
                    seq1, seq2, search_gc_clamp=True
                )
                self.assertEqual(real_gc_matches, int(expected_gc_matches))

    def test_calculate_ambiguous_percentage(self):
        """
        Function to test the calculate_ambiguous_percentage.
        """
        sequences = ["ACGTACGTAG", "ACGTANGATC", "RWNBVDNWSN"]  # 0%  # 10%  # 100%
        ambiguity = []

        for sequence in sequences:
            result = calculate_ambiguous_percentage(sequence)
            ambiguity.append(result)

        self.assertEqual(ambiguity[0], 0.0)

        self.assertEqual(ambiguity[1], 0.1)

        self.assertEqual(ambiguity[2], 1.0)

    def test_filter_sequences_by_ambiguity_not_real_path(self):
        """
        Test if filter_sequences_by_ambiguity detects a non-existent path.
        """
        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity("not_real_path")
        self.assertIn("does not exist", str(context.exception))

    def test_filter_sequences_by_ambiguity_empty_dir(self):
        """
        Test if filter_sequences_by_ambiguity detects a empty folder.
        """
        empty_dir = self.test_directory / "empty"
        empty_dir.mkdir()

        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity(empty_dir)
        self.assertIn("No FASTA files found in", str(context.exception))

        shutil.rmtree(empty_dir)

    def test_filter_sequences_by_ambiguity_file(self):
        """
        Test if filter_sequences_by_ambiguity correctly reads a input file.
        """
        test_file = self.test_directory / "ambiguous_sequences.fasta"

        results = filter_sequences_by_ambiguity(
            test_file, max_ambiguous_percentage=0.25
        )

        output_dir = self.test_directory / "filtered"
        self.assertTrue(output_dir.exists())

        output_file = self.test_directory / "filtered" / "ambiguous_sequences.fasta"
        self.assertTrue(output_file.exists())

        with open(output_file) as f:
            contents = f.read()

            self.assertIn(">seq1", contents)
            self.assertIn("ATGCATGC", contents)
            self.assertIn(">seq2", contents)
            self.assertIn("ATGNNTGC", contents)
            self.assertNotIn(">seq3", contents)
            self.assertNotIn("NNNNNNNN", contents)
            self.assertNotIn(">seq4", contents)
            self.assertNotIn(
                "ACGTACGACGAGCATTCGAAGGTCAGTCGNNNNATTAGCTACTGATCGATCGACTAGCTCCGCATCGATGATGCATGCTAGTCGATGCATGCATCG",
                contents,
            )

    def test_extract_primer_binding_sites(self):
        amplicon_file = self.test_directory / "amplicon-test/primerC.fasta"
        insert_file = self.test_directory / "insert-test/primerA.fasta"

        result = extract_primer_binding_sites(amplicon_file, insert_file)

        self.assertIsInstance(result, pd.DataFrame)

        self.assertEqual(
            list(result.columns),
            ["header", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"],
        )

        self.assertEqual(result["header"].iloc[0], ">abc")
        self.assertEqual(result["header"].iloc[1], ">def")
        self.assertEqual(result["fwd_seq"].iloc[0], "ACGTAGCA")
        self.assertEqual(result["rev_seq"].iloc[0], "ACCATCA")
        self.assertEqual(result["fwd_seq"].iloc[1], "AGTGA")
        self.assertEqual(result["rev_seq"].iloc[1], "GTCGATAGCAT")
        self.assertEqual(result["fwd_seq_len"].iloc[0], 8)
        self.assertEqual(result["rev_seq_len"].iloc[0], 7)
        self.assertEqual(result["fwd_seq_len"].iloc[1], 5)
        self.assertEqual(result["rev_seq_len"].iloc[1], 11)

    def tearDown(self):
        """
        Clean up any files created during tests.
        """
        filtered_dir = self.test_directory / "filtered"
        if filtered_dir.exists():
            shutil.rmtree(filtered_dir)


class TestSequenceCountTracking(unittest.TestCase):
    def setUp(self):
        self.analysis_folder = tempfile.TemporaryDirectory()
        self.original_database = tempfile.NamedTemporaryFile(
            delete=False, suffix=".fasta"
        )

        with open(self.original_database.name, "w") as f:
            f.write(">seq\nAGTGCA\n>seq2\nGTCAGCGA\n>seq3\nGGGGCA")

        self.analysis_subfolder = os.path.join(self.analysis_folder.name, "amplicon")
        os.makedirs(self.analysis_subfolder, exist_ok=True)
        self.analysis_file = tempfile.NamedTemporaryFile(
            delete=False,
            dir=self.analysis_subfolder,
            prefix="amplicon_",
            suffix=".fasta",
        )
        with open(self.analysis_file.name, "w") as f:
            f.write(">seq1\nGTCA\n>seq2\nGTCGGG")

    def tearDown(self):
        self.analysis_folder.cleanup()
        if os.path.exists(self.original_database.name):
            os.remove(self.original_database.name)
        if os.path.exists(self.analysis_file.name):
            os.remove(self.analysis_file.name)

    def test_sequence_count_tracking(self):
        df_pivoted = sequence_count_tracking(
            self.original_database.name, self.analysis_folder.name
        )

        self.assertIsInstance(df_pivoted, pd.DataFrame)

        run_name = os.path.basename(self.analysis_folder.name)
        output_path = os.path.join(
            self.analysis_folder.name, f"{run_name}-sequence_count_track.tsv"
        )
        self.assertTrue(os.path.exists(output_path))

        result_df = pd.read_csv(output_path, sep="\t")
        self.assertIn("original_database", result_df.columns)
        self.assertIn("amplicon", result_df.columns)

        self.assertEqual(df_pivoted.loc["original_database", "original_database"], 3)
        self.assertEqual(df_pivoted.iloc[0, 1], 2)


class TestMultiBarcodeToolsInput(unittest.TestCase):
    def setUp(self):
        """
        Set up a temporary directory and create sample FASTA files for testing.
        """
        self.test_dir = tempfile.mkdtemp()
        self.create_sample_fasta_files()

    def tearDown(self):
        """
        Clean up the temporary directory after tests.
        """
        shutil.rmtree(self.test_dir)

    def create_sample_fasta_files(self):
        """
        Method to create sample FASTA files for testing different scenarios.
        """
        # FASTA with standard header
        with open(os.path.join(self.test_dir, "primer1.fasta"), "w") as f:
            f.write(">seq1|Species Name 1|harmonized_species|rank|kingdom|phylum|order|family|genus|species\n")
            f.write("ATCGATCGATCG\n")
            f.write(">seq2|Species Name 2|harmonized_species|rank|kingdom|phylum|order|family|genus|species\n")
            f.write("GCTAGCTAGCTA\n")

        # FASTA with multiple sequences
        with open(os.path.join(self.test_dir, "primer2.fasta"), "w") as f:
            f.write(">seq3|Species Name 3|harmonized_species|rank|kingdom|phylum|order|family|genus|species\n")
            f.write("TAGCTAGCTAGC\n")
            f.write(">seq4|Species Name 4|harmonized_species|rank|kingdom|phylum|order|family|genus|species\n")
            f.write("CGATCGATCGAT\n")

        # FASTA with problematic header
        with open(os.path.join(self.test_dir, "primer3.fasta"), "w") as f:
            f.write(">seq5\n")  # Missing species name
            f.write("ATCGATCG\n")

    def test_create_multibarcode_tools_input(self):
        """
        Method to test the main function create_MultiBarcodeTools_input.
        """
        output_file = os.path.join(self.test_dir, "output.tsv")

        create_MultiBarcodeTools_input(self.test_dir, self.test_dir, output_file)

        self.assertTrue(os.path.exists(output_file))

        with open(output_file, "r") as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 8)

    def test_process_sequence_invalid_header(self):
        """
        Method to test process_sequence with an invalid header (no species name).
        """
        import io
        import sys

        stdout = io.StringIO()
        sys.stdout = stdout

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as mock_tsv:
            try:
                header = "seq1"
                sequence_lines = ["ATCGATCGATCG"]
                primer_name = "test_primer"

                process_sequence(header, sequence_lines, primer_name, mock_tsv)

                mock_tsv.close()

                error_msg = stdout.getvalue().strip()
                self.assertIn("mozaico WARNING: Unexpected header format", error_msg)

                with open(mock_tsv.name, "r") as f:
                    content = f.read().strip()
                    self.assertEqual(content, "")

            finally:
                sys.stdout = sys.__stdout__
                os.unlink(mock_tsv.name)
