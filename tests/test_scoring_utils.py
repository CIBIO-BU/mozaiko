"""
Unit tests for scoring_utils.py
"""

import shutil
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
                seq1, seq2, expected_mismatches = line.strip().split(" ")
                mismatches = calculate_iupac_mismatches(seq1, seq2)
                self.assertEqual(mismatches, int(expected_mismatches))

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

    def test_extract_primer_binding_sites(self):
        amplicon_file = self.test_directory / "amplicon_test.fasta"
        insert_file = self.test_directory / "insert_test.fasta"

        result = extract_primer_binding_sites(amplicon_file, insert_file)

        self.assertIsInstance(result, pd.DataFrame)

        self.assertEqual(
            list(result.columns),
            ["record", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"],
        )

        self.assertEqual(result["record"].iloc[0], "> abc")
        self.assertEqual(result["record"].iloc[1], "> def")
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
