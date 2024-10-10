"""
Unit tests for scoring_utils.py
"""

import shutil
import unittest
from pathlib import Path
from unittest.mock import mock_open

from src.marker_scoring.scoring_utils import *


class TestScoringUtils(unittest.TestCase):
    def setUp(self):
        self.test_directory = "data/test_data"

    def test_calculate_iupac_mismatches(self):
        """
        Function to test the calculate_iupac_mismatches function.
        """
        mismatches_cases = self.test_directory + "/iupac_mismatches_cases.txt"

        with open(mismatches_cases, "r") as file:
            for line in file:
                seq1, seq2, expected_mismatches = line.strip().split(" ")
                mismatches = calculate_iupac_mismatches(seq1, seq2)
                self.assertEqual(mismatches, int(expected_mismatches))

    def test_filter_sequences_by_ambiguity(self):
        """
        Function to test the filter_sequences_by_ambiguity mehtod.
        Check if sequences that suprass the threshold for ambiguity percentage are filtered from
        the input fasta.
        """
        sequences = ["ACGTACGTAG", "ACGTANGATC", "RWNBVDNWSN"]  # 0%  # 10%  # 100%

    def test_filter_sequences_by_ambiguity_not_real_path(self):
        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity("not_real_path")
        self.assertIn("does not exist", str(context.exception))

    def test_filter_sequences_by_ambiguity_empty_path(self):
        empty_dir = Path(self.test_directory) / "empty"
        empty_dir.mkdir()

        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity(empty_dir)
        self.assertIn("No FASTA files found in", str(context.exception))

        shutil.rmtree(empty_dir)
