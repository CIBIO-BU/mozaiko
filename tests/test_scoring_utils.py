"""
Unit tests for scoring_utils.py
"""
import unittest

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
                seq1, seq2, expected_mismatches = line.strip().split("\t")
                mismatches = calculate_iupac_mismatches(seq1, seq2)
                self.assertEqual(mismatches, int(expected_mismatches))