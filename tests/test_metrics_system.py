"""
Unit tests for metrics_system.py
"""
import unittest

from pathlib import Path
from src.marker_scoring.metrics_system import *

class TestReferenceDatabaseQuality(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path("data/test_data")
        self.inserts_file = str(self.test_directory / "test_amplicon_reffb.fasta")
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.ref_bd_cls = ReferenceDatabaseQuality(self.inserts_file, self.otl)

    def test_calculate_number_of_barcodes_per_taxon(self):
        barcodes = self.ref_bd_cls.calculate_number_of_barcodes_per_taxon()

        expected_barcodes = {
            "speciesA": 3,
            "speciesB": 1,
            "speciesC": 1,
            "speciesD": 2
        }

        self.assertEqual(barcodes, expected_barcodes)

    def test_calculate_percentage_of_taxa_w_1_barcode(self):
        percentage = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(total_taxa=6)

        expected_percentage = 33.33

        self.assertEqual(percentage, expected_percentage)

    def test_calculate_percentage_of_taxa_w_2_barcode(self):
        self.ref_bd_cls.import_otl()
        percentage2 = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(barcode_threshold=2)

        expected_percentage2 = 16.67

        self.assertEqual(percentage2, expected_percentage2)

class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        pass