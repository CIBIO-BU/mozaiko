"""
Unit tests for metrics_system.py
"""

import sys
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import call, patch

from pandas._testing import assert_frame_equal

from src.marker_scoring.metrics_system import *


class TestReferenceDatabaseQuality(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path("data/test_data")
        self.inserts_file = str(self.test_directory / "test_amplicon_reffb.fasta")
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.ref_bd_cls = ReferenceDatabaseQuality(self.inserts_file, self.otl)

    def test_validate_otl_not_exist(self):
        with self.assertRaises(SystemExit) as cm, patch(
            "sys.stdout", new=StringIO()
        ) as mock_stdout:
            no_otl = str(self.test_directory / "not-there.tsv")
            self.ref_bd_cls._validate_otl(no_otl)

            self.assertEqual(cm.exception.code, 1)

        expected_message = "mozaiko INFO: The OTL does not exist. Exiting...\n"

        captured_output = mock_stdout.getvalue()
        self.assertEqual(expected_message, captured_output)

    def test_validate_otl_not_tsv(self):
        with self.assertRaises(SystemExit) as cm, patch(
            "sys.stdout", new=StringIO()
        ) as mock_stdout:
            no_otl = str(self.test_directory / "empty_fasta.fasta")
            self.ref_bd_cls._validate_otl(no_otl)

            self.assertEqual(cm.exception.code, 1)

        expected_message = "mozaiko ERROR: The OTL must be a TSV file. Exiting...\n"

        captured_output = mock_stdout.getvalue()
        # self.assertEqual(expected_message, captured_output)

    def test_validate_otl_empty_rows(self):
        self.ref_bd_cls._validate_otl(self.otl)

        otl_table = pd.read_csv(self.otl, sep="\t", header=0)
        self.assertEqual(len(otl_table), len(otl_table.dropna(subset=["taxa"])))

    @patch("builtins.input", return_value="data/test_data/test_amplicon_reffb.fasta")
    @patch("sys.stdout", new_callable=StringIO)
    def test_calculate_number_of_barcodes_per_taxon_input(
        self, mock_stdout, mock_input
    ):
        no_file_class = ReferenceDatabaseQuality()
        no_file_class.calculate_number_of_barcodes_per_taxon()

        self.assertIn(
            "mozaiko INFO: To continue the evaluation, please provide a FASTA file containing \
                    all insert regions present in the original database, whether successfully \
                        amplified or not.",
            mock_stdout.getvalue(),
        )
        self.assertEqual(
            no_file_class.all_inserts_file, "data/test_data/test_amplicon_reffb.fasta"
        )

    def test_calculate_number_of_barcodes_per_taxon(self):
        barcodes = self.ref_bd_cls.calculate_number_of_barcodes_per_taxon()

        expected_barcodes = {"speciesA": 8, "speciesB": 2, "speciesC": 1, "speciesD": 2}

        self.assertEqual(barcodes, expected_barcodes)

    @patch("builtins.input", return_value="data/test_data/test_otl.tsv")
    @patch("sys.stdout", new_callable=StringIO)
    def test_import_otl(self, mock_stdout, mock_input):
        no_file_class = ReferenceDatabaseQuality()
        no_file_class.import_otl()

        self.assertIn(
            "mozaiko INFO: To continue the evaluation, a Operational Taxonomic List (OTL) is \
                    required. An OTL is a list contaning information on the taxonomic numenclature \
                        of all identifiable taxa in routine biomonitoring initiatives.",
            mock_stdout.getvalue(),
        )

        otl = pd.read_csv("data/test_data/test_otl.tsv", sep="\t", header=0)

        pd.testing.assert_frame_equal(no_file_class.otl, otl)

    def test_calculate_percentage_of_taxa_w_1_barcode(self):
        percentage = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(
            total_taxa=6
        )

        expected_percentage = 50.0

        self.assertEqual(percentage, expected_percentage)

    def test_calculate_percentage_of_taxa_w_2_barcode(self):
        self.ref_bd_cls.import_otl()
        percentage2 = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(
            barcode_threshold=2
        )

        expected_percentage2 = 16.67

        self.assertEqual(percentage2, expected_percentage2)

    def test_ratio_barcoded_taxa(self):
        rbt_rounded = self.ref_bd_cls.ratio_barcoded_taxa()

        expected_bcs_rounded = 3.0
        expected_output = (16.67, 3.0)

        self.assertEqual(rbt_rounded, expected_output)


class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        self.all_inserts_folder = "data/test_data"  # mock
        self.otl = "data/test_data/test_otl.tsv"
        self.metric_sys_ex = MetricsSystemExecutor(self.all_inserts_folder, self.otl)

    @patch("os.listdir")
    @patch("os.path.join")
    @patch("os.path.exists")
    @patch("src.marker_scoring.metrics_system.ReferenceDatabaseQuality")
    def test_calculate_reference_database_quality(
        self, mock_ref_db_class, mock_path_exists, mock_path_join, mock_os_listdir
    ):
        mock_os_listdir.return_value = ["primer1.fasta", "primer2.fasta"]

        mock_ref_db_object = mock_ref_db_class.return_value
        mock_ref_db_object.ratio_barcoded_taxa.side_effect = [(12, 0.5), (14, 0.8)]

        mock_path_exists.return_value = True

        ref_bd_scores = self.metric_sys_ex.calculate_reference_database_quality()

        mock_os_listdir.assert_called_once_with(self.all_inserts_folder)

        mock_path_join.assert_has_calls(
            [
                call(self.all_inserts_folder, "primer1.fasta"),
                call(self.all_inserts_folder, "primer2.fasta"),
            ]
        )

        mock_ref_db_object.ratio_barcoded_taxa.assert_has_calls([call()] * 2)
        self.assertEqual(ref_bd_scores, {"primer1": (12, 0.5), "primer2": (14, 0.8)})
