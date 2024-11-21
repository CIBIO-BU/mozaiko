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
from src.marker_scoring.scoring_utils import *


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

        expected_output = (16.67, 3.0)

        self.assertEqual(rbt_rounded, expected_output)

class TestBinding(unittest.TestCase):
    def setUp(self):
        self.binding = Binding()

    def test_init_with_default_mismatches(self):
        """
        Test Binding class initialization with the default number of mismatches.
        """
        with patch.object(InSilicoAmplification, 'get_number_of_mismatches', return_value = 2) as mock_get_mismatches:
            binding = Binding()
            mock_get_mismatches.assert_called_once()
            self.assertEqual(binding.number_of_mismatches, 2)

    def test_init_with_custom_mismatches(self):
        """
        Test Binding class initialization with custom number of mismatches.
        """
        binding = Binding(number_of_mismatches=3)
        self.assertEqual(binding.number_of_mismatches, 3)

    def test_parse_files_with_same_extension_in_folder(self):
        with patch('os.listdir') as mock_listdir, \
            patch('os.path.join') as mock_path_join, \
            patch('os.path.exists', return_value=True):

            mock_listdir.side_effect = [
                ['file1.fasta', 'file2.txt'],
                ['file1.txt', 'file2.fasta']
            ]
            mock_path_join.side_effect = [
                '/path/A/file1.fasta',
                '/path/B/file1.txt',
                '/path/A/file2.txt',
                '/path/B/file2.fasta'
            ]

            result = self.binding.parse_files_with_same_extension_in_folders('/path/A', '/path/B')

            self.assertEqual(len(result), 2)
            self.assertEqual(result[0], ('/path/A/file1.fasta', '/path/B/file1.txt'))

    def test_get_primer_table(self):
        with patch('pandas.read_csv') as mock_read_csv, \
        patch.object(self.binding.amplification_instance, 'validate_primer_table') as mock_validate:

            mock_primer_table = pd.DataFrame({
                'barcode_region': ['region1'],
                'assay_name': ['assay1'],
                'fwd_seq': ['ACGT'],
                'rev_seq': ['GCTA']
            })
            mock_read_csv.return_value = mock_primer_table

            result = self.binding.get_primer_table('primer_table.tsv')

            mock_validate.assert_called_once_with('primer_table.tsv')
            pd.testing.assert_frame_equal(result, mock_primer_table)
            self.assertTrue(hasattr(self.binding, 'primer_table'))

    def test_pbs_table(self):
        with patch('src.marker_scoring.scoring_utils.extract_primer_binding_sites') as mock_extract, \
            patch('builtins.open', create=True):
            mock_pbs_table = pd.DataFrame({
                'header': ['>seq1|taxon'],
                'fwd_seq': ['ACGT'],
                'rev_seq': ['GCTA']
            })

            mock_extract.return_value = mock_pbs_table

            result = self.binding.get_pbs_table('amplicon.fasta', 'insert.fasta')

            self.assertTrue('primer_name' in result.columns)

    def test_get_number_of_primer_pbs_mismatches(self):
        with patch.object(self.binding, 'get_primer_table'), \
        patch.object(self.binding, 'parse_files_with_same_extension_in_folders', return_value=[('amplicon.fasta', 'insert.fasta')]), \
        patch.object(self.binding, 'get_pbs_table', return_value=pd.DataFrame({
             'header': ['>seq1|taxon1'],
             'fwd_seq': ['ACGT'],
             'rev_seq': ['GCTA']
         })), \
        patch('src.marker_scoring.scoring_utils.calculate_iupac_mismatches', return_value = 1):

            self.binding.primer_table = pd.DataFrame({
                'barcode_region': ['region1'],
                'assay_name': ['assay1'],
                'fwd_seq': ['ACGT'],
                'rev_seq': ['GCTA']
            })

            calculate_ambiguous_percentage.return_value = 1

            result = self.binding.get_number_of_primer_pbs_mismatches(
                'amplicon_folder', 'insert_folder', 'primer_table.tsv'
            )

            self.assertIsNotNone(result)

    def test_get_max_mismatches_per_taxon(self):
        mock_mismatches = {
            'primer_pair1': [
                {'taxon': 'taxon1', 'mismatch_sum': 2},
                {'taxon': 'taxon1', 'mismatch_sum': 3},
                {'taxon': 'taxon2', 'mismatch_sum': 1}
            ]
        }

        result = self.binding.get_max_mismatches_per_taxon(mock_mismatches)
        self.assertEqual(result['primer_pair1']['taxon1'], 3)
        self.assertEqual(result['primer_pair1']['taxon2'], 1)

    def test_get_max_mismatche_across_taxon(self):
        mock_max_mismatches = {
            'primer_pair1': {'taxon1': 3, 'taxon2': 2},
            'primer_pair2': {'taxon3': 1, 'taxon4': 4}
        }

        result = self.binding.get_max_mismatches_count_across_taxon(mock_max_mismatches)

        self.assertEqual(result['primer_pair1'], 5)
        self.assertEqual(result['primer_pair2'], 5)

    def test_get_priming_ratio(self):
        mock_mismatches_all_len = {
            'primer_pair1': {'taxon1': 10, 'taxon2': 5}
        }
        mock_mismatches_three_end = {
            'primer_pair1': {'taxon1': 2, 'taxon2': 1}
        }

        result = self.binding.get_priming_ratio(mock_mismatches_all_len, mock_mismatches_three_end)

        self.assertEqual(result['primer_pair1'], 0.4)

class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        self.all_inserts_folder = "data/test_data"
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