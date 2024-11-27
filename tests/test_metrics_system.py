"""
Unit tests for metrics_system.py
"""

import sys
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, call, patch, mock_open

from pandas._testing import assert_frame_equal

from src.marker_scoring.metrics_system import *
from src.marker_scoring.scoring_utils import *


class TestOtlHandler(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path("data/test_data")
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.handler = OtlHandler(self.otl)

    def test_validate_otl_not_exist(self):
        with self.assertRaises(SystemExit) as cm, patch(
            "sys.stdout", new=StringIO()
        ) as mock_stdout:
            no_otl = str(self.test_directory / "not-there.tsv")
            self.handler.validate_otl(no_otl)

            self.assertEqual(cm.exception.code, 1)

        expected_message = "mozaiko INFO: The OTL does not exist. Exiting...\n"

        captured_output = mock_stdout.getvalue()
        self.assertEqual(expected_message, captured_output)

    def test_validate_otl_not_tsv(self):
        with self.assertRaises(SystemExit) as cm, patch(
            "sys.stdout", new=StringIO()
        ) as mock_stdout:
            no_otl = str(self.test_directory / "empty_fasta.fasta")
            self.handler.validate_otl(no_otl)

            self.assertEqual(cm.exception.code, 1)

        expected_message = "mozaiko ERROR: The OTL must be a TSV file. Exiting...\n"

        captured_output = mock_stdout.getvalue()

    def test_validate_otl_empty_rows(self):
        self.handler.validate_otl(self.otl)

        otl_table = pd.read_csv(self.otl, sep="\t", header=0)
        self.assertEqual(len(otl_table), len(otl_table.dropna(subset=["taxa"])))

    @patch("builtins.input", return_value="data/test_data/test_otl.tsv")
    @patch("sys.stdout", new_callable=StringIO)
    def test_import_otl(self, mock_stdout, mock_input):
        otl_handler = OtlHandler()
        total_taxa_count, unique_taxa = otl_handler.import_otl()

        self.assertIn(
            "mozaiko INFO: To continue the evaluation, a Operational Taxonomic List (OTL) is \
                    required. An OTL is a list contaning information on the taxonomic numenclature \
                        of all identifiable taxa in routine biomonitoring initiatives.",
            mock_stdout.getvalue(),
        )

        otl = pd.read_csv("data/test_data/test_otl.tsv", sep="\t", header=0)

        pd.testing.assert_frame_equal(otl_handler.otl, otl)

        expected_unique_taxa = set(otl["taxa"])
        self.assertEqual(total_taxa_count, len(expected_unique_taxa))
        self.assertEqual(unique_taxa, expected_unique_taxa)


class TestReferenceDatabaseQuality(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path("data/test_data")
        self.inserts_file = str(self.test_directory / "test-folder-metrics")
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.ref_bd_cls = ReferenceDatabaseQuality(self.inserts_file, self.otl)
        self.handler = OtlHandler(self.otl)
        self.handler.import_otl()
        self.total_otl_taxa_count = self.handler.total_taxa

    @patch("builtins.input", return_value="data/test_data/test-folder-metrics")
    @patch("sys.stdout", new_callable=StringIO)
    def test_calculate_number_of_barcodes_per_taxon_input(
        self, mock_stdout, mock_input
    ):
        with patch("src.marker_scoring.metrics_system.OtlHandler") as mock_otl_handler:
            mock_instance = MagicMock()
            mock_otl_handler.return_value = mock_instance

            no_file_class = ReferenceDatabaseQuality()

            with patch("sys.exit") as mock_exit:
                mock_exit.side_effect = SystemExit(1)
                try:
                    no_file_class.calculate_number_of_barcodes_per_taxon()
                except SystemExit:
                    pass

            self.assertIn(
                "mozaiko INFO: Please provide a folder containing FASTA files.",
                mock_stdout.getvalue(),
            )

            self.assertEqual(
                no_file_class.all_inserts_path, "data/test_data/test-folder-metrics"
            )

    def test_calculate_number_of_barcodes_per_taxon(self):
        barcodes = self.ref_bd_cls.calculate_number_of_barcodes_per_taxon()
        expected_barcodes = {
            "test_amplicon_reffb": {
                "speciesA": 8,
                "speciesB": 2,
                "speciesC": 1,
                "speciesD": 2,
            }
        }
        self.assertEqual(barcodes, expected_barcodes)

    def test_calculate_percentage_of_taxa_w_1_barcode(self):
        percentage = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(
            self.total_otl_taxa_count, barcode_threshold=1
        )

        expected_percentage = {"test_amplicon_reffb": 50.0}

        self.assertEqual(percentage, expected_percentage)

    def test_calculate_percentage_of_taxa_w_2_barcode(self):
        percentage2 = self.ref_bd_cls.calculate_percentage_of_taxa_w_x_barcodes(
            self.total_otl_taxa_count, barcode_threshold=2,
        )

        expected_percentage2 = {"test_amplicon_reffb": 16.67}

        self.assertEqual(percentage2, expected_percentage2)

    def test_ratio_barcoded_taxa(self):
        rbt_rounded = self.ref_bd_cls.barcoded_taxa_ratio(self.total_otl_taxa_count)
        print(rbt_rounded)

        expected_output = {
            "test_amplicon_reffb": {
                "barcoded_taxa_five_plus": 16.67,
                "ratio_barcoded_taxa": 0.33,
            }
        }

        self.assertEqual(rbt_rounded, expected_output)


class TestBinding(unittest.TestCase):
    def setUp(self):
        self.binding = Binding()

    def test_init_with_default_mismatches(self):
        """
        Test Binding class initialization with the default number of mismatches.
        """
        with patch.object(
            InSilicoAmplification, "get_number_of_mismatches", return_value=2
        ) as mock_get_mismatches:
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
        with patch("os.listdir") as mock_listdir, patch(
            "os.path.join"
        ) as mock_path_join, patch("os.path.exists", return_value=True):

            mock_listdir.side_effect = [
                ["file1.fasta", "file2.txt"],
                ["file1.txt", "file2.fasta"],
            ]
            mock_path_join.side_effect = [
                "/path/A/file1.fasta",
                "/path/B/file1.txt",
                "/path/A/file2.txt",
                "/path/B/file2.fasta",
            ]

            result = self.binding.parse_files_with_same_extension_in_folders(
                "/path/A", "/path/B"
            )

            self.assertEqual(len(result), 2)
            self.assertEqual(result[0], ("/path/A/file1.fasta", "/path/B/file1.txt"))

    def test_get_primer_table(self):
        with patch("pandas.read_csv") as mock_read_csv, patch.object(
            self.binding.amplification_instance, "validate_primer_table"
        ) as mock_validate:

            mock_primer_table = pd.DataFrame(
                {
                    "barcode_region": ["region1"],
                    "assay_name": ["assay1"],
                    "fwd_seq": ["ACGT"],
                    "rev_seq": ["GCTA"],
                }
            )
            mock_read_csv.return_value = mock_primer_table

            result = self.binding.get_primer_table("primer_table.tsv")

            mock_validate.assert_called_once_with("primer_table.tsv")
            pd.testing.assert_frame_equal(result, mock_primer_table)
            self.assertTrue(hasattr(self.binding, "primer_table"))

    def test_pbs_table(self):
        with patch(
            "src.marker_scoring.scoring_utils.extract_primer_binding_sites"
        ) as mock_extract, patch("builtins.open", create=True):
            mock_pbs_table = pd.DataFrame(
                {"header": [">seq1|taxon"], "fwd_seq": ["ACGT"], "rev_seq": ["GCTA"]}
            )

            mock_extract.return_value = mock_pbs_table

            result = self.binding.get_pbs_table("amplicon.fasta", "insert.fasta")

            self.assertTrue("primer_name" in result.columns)

    @patch('src.marker_scoring.metrics_system.gc_fraction')
    @patch('src.marker_scoring.metrics_system.MeltingTemp.Tm_GC')
    @patch('src.marker_scoring.metrics_system.calculate_iupac_mismatches')
    @patch('src.marker_scoring.metrics_system.os.path.splitext')
    @patch('src.marker_scoring.metrics_system.os.path.basename')
    def test_primer_pbs_analysis(self,
        mock_basename,
        mock_splitext,
        mock_calculate_mismatches,
        mock_melting_temp,
        mock_gc_fraction
    ):
        mock_basename.return_value = 'rbcL_assay1'
        mock_splitext.return_value = ('rbcL_assay1', '.fasta')
        mock_melting_temp.return_value = 60.0
        mock_gc_fraction.return_value = 0.5

        def mock_calculate_mismatches_side_effect(*args, **kwargs):
            if kwargs.get('search_gc_clamp', False):
                return (1, 1)  # Return a tuple when search_gc_clamp is True
            return 1  # Return an integer if False

        mock_calculate_mismatches.side_effect = mock_calculate_mismatches_side_effect

        primer_table = pd.DataFrame({
            'barcode_region': ['rbcL'],
            'assay_name': ['assay1'],
            'fwd_seq': ['ACGTACGT'],
            'rev_seq': ['TGCATGCA']
        })

        pbs_table = pd.DataFrame({
            'header': ['>seq1|Taxon1'],
            'fwd_seq': ['ACGTACfGT'],
            'rev_seq': ['TGCATGCA']
        })

        with patch.object(self.binding, 'get_primer_table') as mock_get_primer_table, \
             patch.object(self.binding, 'parse_files_with_same_extension_in_folders',
                          return_value=[('rbcL_assay1.fasta', 'insert_rbcL_assay1.fasta')]) as mock_parse_files, \
             patch.object(self.binding, 'get_pbs_table', return_value=pbs_table) as mock_get_pbs_table, \
             patch('builtins.open', mock_open()) as mock_file:

            self.binding.primer_table = primer_table

            result = self.binding.primer_pbs_analysis(
                'amplicon_folder',
                'insert_folder',
                'primer_table.tsv',
                save_results=True
            )
        self.assertIsNotNone(result)
        self.assertIn('rbcL_assay1', result)

        primer_props = result['rbcL_assay1']['primer_properties']
        self.assertEqual(primer_props['forward_primer']['gc_fraction'], 0.5)
        self.assertEqual(primer_props['forward_primer']['melting_temp'], 60.0)

        self.assertTrue('full_mismatches' in result['rbcL_assay1'])
        self.assertTrue('three_end_mismatches' in result['rbcL_assay1'])
        self.assertTrue('three_end_gc_matches' in result['rbcL_assay1'])

    def test_primer_pbs_analysis_no_matching_files(self):
        # Patch to return no matching files between folders
        with patch.object(self.binding, 'parse_files_with_same_extension_in_folders',
                          return_value=[]) as mock_parse_files:

            result = self.binding.primer_pbs_analysis(
                'amplicon_folder',
                'insert_folder',
                'data/test_data/test_primer_table.tsv'
            )
            self.assertIsNone(result)

    @patch('json.dump')
    def test_primer_pbs_analysis_save_results(self, mock_json_dump):
        primer_table = pd.DataFrame({
            'barcode_region': ['rbcL'],
            'assay_name': ['assay1'],
            'fwd_seq': ['ACGTACGT'],
            'rev_seq': ['TGCATGCA']
        })

        pbs_table = pd.DataFrame({
            'header': ['>seq1|Taxon1'],
            'fwd_seq': ['ACGTACGT'],
            'rev_seq': ['TGCATGCA']
        })

        with patch.object(self.binding, 'get_primer_table') as mock_get_primer_table, \
             patch.object(self.binding, 'parse_files_with_same_extension_in_folders',
                          return_value=[('rbcL_assay1.fasta', 'insert_rbcL_assay1.fasta')]) as mock_parse_files, \
             patch.object(self.binding, 'get_pbs_table', return_value=pbs_table) as mock_get_pbs_table, \
             patch('builtins.open', mock_open()) as mock_file:

            self.binding.primer_table = primer_table

            result = self.binding.primer_pbs_analysis(
                'amplicon_folder',
                'insert_folder',
                'data/test_data/test_primer_table.tsv',
                save_results=True
            )

            mock_json_dump.assert_called_once()

    def test_get_max_mismatches_per_taxon(self):
        mock_mismatches = {
            "primer_pair1": [
                {"taxon": "taxon1", "mismatch_sum": 2},
                {"taxon": "taxon1", "mismatch_sum": 3},
                {"taxon": "taxon2", "mismatch_sum": 1},
            ]
        }

        result = self.binding.get_max_mismatches_per_taxon(mock_mismatches)
        self.assertEqual(result["primer_pair1"]["taxon1"], 3)
        self.assertEqual(result["primer_pair1"]["taxon2"], 1)

    def test_get_max_mismatche_across_taxon(self):
        mock_max_mismatches = {
            "primer_pair1": {"taxon1": 3, "taxon2": 2},
            "primer_pair2": {"taxon3": 1, "taxon4": 4},
        }

        result = self.binding.get_max_mismatches_count_across_taxon(mock_max_mismatches)

        self.assertEqual(result["primer_pair1"], 5)
        self.assertEqual(result["primer_pair2"], 5)

    def test_get_priming_ratio(self):
        mock_mismatches_all_len = {"primer_pair1": {"taxon1": 10, "taxon2": 5}}
        mock_mismatches_three_end = {"primer_pair1": {"taxon1": 2, "taxon2": 1}}

        result = self.binding.get_priming_ratio(
            mock_mismatches_all_len, mock_mismatches_three_end
        )

        self.assertEqual(result["primer_pair1"], 0.4)


class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        self.all_inserts_folder = "data/test_data/test-folder-metrics"
        self.otl = "data/test_data/test_otl.tsv"
        self.metric_sys_ex = MetricsSystemExecutor(self.all_inserts_folder, self.otl)

    def test_calculate_reference_database_quality(self):
        ref_bd_scores = self.metric_sys_ex.calculate_reference_database_quality()

        expected_scores = {
            'test_amplicon_reffb':
                {'barcoded_taxa_five_plus': 16.67,
                'ratio_barcoded_taxa': 0.33}
             }
        self.assertEqual(ref_bd_scores, expected_scores)
