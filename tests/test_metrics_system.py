"""
Unit tests for metrics_system.py
"""

import sys
import unittest
import tempfile
import shutil
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, call, mock_open, patch

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

    @patch("builtins.print")
    @patch("os.replace")
    def test_filter_fasta_no_pipe_header(self, mock_replace, mock_print):
        """
        Method to test handling of a header without a pipe character.
        """
        test_fasta_content = """>no_pipe_header
                                ATGCATGCATGC"""

        otl_taxa_ex = set()

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            self.handler.filter_fasta_for_species_not_in_otl(temp_file_path, otl_taxa_ex)
            mock_print.assert_any_call("mozaico WARNING: No '|' found in header - >no_pipe_header")
        finally:
            os.unlink(temp_file_path)

    @patch("builtins.print")
    @patch("os.replace")
    def test_filter_fasta_no_taxa(self, mock_replace, mock_print):
        """
        Method to test handling of a header without a pipe character.
        """
        test_fasta_content = """>notax | \nATGCATGCATGC"""

        otl_taxa_ex = set()

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            self.handler.filter_fasta_for_species_not_in_otl(temp_file_path, otl_taxa_ex)
            mock_print.assert_any_call("mozaico WARNING: Taxonomy seems to not be present for - >notax |")
        finally:
            os.unlink(temp_file_path)

    @patch("os.replace")
    def test_filter_fasta_overwrite_mode(self, mock_replace):
        """
        Method to test the overwrite functionality.
        """
        test_fasta_content = """>valid_header | taxa1\nATGCATGCATGC"""

        otl_taxa_ex = {'taxa1'}

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            total_count, kept_count, output_file = self.handler.filter_fasta_for_species_not_in_otl(
                temp_file_path, otl_taxa_ex, overwrite=True
            )

            mock_replace.assert_called_once_with(
                temp_file_path + ".temp",
                temp_file_path
            )

            self.assertEqual(total_count, 1)
            self.assertEqual(kept_count, 1)
            self.assertEqual(output_file, temp_file_path)
        finally:
            os.unlink(temp_file_path)

    def test_filter_fasta_non_overwrite_mode(self):
        """
        Test function when overwrite is set to False.
        """
        test_fasta_content = """>valid_header | taxa1\nATGCATGCATGC"""
        otl_taxa_ex = {'taxa1'}

        test_dir = os.path.join(os.path.dirname(__file__), 'test_data')
        os.makedirs(test_dir, exist_ok=True)

        input_file_path = os.path.join(test_dir, "test_input.fasta")
        with open(input_file_path, 'w') as f:
            f.write(test_fasta_content)

        try:
            total_count, kept_count, output_file = self.handler.filter_fasta_for_species_not_in_otl(
                input_file_path, otl_taxa_ex, overwrite=False
            )
            self.assertEqual(total_count, 1)
            self.assertEqual(kept_count, 1)

            assert os.path.exists(output_file)
            assert output_file.endswith("test_input.fasta")
            assert "_otl_filtered" in output_file

            with open(output_file, 'r') as f:
                content = f.read().strip()
            assert content == test_fasta_content
            assert os.path.dirname(output_file) != os.path.dirname(input_file_path), "Output file should be in a different directory"

        finally:
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            if 'output_file' in locals() and os.path.exists(output_file):
                os.remove(output_file)
            if os.path.exists(test_dir) and not os.listdir(test_dir):
                os.rmdir(test_dir)

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
            self.total_otl_taxa_count,
            barcode_threshold=2,
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

    def test_parse_files_with_no_matching_base_names(self):
        with patch("os.listdir") as mock_listdir, patch("builtins.print") as mock_print:
            mock_listdir.side_effect = [
                ["file1.fasta", "file2.txt"],
                ["file3.fasta", "file4.txt"],
            ]

            result = self.binding.parse_files_with_same_extension_in_folders(
                "/path/A", "/path/B"
            )

            self.assertIsNone(result)

            mock_print.assert_called_once_with(
                "mozaiko ERROR: No matching files found between the two folders. Check folder directory and content before re-running."
            )

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

    def test_primer_pbs_analysis_no_matching_files(self):
        with patch.object(
            self.binding, "parse_files_with_same_extension_in_folders", return_value=[]
        ) as mock_parse_files:

            result = self.binding.primer_pbs_analysis(
                "amplicon_folder",
                "insert_folder",
                "data/test_data/test_primer_table.tsv",
            )
            self.assertEqual(result, (None, None))

    def test_add_missing_otl_taxa_to_df_with_values_of_zero(self):
        input_df = pd.DataFrame({
            "taxon": ["taxon1", "taxon2"],
            "seq_id": ["seq1", "seq2"],
            "analysis1": [10, 20],
            "analysis2": [5, 15]
        })

        otl_taxa_set = {"taxon1", "taxon2", "taxon3", "taxon4"}

        expected_df = pd.DataFrame({
            "taxon": ["taxon1", "taxon2", "taxon3", "taxon4"],
            "seq_id": ["seq1", "seq2", "otl-import", "otl-import"],
            "analysis1": [10, 20, 0, 0],
            "analysis2": [5, 15, 0, 0]
        })

        otl_populated_df = self.binding.add_missing_otl_taxa_to_df_with_values_of_zero(input_df, otl_taxa_set)

        pd.testing.assert_frame_equal(otl_populated_df.sort_values(by="taxon").reset_index(drop=True),
                                       expected_df.sort_values(by="taxon").reset_index(drop=True))

    def test_iterate_over_primer_pbs_df_with_otl_taxa(self):
        comprehensive_primer_dfs = {
            "primer_1": pd.DataFrame({
                "taxon": ["taxon1", "taxon2"],
                "seq_id": ["seq1", "seq2"],
                "analysis1": [10, 20],
                "analysis2": [5, 15],
            }),
            "primer_2": pd.DataFrame({
                "taxon": ["taxon3", "taxon4"],
                "seq_id": ["seq3", "seq4"],
                "analysis1": [30, 40],
                "analysis2": [25, 35],
            }),
        }

        otl_taxa_set = {"taxon1", "taxon2", "taxon3", "taxon4", "taxon5"}

        expected_primer_1_df = pd.DataFrame({
            "taxon": ["taxon1", "taxon2", "taxon5"],
            "seq_id": ["seq1", "seq2", "otl-import"],
            "analysis1": [10, 20, 0],
            "analysis2": [5, 15, 0],
        })
        expected_primer_2_df = pd.DataFrame({
            "taxon": ["taxon3", "taxon4", "taxon5"],
            "seq_id": ["seq3", "seq4", "otl-import"],
            "analysis1": [30, 40, 0],
            "analysis2": [25, 35, 0],
        })

        def side_effect_check(df, otl_taxa):
            if df.equals(comprehensive_primer_dfs["primer_1"]) and otl_taxa == otl_taxa_set:
                return expected_primer_1_df
            elif df.equals(comprehensive_primer_dfs["primer_2"]) and otl_taxa == otl_taxa_set:
                return expected_primer_2_df
            raise ValueError("Unexpected call to mocked method.")

        with patch.object(
            self.binding,
            "add_missing_otl_taxa_to_df_with_values_of_zero",
            side_effect=side_effect_check
        ) as mock_add_otl_taxa:
            processed_primers = self.binding.iterate_over_primer_pbs_df(
                comprehensive_primer_dfs,
                add_otl_taxa=True,
                otl_taxa_set=otl_taxa_set
            )

            self.assertIn("primer_1", processed_primers.keys())
            self.assertIn("primer_2", processed_primers.keys())

            pd.testing.assert_frame_equal(
                processed_primers["primer_1"].sort_values(by="taxon").reset_index(drop=True),
                expected_primer_1_df.sort_values(by="taxon").reset_index(drop=True),
            )
            pd.testing.assert_frame_equal(
                processed_primers["primer_2"].sort_values(by="taxon").reset_index(drop=True),
                expected_primer_2_df.sort_values(by="taxon").reset_index(drop=True),
            )

    def test_iterate_over_primer_pbs_df_without_otl_taxa(self):
        comprehensive_primer_dfs = {
            "primer_1": pd.DataFrame({
                "taxon": ["taxon1", "taxon2"],
                "seq_id": ["seq1", "seq2"],
                "analysis1": [10, 20],
                "analysis2": [5, 15],
            }),
            "primer_2": pd.DataFrame({
                "taxon": ["taxon3", "taxon4"],
                "seq_id": ["seq3", "seq4"],
                "analysis1": [30, 40],
                "analysis2": [25, 35],
            }),
        }

        processed_primers = self.binding.iterate_over_primer_pbs_df(
            comprehensive_primer_dfs,
            add_otl_taxa=False,
        )

        self.assertIn("primer_1", processed_primers)
        self.assertIn("primer_2", processed_primers)

        pd.testing.assert_frame_equal(
            processed_primers["primer_1"].sort_values(by="taxon").reset_index(drop=True),
            comprehensive_primer_dfs["primer_1"].sort_values(by="taxon").reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            processed_primers["primer_2"].sort_values(by="taxon").reset_index(drop=True),
            comprehensive_primer_dfs["primer_2"].sort_values(by="taxon").reset_index(drop=True),
        )


class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        self.all_inserts_folder = "data/test_data/test-folder-metrics"
        self.otl = "data/test_data/test_otl.tsv"
        self.metric_sys_ex = MetricsSystemExecutor(self.all_inserts_folder, self.otl)

    def test_calculate_reference_database_quality(self):
        ref_bd_scores = self.metric_sys_ex.calculate_reference_database_quality()

        expected_scores = {
            "test_amplicon_reffb": {
                "barcoded_taxa_five_plus": 16.67,
                "ratio_barcoded_taxa": 0.33,
            }
        }
        self.assertEqual(ref_bd_scores, expected_scores)
