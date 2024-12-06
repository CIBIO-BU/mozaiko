"""
Unit tests for metrics_system.py
"""

import shutil
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, call, mock_open, patch

from Bio.SeqRecord import SeqRecord
from pandas._testing import assert_frame_equal

import src
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

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            self.handler.filter_fasta_for_species_not_in_otl(
                temp_file_path, otl_taxa_ex
            )
            mock_print.assert_any_call(
                "mozaico WARNING: No '|' found in header - >no_pipe_header"
            )
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

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            self.handler.filter_fasta_for_species_not_in_otl(
                temp_file_path, otl_taxa_ex
            )
            mock_print.assert_any_call(
                "mozaico WARNING: Taxonomy seems to not be present for - >notax |"
            )
        finally:
            os.unlink(temp_file_path)

    @patch("os.replace")
    def test_filter_fasta_overwrite_mode(self, mock_replace):
        """
        Method to test the overwrite functionality.
        """
        test_fasta_content = """>valid_header | taxa1\nATGCATGCATGC"""

        otl_taxa_ex = {"taxa1"}

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_file:
            temp_file.write(test_fasta_content)
            temp_file_path = temp_file.name

        try:
            total_count, kept_count, output_file = (
                self.handler.filter_fasta_for_species_not_in_otl(
                    temp_file_path, otl_taxa_ex, overwrite=True
                )
            )

            mock_replace.assert_called_once_with(
                temp_file_path + ".temp", temp_file_path
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
        otl_taxa_ex = {"taxa1"}

        test_dir = os.path.join(os.path.dirname(__file__), "test_data")
        os.makedirs(test_dir, exist_ok=True)

        input_file_path = os.path.join(test_dir, "test_input.fasta")
        with open(input_file_path, "w") as f:
            f.write(test_fasta_content)

        try:
            total_count, kept_count, output_file = (
                self.handler.filter_fasta_for_species_not_in_otl(
                    input_file_path, otl_taxa_ex, overwrite=False
                )
            )
            self.assertEqual(total_count, 1)
            self.assertEqual(kept_count, 1)

            assert os.path.exists(output_file)
            assert output_file.endswith("test_input.fasta")
            assert "_otl_filtered" in output_file

            with open(output_file, "r") as f:
                content = f.read().strip()
            assert content == test_fasta_content
            assert os.path.dirname(output_file) != os.path.dirname(
                input_file_path
            ), "Output file should be in a different directory"

        finally:
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            if "output_file" in locals() and os.path.exists(output_file):
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
        self.binding.processed_primers = {
            "primer_1": pd.DataFrame({"taxon": ["A", "B"], "value": [1, 2]}),
            "primer_2": pd.DataFrame({"taxon": ["C", "D"], "value": [3, 4]}),
        }
        self.created_files = []
        self.mock_test_dir = tempfile.mkdtemp()

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
        input_df = pd.DataFrame(
            {
                "taxon": ["taxon1", "taxon2"],
                "seq_id": ["seq1", "seq2"],
                "analysis1": [10, 20],
                "analysis2": [5, 15],
            }
        )

        otl_taxa_set = {"taxon1", "taxon2", "taxon3", "taxon4"}

        expected_df = pd.DataFrame(
            {
                "taxon": ["taxon1", "taxon2", "taxon3", "taxon4"],
                "seq_id": ["seq1", "seq2", "otl-import", "otl-import"],
                "analysis1": [10, 20, 0, 0],
                "analysis2": [5, 15, 0, 0],
            }
        )

        otl_populated_df = self.binding.add_missing_otl_taxa_to_df_with_values_of_zero(
            input_df, otl_taxa_set
        )

        pd.testing.assert_frame_equal(
            otl_populated_df.sort_values(by="taxon").reset_index(drop=True),
            expected_df.sort_values(by="taxon").reset_index(drop=True),
        )

    def test_iterate_over_primer_pbs_df_with_otl_taxa(self):
        comprehensive_primer_dfs = {
            "primer_1": pd.DataFrame(
                {
                    "taxon": ["taxon1", "taxon2"],
                    "seq_id": ["seq1", "seq2"],
                    "analysis1": [10, 20],
                    "analysis2": [5, 15],
                }
            ),
            "primer_2": pd.DataFrame(
                {
                    "taxon": ["taxon3", "taxon4"],
                    "seq_id": ["seq3", "seq4"],
                    "analysis1": [30, 40],
                    "analysis2": [25, 35],
                }
            ),
        }

        otl_taxa_set = {"taxon1", "taxon2", "taxon3", "taxon4", "taxon5"}

        expected_primer_1_df = pd.DataFrame(
            {
                "taxon": ["taxon1", "taxon2", "taxon5"],
                "seq_id": ["seq1", "seq2", "otl-import"],
                "analysis1": [10, 20, 0],
                "analysis2": [5, 15, 0],
            }
        )
        expected_primer_2_df = pd.DataFrame(
            {
                "taxon": ["taxon3", "taxon4", "taxon5"],
                "seq_id": ["seq3", "seq4", "otl-import"],
                "analysis1": [30, 40, 0],
                "analysis2": [25, 35, 0],
            }
        )

        def side_effect_check(df, otl_taxa):
            if (
                df.equals(comprehensive_primer_dfs["primer_1"])
                and otl_taxa == otl_taxa_set
            ):
                return expected_primer_1_df
            elif (
                df.equals(comprehensive_primer_dfs["primer_2"])
                and otl_taxa == otl_taxa_set
            ):
                return expected_primer_2_df
            raise ValueError("Unexpected call to mocked method.")

        with patch.object(
            self.binding,
            "add_missing_otl_taxa_to_df_with_values_of_zero",
            side_effect=side_effect_check,
        ) as mock_add_otl_taxa:
            processed_primers = self.binding.iterate_over_primer_pbs_df(
                comprehensive_primer_dfs, add_otl_taxa=True, otl_taxa_set=otl_taxa_set
            )

            self.assertIn("primer_1", processed_primers.keys())
            self.assertIn("primer_2", processed_primers.keys())

            pd.testing.assert_frame_equal(
                processed_primers["primer_1"]
                .sort_values(by="taxon")
                .reset_index(drop=True),
                expected_primer_1_df.sort_values(by="taxon").reset_index(drop=True),
            )
            pd.testing.assert_frame_equal(
                processed_primers["primer_2"]
                .sort_values(by="taxon")
                .reset_index(drop=True),
                expected_primer_2_df.sort_values(by="taxon").reset_index(drop=True),
            )

    def test_iterate_over_primer_pbs_df_without_otl_taxa(self):
        comprehensive_primer_dfs = {
            "primer_1": pd.DataFrame(
                {
                    "taxon": ["taxon1", "taxon2"],
                    "seq_id": ["seq1", "seq2"],
                    "analysis1": [10, 20],
                    "analysis2": [5, 15],
                }
            ),
            "primer_2": pd.DataFrame(
                {
                    "taxon": ["taxon3", "taxon4"],
                    "seq_id": ["seq3", "seq4"],
                    "analysis1": [30, 40],
                    "analysis2": [25, 35],
                }
            ),
        }

        processed_primers = self.binding.iterate_over_primer_pbs_df(
            comprehensive_primer_dfs,
            add_otl_taxa=False,
        )

        self.assertIn("primer_1", processed_primers)
        self.assertIn("primer_2", processed_primers)

        pd.testing.assert_frame_equal(
            processed_primers["primer_1"]
            .sort_values(by="taxon")
            .reset_index(drop=True),
            comprehensive_primer_dfs["primer_1"]
            .sort_values(by="taxon")
            .reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            processed_primers["primer_2"]
            .sort_values(by="taxon")
            .reset_index(drop=True),
            comprehensive_primer_dfs["primer_2"]
            .sort_values(by="taxon")
            .reset_index(drop=True),
        )

    def test_get_primer_df_exists(self):
        result = self.binding.get_primer_df("primer_1")
        expected = pd.DataFrame({"taxon": ["A", "B"], "value": [1, 2]})
        pd.testing.assert_frame_equal(result, expected)

    def test_get_primer_df_not_exists(self):
        result = self.binding.get_primer_df("primer_3")
        self.assertIsNone(result)

    def test_process_analysis_per_taxon(self):
        primer_df = pd.DataFrame(
            {
                "seq-id": ["abc", "def", "gh"],
                "taxon": ["A", "A", "B"],
                "value": [5, 10, 20],
            }
        )
        # Test mean operation
        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="mean", analysis_name="value"
        )
        expected = pd.DataFrame({"value": {"A": 7.5, "B": 20.0}})
        expected.index.name = "taxon"
        pd.testing.assert_frame_equal(result, expected)

        # Test min operation
        result_min = self.binding.process_analysis_per_taxon(
            primer_df, operation="min", analysis_name="value"
        )
        expected_min = pd.DataFrame({"value": {"A": 5.0, "B": 20.0}})
        expected_min.index.name = "taxon"
        pd.testing.assert_frame_equal(result_min, expected_min)

        # Test max operation
        result_max = self.binding.process_analysis_per_taxon(
            primer_df, operation="max", analysis_name="value"
        )
        expected_max = pd.DataFrame({"value": {"A": 10.0, "B": 20.0}})
        expected_max.index.name = "taxon"
        pd.testing.assert_frame_equal(result_max, expected_max)

        # Test sum operation
        result_sum = self.binding.process_analysis_per_taxon(
            primer_df, operation="sum", analysis_name="value"
        )
        expected_sum = pd.DataFrame({"value": {"A": 15.0, "B": 20.0}})
        expected_sum.index.name = "taxon"
        pd.testing.assert_frame_equal(result_sum, expected_sum)

        # Test coef_var operation
        result_coef_var = self.binding.process_analysis_per_taxon(
            primer_df, operation="coef_var", analysis_name="value"
        )
        expected_coef_var = pd.DataFrame({"value": {"A": 47.14, "B": np.nan}})
        expected_coef_var.index.name = "taxon"
        pd.testing.assert_frame_equal(result_coef_var, expected_coef_var)

    def test_process_analysis_per_taxon_single_row(self):
        primer_df = pd.DataFrame({"seq-id": ["abc"], "taxon": ["A"], "value": [5]})
        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="mean", analysis_name="value"
        )
        expected = pd.DataFrame({"value": {"A": 5.0}})
        expected.index.name = "taxon"
        pd.testing.assert_frame_equal(result, expected)

    def test_process_analysis_per_taxon_missing_values(self):
        primer_df = pd.DataFrame(
            {
                "seq-id": ["abc", "def", "gh"],
                "taxon": ["A", "A", "B"],
                "value": [5, None, 20],
            }
        )
        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="coef_var", analysis_name="value"
        )
        expected = pd.DataFrame({"value": {"A": np.nan, "B": np.nan}})
        expected.index.name = "taxon"
        pd.testing.assert_frame_equal(result, expected)

    def test_process_analysis_per_taxon_invalid_taxon(self):
        primer_df = pd.DataFrame({"value": [10, 20]})
        with self.assertRaises(ValueError):
            self.binding.process_analysis_per_taxon(
                primer_df, operation="mean", analysis_name="value"
            )

    def test_process_analysis_per_taxon_invalid_analysis_name(self):
        primer_df = pd.DataFrame(
            {"seq-id": ["abc", "def"], "taxon": ["A", "B"], "value": [10, 20]}
        )
        with self.assertRaises(ValueError):
            self.binding.process_analysis_per_taxon(
                primer_df, operation="mean", analysis_name="missing_column"
            )

    def test_process_analysis_across_taxon_sum(self):
        tax_grouped_df = pd.DataFrame({"value": [10, 20]})
        result = self.binding.process_analysis_across_taxon(
            tax_grouped_df, operation="sum"
        )
        self.assertEqual(result, 30)

    def test_process_analysis_across_taxon_coef_var(self):
        tax_grouped_df = pd.DataFrame({"value": [30, 6]})
        result = self.binding.process_analysis_across_taxon(
            tax_grouped_df, operation="coef_var"
        )
        self.assertEqual(result, 94.28)

    def test_process_analysis_across_taxon_empty_df(self):
        tax_grouped_df = pd.DataFrame(columns=["value"])
        with self.assertRaises(ValueError):
            self.binding.process_analysis_across_taxon(tax_grouped_df, operation="sum")

    def test_process_analysis_across_taxon_invalid_operation(self):
        tax_grouped_df = pd.DataFrame({"value": [10, 20]})
        with self.assertRaises(ValueError):
            self.binding.process_analysis_across_taxon(
                tax_grouped_df, operation="invalid_op"
            )

    @patch("src.marker_scoring.scoring_utils.calculate_iupac_mismatches")
    @patch("Bio.SeqUtils.MeltingTemp.Tm_GC", return_value=60.0)
    @patch("Bio.SeqUtils.gc_fraction", return_value=0.5)
    @patch(
        "src.marker_scoring.metrics_system.Binding.get_primer_table",
        return_value=pd.DataFrame([...]),
    )
    @patch(
        "src.marker_scoring.metrics_system.Binding.parse_files_with_same_extension_in_folders",
        return_value=[...],
    )
    @patch(
        "src.marker_scoring.metrics_system.Binding.get_pbs_table",
        return_value=pd.DataFrame([...]),
    )
    def test_primer_pbs_analysis(
        self,
        mock_get_pbs_table,
        mock_parse_files,
        mock_get_primer_table,
        mock_gc_fraction,
        mock_tm_gc,
        mock_calculate_mismatches,
    ):
        mock_primer_table = pd.DataFrame(
            [
                {
                    "barcode_region": "COI",
                    "assay_name": "TestAssay",
                    "fwd_seq": "AGCTTAGCTA",
                    "rev_seq": "TCGATCGATC",
                }
            ]
        )
        mock_get_primer_table.return_value = mock_primer_table

        mock_parse_files.return_value = [("amplicon_file.fasta", "insert_file.fasta")]

        mock_pbs_table = pd.DataFrame(
            [{"header": ">seq1|taxon1", "fwd_seq": "AGCTT", "rev_seq": "TCGAT"}]
        )
        mock_get_pbs_table.return_value = mock_pbs_table

        mock_gc_fraction.side_effect = lambda seq: len(
            [c for c in seq if c in "GC"]
        ) / len(seq)
        mock_tm_gc.side_effect = lambda seq, valueset, strict: 60.0
        mock_calculate_mismatches.side_effect = (
            lambda seq1, seq2, search_gc_clamp=False: (
                1 if seq1 != seq2 else 0,
                2 if search_gc_clamp else 0,
            )
        )

        primer_pbs_df, primer_gc_df = self.binding.primer_pbs_analysis(
            "mock_amplicon_folder",
            "mock_insert_folder",
            "mock_primer_table",
            save_results=True,
        )

        self.assertIsInstance(primer_pbs_df, dict)
        self.assertIsInstance(primer_gc_df, pd.DataFrame)

        self.assertEqual(primer_gc_df.iloc[0]["barcode_region"], "COI")
        self.assertEqual(primer_gc_df.iloc[0]["assay_name"], "TestAssay")

    @patch(
        "src.in_silico_analysis.amplification.InSilicoAmplification.validate_primer_table"
    )
    @patch("src.marker_scoring.scoring_utils.calculate_iupac_mismatches")
    @patch("Bio.SeqUtils.MeltingTemp.Tm_GC")
    @patch("Bio.SeqUtils.gc_fraction")
    @patch("src.marker_scoring.metrics_system.Binding.get_primer_table")
    @patch(
        "src.marker_scoring.metrics_system.Binding.parse_files_with_same_extension_in_folders"
    )
    @patch("src.marker_scoring.metrics_system.Binding.get_pbs_table")
    def test_primer_pbs_analysis_single_primer(
        self,
        mock_get_pbs_table,
        mock_parse_files,
        mock_get_primer_table,
        mock_gc_fraction,
        mock_tm_gc,
        mock_calculate_mismatches,
        mock_validate_primer_table,
    ):
        mock_validate_primer_table.return_value = None

        mock_primer_table = pd.DataFrame(
            [
                {
                    "barcode_region": "COI",
                    "assay_name": "TestAssay",
                    "fwd_seq": "AGCTTAGCTA",
                    "rev_seq": "TCGATCGATC",
                }
            ]
        )
        mock_get_primer_table.return_value = mock_primer_table

        mock_parse_files.return_value = [
            ("COI_TestAssay.fasta", "COI_TestAssay_insert.fasta")
        ]

        mock_pbs_table = pd.DataFrame(
            [
                {"header": ">seq1|taxon1", "fwd_seq": "AGCTT", "rev_seq": "TCGAT"},
                {"header": ">seq2|taxon2", "fwd_seq": "AGCTA", "rev_seq": "TCGAC"},
            ]
        )
        mock_get_pbs_table.return_value = mock_pbs_table

        mock_gc_fraction.side_effect = lambda seq: len(
            [c for c in seq if c in "GC"]
        ) / len(seq)
        mock_tm_gc.return_value = 60.0
        mock_calculate_mismatches.side_effect = (
            lambda seq1, seq2, search_gc_clamp=False: (
                (1 if seq1 != seq2 else 0, 2 if search_gc_clamp else 0)
            )
        )

        primer_pbs_df, primer_gc_df = self.binding.primer_pbs_analysis(
            "mock_amplicon_folder",
            "mock_insert_folder",
            "mock_primer_table",
            save_results=True,
        )

        self.assertIsInstance(primer_pbs_df, dict)
        self.assertIsInstance(primer_gc_df, pd.DataFrame)

        self.assertIn("COI_TestAssay", primer_pbs_df)

        comprehensive_df = primer_pbs_df["COI_TestAssay"]
        self.assertEqual(len(comprehensive_df), 2)
        self.assertTrue(
            all(
                col in comprehensive_df.columns
                for col in [
                    "seq_id",
                    "taxon",
                    "full_len_mismatch_sum",
                    "three_end_mismatch_sum",
                    "gc_matches_fwd",
                    "gc_matches_rev",
                    "min_tm",
                    "delta_tm",
                ]
            )
        )

        self.created_files.extend(
            ["COI_TestAssay_comprehensive.csv", "primer_gc_fractions.csv"]
        )

    def test_get_priming_ratio(self):
        max_mismatch_full_len = pd.DataFrame(
            {"full_len_mismatch_sum": [10, 20]}, index=["A", "B"]
        )
        max_mismatch_three_end = pd.DataFrame(
            {"three_end_mismatch_sum": [5, 15]}, index=["A", "B"]
        )
        result = self.binding.get_priming_ratio(
            max_mismatch_full_len, max_mismatch_three_end
        )
        self.assertEqual(result, 1.25)

    def test_get_total_gc_matches(self):
        primer_pbs_df = pd.DataFrame(
            {
                "gc_matches_fwd": [1, 3, 4],
                "gc_matches_rev": [0, 3, 2],
            }
        )
        result = self.binding.get_total_gc_matches(primer_pbs_df)
        expected = pd.DataFrame(
            {
                "gc_matches_fwd": [1, 1, 0],
                "gc_matches_rev": [0, 1, 1],
                "gc_matches_score": [1, 2, 1],
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_tm_score(self):
        test_data = {"delta_tm": [1.5, 2.0, 1.8, 3.0, 1.2, 0.5]}
        primer_pbs_df = pd.DataFrame(test_data)

        temp_threshold = 2.0
        result = self.binding.tm_score(primer_pbs_df, temp_threshold=temp_threshold)

        expected = 66.67  # 4 out of 6 entries

        self.assertEqual(result, expected)

    def create_mock_fasta(self, filename, records):
        """
        Create a mock FASTA file with given records

        Args:
            filename (str): Name of the file to create
            records (list): List of SeqRecord objects

        Returns:
            str: Full path to the created FASTA file
        """
        filepath = os.path.join(self.mock_test_dir, filename)
        SeqIO.write(records, filepath, "fasta")

        return filepath

    def test_count_unique_taxa(self):
        """
        Method to test the count_unique_taxa method.
        """
        records = [
            SeqRecord(Seq("ATCG"), id="seq1", description="seq1|Taxon1"),
            SeqRecord(Seq("GCTA"), id="seq2", description="seq2|Taxon1"),
            SeqRecord(Seq("TAGC"), id="seq3", description="seq3|Taxon2"),
            SeqRecord(Seq("CGAT"), id="seq4", description="seq4|Taxon3"),
        ]

        fasta_file = self.create_mock_fasta("test_taxa.fasta", records)
        result = self.binding.count_unique_taxa(fasta_file)
        self.assertEqual(result, 3)

    def test_count_unique_taxa_no_taxa(self):
        """
        Method to test the count_unique_taxa method when there is no taxa.
        """
        records = [
            SeqRecord(Seq("ATCG"), id="seq1", description="seq1"),
            SeqRecord(Seq("GCTA"), id="seq2", description="seq2"),
        ]

        fasta_file = self.create_mock_fasta("test_no_taxa.fasta", records)

        result = self.binding.count_unique_taxa(fasta_file)
        self.assertEqual(result, 0)

    def test_get_outputs_taxa_counts_amplification_success_score(self):
        """
        Method to test the the output generated by get_outputs_taxa_counts and the
        calculate_amplification_success_score.
        """
        results_folder = self.mock_test_dir
        os.makedirs(os.path.join(results_folder, "all_inserts"))
        os.makedirs(os.path.join(results_folder, "all_complete_pbs/filtered"))
        os.makedirs(os.path.join(results_folder, "incomplete_pbs/filtered"))

        in_silico_records = [
            SeqRecord(Seq("ATCG"), id="seq1", description="seq1|Taxon1"),
            SeqRecord(Seq("GCTA"), id="seq2", description="seq2|Taxon2"),
        ]
        pbs_records = [
            SeqRecord(Seq("TAGC"), id="seq3", description="seq3|Taxon3"),
            SeqRecord(Seq("CGAT"), id="seq4", description="seq4|Taxon1"),
        ]
        incomplete_pbs_records = [
            SeqRecord(Seq("AAAA"), id="seq5", description="seq5|Taxon4"),
        ]

        self.create_mock_fasta("all_inserts/primer1.fasta", in_silico_records)
        self.create_mock_fasta("all_complete_pbs/filtered/primer1.fasta", pbs_records)
        self.create_mock_fasta(
            "incomplete_pbs/filtered/primer1.fasta", incomplete_pbs_records
        )

        result = self.binding.calculate_amplification_success_score(results_folder)

        self.assertIn("taxa_in_silico_amplified", result.columns)
        self.assertIn("taxa_with_pbs", result.columns)
        self.assertIn("taxa_with_incomplete_pbs", result.columns)

        primer1_row = result.loc[result.index == "primer1"]
        self.assertEqual(primer1_row["taxa_in_silico_amplified"].values[0], 2)
        self.assertEqual(primer1_row["taxa_with_incomplete_pbs"].values[0], 1)
        self.assertEqual(primer1_row["taxa_with_pbs"].values[0], 4)

        self.assertIn("amplification_sucess_percent", result.columns)
        self.assertEqual(
            primer1_row["amplification_sucess_percent"].values[0], float(50.0)
        )

    def test_calculate_amplification_success_score(self):
        """
        Method to test the calculate_amplification_success_score.
        """
        results_folder = self.mock_test_dir

        os.makedirs(os.path.join(results_folder, "all_inserts"))
        os.makedirs(os.path.join(results_folder, "all_complete_pbs/filtered"))
        os.makedirs(os.path.join(results_folder, "incomplete_pbs/filtered"))

        result = self.binding.calculate_amplification_success_score(results_folder)

        self.assertIn("amplification_sucess_percent", result.columns)

    def tearDown(self):
        for file_path in self.created_files:
            if os.path.exists(file_path):
                os.remove(file_path)
        self.created_files.clear()

        shutil.rmtree(self.mock_test_dir)


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
