"""
Unit tests for metrics_system.py
"""

import shutil
import sys
import tempfile
import unittest
from pathlib import Path
import os
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, call, mock_open, patch, Mock

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

        expected_output = pd.DataFrame({
            "barcoded_taxa_one_plus": [50.0],
            "ratio_barcoded_taxa": [0.33],
        }, index=["test_amplicon_reffb"])

        print(rbt_rounded)


        pd.testing.assert_frame_equal(rbt_rounded, expected_output)


class TestBinding(unittest.TestCase):
    def setUp(self):
        self.otl = "data/test_data/test_otl.tsv"
        self.binding = Binding(self.otl)
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
            binding = Binding(self.otl)
            mock_get_mismatches.assert_called_once()
            self.assertEqual(binding.number_of_mismatches, 2)

    def test_init_with_custom_mismatches(self):
        """
        Test Binding class initialization with custom number of mismatches.
        """
        binding = Binding(self.otl, number_of_mismatches=3)
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
            "mock_primer_table"
        )

        self.assertIsInstance(primer_pbs_df, dict)
        self.assertIsInstance(primer_gc_df, pd.DataFrame)

        self.assertEqual(primer_gc_df.index.name, "primer_name")
        self.assertEqual(primer_gc_df.index[0], "COI_TestAssay")

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
            "mock_primer_table"
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
        os.makedirs(os.path.join(results_folder, "insert/filtered"))
        os.makedirs(os.path.join(results_folder, "all_complete_pbs/filtered/filtered_intersection"))
        os.makedirs(os.path.join(results_folder, "incomplete_pbs/filtered/filtered_intersection"))

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

        self.create_mock_fasta("insert/filtered/primer1.fasta", in_silico_records)
        self.create_mock_fasta("all_complete_pbs/filtered/filtered_intersection/primer1.fasta", pbs_records)
        self.create_mock_fasta(
            "incomplete_pbs/filtered/filtered_intersection/primer1.fasta", incomplete_pbs_records
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

class TestTraitsAndResolution(unittest.TestCase):
    def setUp(self):
        """
        Set up test parameters using existing directories
        """
        self.test_dir = "data/test_data"
        self.amplicon_dir = self.test_dir + "/amplicon-test"
        self.insert_dir = self.test_dir + "/insert-test"
        self.incomplete_pbs_dir= self.test_dir + "/insert-test"
        self.multibarcodetools_input = self.test_dir + "/insert-test/multibarcodetools-input.tsv"
        self.traits = TraitsAndResolution(
            insert_folder_path=self.insert_dir,
            amplicon_folder_path=self.amplicon_dir,
            incomplete_pbs_folder_path=self.incomplete_pbs_dir
        )
        self.created_files = [self.multibarcodetools_input]
        self.traits.multibarcode_output_folder = os.path.join(os.path.dirname(self.traits.insert_folder_path), 'multibarcode')
        self.primer_resolv_species = pd.DataFrame({
            'primer': ['primerA', 'primerB', 'primerC'],
            'additional_resolved_species': [90, 8, 2],
            'cumulative_resolved_species': [90, 98, 100]
        })
        self.traits.primer_resolv_species = self.primer_resolv_species

    def test_init_with_results_folder(self):
        results_folder = "path/to/results"
        expected_insert_path = os.path.join(results_folder, 'insert/filtered')
        expected_amplicon_path = os.path.join(results_folder, 'amplicon/filtered')
        expected_incomplete_pbs_path = os.path.join(results_folder, 'incomplete_pbs/filtered/filtered_intersection')

        traits = TraitsAndResolution(results_folder=results_folder)

        self.assertEqual(traits.insert_folder_path, expected_insert_path)
        self.assertEqual(traits.amplicon_folder_path, expected_amplicon_path)
        self.assertEqual(traits.incomplete_pbs_path, expected_incomplete_pbs_path)

    def test_init_with_individual_paths(self):
        insert_path = "path/to/insert"
        amplicon_path = "path/to/amplicon"
        incomplete_pbs_path = "path/to/incomplete_pbs"

        traits = TraitsAndResolution(
            insert_folder_path=insert_path,
            amplicon_folder_path=amplicon_path,
            incomplete_pbs_folder_path=incomplete_pbs_path
        )

        self.assertEqual(traits.insert_folder_path, insert_path)
        self.assertEqual(traits.amplicon_folder_path, amplicon_path)
        self.assertEqual(traits.incomplete_pbs_path, incomplete_pbs_path)

    def test_init_with_missing_arguments(self):
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution()
        self.assertIn("Either provide a path to the in-silico amplification results folder", str(context.exception))

    def test_init_with_partial_arguments(self):
        # Test missing amplicon_folder_path
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution(
                insert_folder_path="path/to/insert",
                incomplete_pbs_folder_path="path/to/incomplete_pbs"
            )
        self.assertIn("Either provide a path to the in-silico amplification results folder", str(context.exception))

        # Test missing insert_folder_path
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution(
                amplicon_folder_path="path/to/amplicon",
                incomplete_pbs_folder_path="path/to/incomplete_pbs"
            )
        self.assertIn("Either provide a path to the in-silico amplification results folder", str(context.exception))

    def test_get_min_max_avg_seq_length_in_a_fasta(self):
        """
        Test sequence length calculation for individual FASTA files
        """
        amplicon_primer_a_path = os.path.join(self.amplicon_dir, 'primerA.fasta')
        min_len, max_len, avg_len = self.traits.get_min_max_avg_seq_length_in_a_fasta(amplicon_primer_a_path)

        self.assertEqual(min_len, 6)
        self.assertEqual(max_len, 26)
        self.assertAlmostEqual(avg_len, 18.33)

    def test_get_length_stats_for_amplicon_and_insert(self):
        """
        Test getting length statistics for both amplicon and insert
        """
        length_stats = self.traits.get_length_stats_for_amplicon_and_insert()

        self.assertIsInstance(length_stats, pd.DataFrame)

        self.assertListEqual(list(length_stats.index), ['primerA', 'primerC'])

        # Insert average length
        self.assertAlmostEqual(length_stats.loc['primerA', 'insert_avg_length'], 9.0)
        self.assertTrue(np.isnan(length_stats.loc['primerC', 'insert_avg_length']))

        # Amplicon min length
        self.assertAlmostEqual(length_stats.loc['primerA', 'amplicon_min_length'], 6)
        self.assertAlmostEqual(length_stats.loc['primerC', 'amplicon_min_length'], 23)

        # Amplicon max length
        self.assertAlmostEqual(length_stats.loc['primerA', 'amplicon_max_length'], 26)
        self.assertAlmostEqual(length_stats.loc['primerC', 'amplicon_max_length'], 26)

        # Amplicon average length
        self.assertAlmostEqual(length_stats.loc['primerA', 'amplicon_avg_length'], 18.33, places=2)
        self.assertAlmostEqual(length_stats.loc['primerC', 'amplicon_avg_length'], 24.50, places=2)

        expected_columns = [
            'amplicon_min_length',
            'amplicon_max_length',
            'amplicon_avg_length',
            'insert_avg_length'
        ]
        for column in expected_columns:
            self.assertIn(column, length_stats.columns, f"{column} column missing")

    @patch('subprocess.run')
    def test_run_multibarcode_pipeline_success(self, mock_subprocess_run):
        mock_result = MagicMock()
        mock_result.stdout = """
        1: primerA, resolve  3 species
        2: primerC, resolve additional 1 species
        """
        mock_subprocess_run.return_value = mock_result

        result = self.traits.run_multibarcode_pipeline()

        mock_subprocess_run.assert_called_once()
        self.assertIsNotNone(result)

        expected_multibarcode_output_folder = os.path.join(self.test_dir, 'multibarcode')
        self.assertTrue(os.path.exists(expected_multibarcode_output_folder))

    def test_parse_multibarcode_output(self):
        test_stdout = """
        1: primerA, resolve  3 species
        2: primerC, resolve additional 1 species
        3: primerD, resolve additional 1 species
        """

        result_df = self.traits.parse_multibarcode_output(test_stdout)

        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(len(result_df), 3)

        expected_columns = [
            'primer',
            'additional_resolved_species',
            'cumulative_resolved_species'
        ]
        self.assertListEqual(list(result_df.columns), expected_columns)

        self.assertEqual(result_df.iloc[0]['primer'], 'primerA')
        self.assertEqual(result_df.iloc[0]['additional_resolved_species'], 3)
        self.assertEqual(result_df.iloc[0]['cumulative_resolved_species'], 3)

        self.assertEqual(result_df.iloc[1]['primer'], 'primerC')
        self.assertEqual(result_df.iloc[1]['additional_resolved_species'], 1)
        self.assertEqual(result_df.iloc[1]['cumulative_resolved_species'], 4)

    def test_parse_multibarcode_output_empty(self):
        result_df = self.traits.parse_multibarcode_output("")

        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(len(result_df), 0)

    @patch('subprocess.run')
    def test_run_multibarcode_pipeline_failure(self, mock_subprocess_run):
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd=['multi-barcode']
        )

        with self.assertRaises(subprocess.CalledProcessError):
            self.traits.run_multibarcode_pipeline()

    # def test_get_taxonomic_resolution(self):
    #     tax_rex_df = self.traits.get_taxonomic_resolution(otl_total_taxa_count=100)

    #     self.assertIn('taxonomic_resolution_percentage', tax_rex_df.columns)

    #     self.assertEqual(tax_rex_df.iloc[0]['taxonomic_resolution_percentage'], 90.0)
    #     self.assertEqual(tax_rex_df.iloc[1]['taxonomic_resolution_percentage'], 98.0)
    #     self.assertEqual(tax_rex_df.iloc[-1]['taxonomic_resolution_percentage'], 100.0)

    @patch("pandas.read_excel")
    def test_load_nucleotide_distance(self, mock_read_excel):
        mock_df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        mock_read_excel.return_value = mock_df

        result = self.traits.load_nucleotide_distance()

        mock_read_excel.assert_called_once_with("data/test_data/multibarcode/matrix.xlsx")
        pd.testing.assert_frame_equal(result, mock_df)

    @patch.object(TraitsAndResolution, "load_nucleotide_distance")
    @patch.object(TraitsAndResolution, "get_length_stats_for_amplicon_and_insert")
    def test_compute_genetic_divergence_per_taxon(self, mock_get_length_stats, mock_load_nucleotide_distance):
        # Mock data for load_nucleotide_distance
        mock_nuc_dist_matrix = pd.DataFrame({
            'Species': ['Species1', 'Species2'],
            'primerA': [10, 20],
            'primerB': [30, 40]
        })
        mock_load_nucleotide_distance.return_value = mock_nuc_dist_matrix

        # Mock data for get_length_stats_for_amplicon_and_insert
        mock_length_stats = {
            'insert_avg_length': pd.Series([100, 200], index=['primerA', 'primerB'])
            }
        mock_get_length_stats.return_value = mock_length_stats

        result = self.traits.compute_genetic_divergence_per_taxon()

        expected_result = pd.DataFrame({
            'Species': ['Species1', 'Species2'],
            'primerA': [10.0, 20.0],  # (10/100)*100, (20/100)*100
            'primerB': [15.0, 20.0]   # (30/200)*100, (40/200)*100
        })

        pd.testing.assert_frame_equal(result, expected_result)

        mock_load_nucleotide_distance.assert_called_once()
        mock_get_length_stats.assert_called_once()

    @patch.object(TraitsAndResolution, "load_nucleotide_distance")
    @patch.object(TraitsAndResolution, "get_length_stats_for_amplicon_and_insert")
    def test_compute_genetic_divergence_per_taxon_nan_values(self, mock_get_length_stats, mock_load_nucleotide_distance):
        # Mock data for load_nucleotide_distance
        mock_nuc_dist_matrix = pd.DataFrame({
            'Species': ['Species1', 'Species2'],
            'primerA': [10, 20],
            'primerB': ['-', 40]
        })
        mock_load_nucleotide_distance.return_value = mock_nuc_dist_matrix

        # Mock data for get_length_stats_for_amplicon_and_insert
        mock_length_stats = {
            'insert_avg_length': pd.Series([100, 0], index=['primerA', 'primerB'])
            }
        mock_get_length_stats.return_value = mock_length_stats

        result = self.traits.compute_genetic_divergence_per_taxon()

        expected_result = pd.DataFrame({
            'Species': ['Species1', 'Species2'],
            'primerA': [10.0, 20.0],  # (10/100)*100, no length insert results in NaN
            'primerB': [np.nan, np.nan]   # 'invalid' results in NaN, no length insert results in NaN
        })

        pd.testing.assert_frame_equal(result, expected_result)

        mock_load_nucleotide_distance.assert_called_once()
        mock_get_length_stats.assert_called_once()

    def test_get_divergence_score_normal_case(self):
        mock_divergence_df = pd.DataFrame({
            'Species': ['Species1', 'Species2', 'Species3', 'Species4'],
            'primerA': [1.5, 2.5, 3.0, 1.0],
            'primerB': [1.0, 1.5, 2.0, 3.0]
        })

        with patch.object(self.traits, 'compute_genetic_divergence_per_taxon',
                        return_value=mock_divergence_df):

            result = self.traits.get_divergence_score(4)

            self.assertIsInstance(result, pd.DataFrame)
            self.assertEqual(len(result), 2)

            expected_columns = ['primer', 'total_taxa', 'n_taxa_above_cutoff', 'divergence_score']
            self.assertListEqual(list(result.columns), expected_columns)

    def test_get_divergence_score_zero_total_taxa(self):
        with self.assertRaises(ValueError):
            self.traits.get_divergence_score(0)

    def test_get_divergence_score_negative_total_taxa(self):
        with self.assertRaises(ValueError):
            self.traits.get_divergence_score(-5)

    def test_get_divergence_score_custom_cutoff(self):
        mock_divergence_df = pd.DataFrame({
            'Species': ['Species1', 'Species2', 'Species3', 'Species4'],
            'primerA': [1.5, 2.5, 3.0, 1.0],
            'primerB': [1.0, 1.5, 2.0, 3.0]
        })

        with patch.object(self.traits, 'compute_genetic_divergence_per_taxon',
                        return_value=mock_divergence_df):

            result = self.traits.get_divergence_score(4, cutoff=1.8)

            for _, row in result.iterrows():
                expected_taxa_above_cutoff = sum(
                    mock_divergence_df[row['primer']] > 1.8
                )
                self.assertEqual(row['n_taxa_above_cutoff'], expected_taxa_above_cutoff)

    def test_get_divergence_score_no_taxa_above_cutoff(self):
        mock_divergence_df = pd.DataFrame({
            'Species': ['Species1', 'Species2'],
            'primerA': [1.0, 1.5],
            'primerB': [0.5, 1.0]
        })

        with patch.object(self.traits, 'compute_genetic_divergence_per_taxon',
                        return_value=mock_divergence_df):

            result = self.traits.get_divergence_score(2, cutoff=2.0)

            self.assertTrue((result['divergence_score'] == 0).all())

    def test_get_divergence_score_dataframe_structure(self):
        mock_divergence_df = pd.DataFrame({
            'Species': ['Species1', 'Species2', 'Species3', 'Species4'],
            'primerA': [1.5, 2.5, 3.0, 1.0],
            'primerB': [1.0, 1.5, 2.0, 3.0]
        })

        with patch.object(self.traits, 'compute_genetic_divergence_per_taxon',
                        return_value=mock_divergence_df):

            result = self.traits.get_divergence_score(4)

            self.assertTrue(result['primer'].is_unique)

            self.assertTrue(((result['divergence_score'] >= 0) &
                            (result['divergence_score'] <= 100)).all())

    def tearDown(self):
        for file_path in self.created_files:
            if os.path.exists(file_path):
                os.remove(file_path)
        self.created_files.clear()

        if os.path.exists(self.traits.multibarcode_output_folder):
            shutil.rmtree(self.traits.multibarcode_output_folder)

class TestMetricsSystemExecutor(unittest.TestCase):
    def setUp(self):
        """
        Set up test fixtures before each test method.
        """
        self.results_folder = "/fake/path/results"
        self.otl = "/fake/path/otl.tsv"
        self.primer_table = "/fake/path/primer_table.tsv"

        # Mock OTLHandler
        self.mock_otl_handler = Mock()
        self.mock_otl_handler.total_taxa = 100
        self.mock_otl_handler.otl_taxa_set = set(['species1', 'species2'])

        # Set up path patches
        self.mock_paths = {
            'insert': '/fake/path/results/insert/filtered',
            'amplicon': '/fake/path/results/amplicon/filtered',
            'incomplete_pbs': '/fake/path/results/incomplete_pbs/filtered/filtered_intersection'
        }

    @patch('os.path.join')
    @patch('src.marker_scoring.metrics_system.OtlHandler')
    def test_initialization(self, mock_otl_handler_class, mock_join):
        """
        Test the initialization of MetricsSystemExecutor.
        """
        mock_otl_handler_instance = Mock()
        mock_otl_handler_class.return_value = mock_otl_handler_instance
        mock_join.side_effect = lambda *args: '/'.join(args)

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table
        )

        self.assertEqual(executor.results_folder, self.results_folder)
        self.assertEqual(executor.otl, self.otl)
        self.assertEqual(executor.primer_table, self.primer_table)
        mock_otl_handler_instance.import_otl.assert_called_once()

    @patch('src.marker_scoring.metrics_system.ReferenceDatabaseQuality')
    @patch('os.path.exists')
    @patch('os.path.join')
    @patch('src.marker_scoring.metrics_system.OtlHandler')
    def test_get_reference_database_quality(self, mock_otl_handler_class, mock_join, mock_exists, mock_ref_db_class):
        """
        Test reference database quality calculation.
        """
        mock_exists.return_value = True
        mock_join.side_effect = lambda *args: '/'.join(args)

        # Mock OtlHandler
        mock_otl_instance = Mock()
        mock_otl_instance.total_taxa = 100
        mock_otl_instance.otl_taxa_set = set(['species1', 'species2'])
        mock_otl_handler_class.return_value = mock_otl_instance

        # set ReferenceDatabaseQuality mock
        mock_ref_db_instance = Mock()
        mock_ref_db_class.return_value = mock_ref_db_instance
        expected_result = pd.DataFrame({
            'barcoded_taxa_one_plus': [80.0],
            'ratio_barcoded_taxa': [0.8]
        })
        mock_ref_db_instance.barcoded_taxa_ratio.return_value = expected_result

        with patch('builtins.open', create=True):
            executor = MetricsSystemExecutor(
                results_folder=self.results_folder,
                otl=self.otl,
                primer_table=self.primer_table
            )

        result = executor.get_reference_database_quality()

        # Verify results
        pd.testing.assert_frame_equal(result, expected_result)
        expected_calls = [
            call(otl='/fake/path/otl.tsv', all_inserts_path='/fake/path/results/insert/filtered'),
            call().barcoded_taxa_ratio(total_taxa_count=100)
        ]
        mock_ref_db_class.assert_has_calls(expected_calls)

        # Verify ReferenceDatabaseQuality was initialized correctly
        mock_ref_db_instance.barcoded_taxa_ratio.assert_called_once_with(
            total_taxa_count=executor.total_otl_taxa_count
        )

    @patch('src.marker_scoring.metrics_system.Binding')
    @patch('src.marker_scoring.metrics_system.OtlHandler')
    @patch('os.path.exists')
    @patch('os.path.join')
    def test_get_primer_pbs_analysis(self, mock_join, mock_exists, mock_otl_handler_class, mock_binding_class):
        """
        Test primer PBS analysis.
        """
        # Configure mock
        mock_binding_instance = Mock()
        mock_binding_class.return_value = mock_binding_instance

        expected_primer_pbs = {
            'primer1': {'data': 'value1'},
            'primer2': {'data': 'value2'}
        }
        expected_gc_df = pd.DataFrame({
            'gc_content': [0.5, 0.6]
        }, index=['primer1', 'primer2'])

        mock_binding_instance.primer_pbs_analysis.return_value = (
            expected_primer_pbs,
            expected_gc_df
        )

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table
        )

        primer_pbs_dict, gc_df = executor.get_primer_pbs_analysis()

        self.assertEqual(primer_pbs_dict, expected_primer_pbs)
        pd.testing.assert_frame_equal(gc_df, expected_gc_df)
        mock_binding_instance.primer_pbs_analysis.assert_called_once_with(
            insert_folder=executor.insert_folder_path,
            amplicon_folder=executor.amplicon_folder_path,
            primer_table=executor.primer_table
        )

    @patch('src.marker_scoring.metrics_system.TraitsAndResolution')
    @patch('src.marker_scoring.metrics_system.OtlHandler')
    @patch('os.path.exists')
    @patch('os.path.join')
    def test_get_traits_and_resolution(self,  mock_join, mock_exists, mock_otl_handler_class, mock_traits_class):
        """
        Test traits and resolution calculation.
        """
        mock_exists.return_value = True
        mock_join.side_effect = lambda *args: '/'.join(args)

        # Configure mock
        mock_otl_handler_instance = Mock()
        mock_otl_handler_instance.get_taxa_count.return_value = 100
        mock_otl_handler_class.return_value = mock_otl_handler_instance

        mock_traits_instance = Mock()
        mock_traits_class.return_value = mock_traits_instance

        # Mock expected DataFrame with 'primer' index
        expected_result = pd.DataFrame({
            'primer': ['primer1'],
            'divergence_score': [0.75]
        }).set_index('primer')
        mock_traits_instance.get_divergence_score.return_value = expected_result

        # Create instance with mocked dependencies
        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table
        )
        executor.total_otl_taxa_count = 100

        # Run test
        result = executor.get_traits_and_resolution()

        # Verify results
        pd.testing.assert_frame_equal(result, expected_result)
        mock_traits_instance.get_divergence_score.assert_called_once_with(
            total_otl_taxa_count=100
        )


    @patch('src.marker_scoring.metrics_system.OtlHandler')
    @patch('os.path.exists')
    @patch('os.path.join')
    def test_rank_primers(self, mock_join, mock_exists, mock_otl_handler_class):
        """
        Test primer ranking functionality.
        """
        mock_analysis_results = pd.DataFrame({
            'barcoded_taxa_one_plus': [90, 80],
            'ratio_barcoded_taxa': [0.9, 0.8],
            'mismatch_score': [2, 3],
            'priming_ratio_sum': [0.8, 0.7],
            'gc_matches_across_taxon': [15, 12],
            'min_tm_cv': [0.1, 0.2],
            'tm_score': [0.9, 0.8],
            'amplification_success_percent': [95, 85],
            'divergence_score': [0.2, 0.3]
        }, index=['primer1', 'primer2'])

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table
        )
        executor.join_analysis_results = Mock(return_value=mock_analysis_results)

        result = executor.rank_primers()
        result = result.set_index('index')
        print(result)

        self.assertTrue('final_rank' in result.columns)
        self.assertEqual(len(result), 2)
        primer1_rank = result.loc['primer1', 'final_rank']
        primer2_rank = result.loc['primer2', 'final_rank']
        self.assertLess(primer1_rank, primer2_rank)

