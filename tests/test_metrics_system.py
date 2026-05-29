"""
Unit tests for metrics_system.py
"""

import os
import shutil
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, Mock, call, mock_open, patch

import numpy as np
import pandas as pd

from src.mozaiko.marker_scoring.metrics_system import *
from src.mozaiko.marker_scoring.scoring_utils import *


class TestOtlHandler(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path(__file__).resolve().parent / "data/test_data"
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.handler = OtlHandler(self.otl)
        self.maxDiff = None

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

    @patch("builtins.input", return_value=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv")
    @patch("sys.stdout", new_callable=StringIO)
    def test_import_otl(self, mock_stdout, mock_input):
        otl_handler = OtlHandler()
        total_taxa_count = otl_handler.import_otl()

        self.assertIn(
            "mozaiko INFO: To continue the evaluation, a Operational Taxonomic List (OTL) is \
                    required. An OTL is a list contaning information on the taxonomic numenclature \
                        of all identifiable taxa in routine biomonitoring initiatives.",
            mock_stdout.getvalue(),
        )

        otl = pd.read_csv((Path(__file__).resolve().parent / "data/test_data/test_otl.tsv"), sep="\t", header=0)

        # Check OTL loaded correctly
        self.assertIsInstance(otl_handler.otl, pd.DataFrame)

        # Ensure expected columns exist
        expected_columns = [
            "taxa",
            "scientificName",
            "rank",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]

        self.assertEqual(list(otl_handler.otl.columns), expected_columns)

        # Validate taxa count
        self.assertEqual(total_taxa_count, 6)


class TestReferenceDatabaseQuality(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path(__file__).resolve().parent / "data/test_data"
        self.inserts_file = str(self.test_directory / "test-folder-metrics")
        self.otl = str(self.test_directory / "test_otl.tsv")
        self.ref_bd_cls = ReferenceDatabaseQuality(
            all_inserts_path=self.inserts_file, otl=self.otl
        )
        self.handler = OtlHandler(self.otl)
        self.handler.import_otl()
        self.total_otl_taxa_count = self.handler.total_taxa

    @patch("builtins.input", return_value=Path(__file__).resolve().parent / "data/test_data/test-folder-metrics")
    @patch("sys.stdout", new_callable=StringIO)
    def test_calculate_number_of_barcodes_per_taxon_input(
        self, mock_stdout, mock_input
    ):
        with patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler") as mock_otl_handler:
            mock_instance = MagicMock()
            mock_instance.taxa_hierarchy = {
                "familyA": {"genera": {"genusA": {"species": {"speciesA": {}}}}},
                "familyB": {"genera": {"genusB": {"species": {"speciesB": {}}}}},
            }
            mock_otl_handler.return_value = mock_instance

            no_file_class = ReferenceDatabaseQuality(otl=(Path(__file__).resolve().parent / "data/test_data/test_otl.tsv"))

            with patch("sys.exit") as mock_exit:
                mock_exit.side_effect = SystemExit(1)
                try:
                    barcodes_per_entry = (
                        no_file_class.calculate_number_of_barcodes_per_fasta_entry()
                    )
                    no_file_class.calculate_number_of_barcodes_per_otl_taxonomy(
                        barcodes_per_entry, no_file_class.taxa_hierarchy
                    )
                except SystemExit:
                    pass

            self.assertIn(
                "mozaiko INFO: Please provide a folder containing FASTA files.",
                mock_stdout.getvalue(),
            )

            expected = str(Path(__file__).resolve().parent / "data/test_data/test-folder-metrics")

            self.assertEqual(str(no_file_class.all_inserts_path), expected)

    def test_calculate_number_of_barcodes_per_taxon(self):
        barcodes_per_entry = (
            self.ref_bd_cls.calculate_number_of_barcodes_per_fasta_entry()
        )
        print("Barcodes per entry:", barcodes_per_entry)
        taxa_hierarchy = self.ref_bd_cls.taxa_hierarchy
        barcodes = self.ref_bd_cls.calculate_number_of_barcodes_per_otl_taxonomy(
            barcodes_per_entry, taxa_hierarchy
        )
        print("Barcodes per taxon:", barcodes)
        print(taxa_hierarchy)

        expected_barcodes = {
            "test_amplicon_reffb": {
                "familyA": {
                    "genera": {
                        "genusA": {"species": {"speciesA": {"count": 8}}, "count": 8}
                    },
                    "count": 8,
                },
                "familyB": {
                    "genera": {
                        "genusB": {"species": {"speciesB": {"count": 2}}, "count": 2}
                    },
                    "count": 2,
                },
                "familyC": {
                    "genera": {
                        "genusC": {"species": {"speciesC": {"count": 1}}, "count": 1}
                    },
                    "count": 1,
                },
                "familyD": {
                    "genera": {
                        "genusD": {"species": {"speciesD": {"count": 2}}, "count": 2}
                    },
                    "count": 2,
                },
                "familyE": {
                    "genera": {
                        "genusE": {"species": {"speciesE": {"count": 0}}, "count": 0}
                    },
                    "count": 0,
                },
                "familyZ": {
                    "genera": {
                        "genusZ": {"species": {"speciesZ": {"count": 0}}, "count": 0}
                    },
                    "count": 0,
                },
            }
        }
        self.assertEqual(barcodes, expected_barcodes)

    def test_calculate_percentage_of_taxa_w_1_barcode(self):
        percentage = self.ref_bd_cls.calculate_proportion_of_taxa_w_x_barcodes(
            self.total_otl_taxa_count, barcode_threshold=1
        )

        expected_percentage = {"test_amplicon_reffb": 0.67}

        self.assertEqual(percentage, expected_percentage)

    def test_calculate_percentage_of_taxa_w_2_barcode(self):
        percentage2 = self.ref_bd_cls.calculate_proportion_of_taxa_w_x_barcodes(
            self.total_otl_taxa_count,
            barcode_threshold=2,
        )

        expected_percentage2 = {"test_amplicon_reffb": 0.5}

        self.assertEqual(percentage2, expected_percentage2)

    def test_ratio_barcoded_taxa(self):
        rbt_rounded = self.ref_bd_cls.barcoded_taxa_ratio(self.total_otl_taxa_count)

        expected_output = pd.DataFrame(
            {
                "barcoded_taxa": [0.67],
                "ratio_barcoded_taxa": [0.25],
            },
            index=["test_amplicon_reffb"],
        )

        print(rbt_rounded)

        pd.testing.assert_frame_equal(rbt_rounded, expected_output)


class TestBinding(unittest.TestCase):
    def setUp(self):
        self.data_dir = Path(__file__).resolve().parent / "data/test_data"
        self.otl =  self.data_dir / "test_otl.tsv"
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
            "src.mozaiko.marker_scoring.scoring_utils.extract_primer_binding_sites"
        ) as mock_extract, patch("builtins.open", create=True):
            mock_pbs_table = pd.DataFrame(
                {"header": [">seq1|taxon"], "fwd_seq": ["ACGT"], "rev_seq": ["GCTA"]}
            )

            mock_extract.return_value = mock_pbs_table

            result = self.binding.get_pbs_table("amplicon.fasta", "insert.fasta")

            self.assertTrue("primer_name" in result.columns)

    def test_primer_pbs_analysis_no_matching_files(self):
        with patch.object(
            Binding, "parse_files_with_same_extension_in_folders", return_value=[]
        ) as mock_parse_files:

            otl_path = Path(__file__).resolve().parent / "data/test_data/test_primer_table.tsv"

            result = self.binding.primer_pbs_analysis(
                "amplicon_folder",
                "insert_folder",
                otl_path,
            )
            self.assertEqual(result, (None, None))

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
                "family": ["familyA", "familyA", "familyB"],
                "genus": ["genusA", "genusA", "genusB"],
                "species": ["speciesA", "speciesA", "speciesB"],
            }
        )

        # Test mean operation
        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="mean", analysis_name="value"
        )

        expected = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [7.5, 20.0, np.nan, np.nan, np.nan, np.nan],
            }
        )

        pd.testing.assert_frame_equal(result, expected)

        # Test min operation
        result_min = self.binding.process_analysis_per_taxon(
            primer_df, operation="min", analysis_name="value"
        )

        expected_min = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [5.0, 20.0, np.nan, np.nan, np.nan, np.nan],
            }
        )

        pd.testing.assert_frame_equal(result_min, expected_min)

        # Test max operation
        result_max = self.binding.process_analysis_per_taxon(
            primer_df, operation="max", analysis_name="value"
        )

        expected_max = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [10.0, 20.0, np.nan, np.nan, np.nan, np.nan],
            }
        )

        pd.testing.assert_frame_equal(result_max, expected_max)

        # Test sum operation
        result_sum = self.binding.process_analysis_per_taxon(
            primer_df, operation="sum", analysis_name="value"
        )

        expected_sum = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [15.0, 20.0, np.nan, np.nan, np.nan, np.nan],
            }
        )

        pd.testing.assert_frame_equal(result_sum, expected_sum)

        # Test coef_var operation
        result_coef_var = self.binding.process_analysis_per_taxon(
            primer_df, operation="coef_var", analysis_name="value"
        )

        expected_coef_var = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [47.14, np.nan, np.nan, np.nan, np.nan, np.nan],
            }
        )

        print(expected_coef_var)
        print(result_coef_var)

        pd.testing.assert_frame_equal(result_coef_var, expected_coef_var)

    def test_process_analysis_per_taxon_single_row(self):
        primer_df = pd.DataFrame(
            {
                "seq-id": ["abc"],
                "taxon": ["A"],
                "value": [5],
                "family": ["familyA"],
                "genus": ["genusA"],
                "species": ["speciesA"],
            }
        )

        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="mean", analysis_name="value"
        )

        expected = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [5.0, np.nan, np.nan, np.nan, np.nan, np.nan],
            }
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_process_analysis_per_taxon_missing_values(self):
        primer_df = pd.DataFrame(
            {
                "seq-id": ["abc", "def", "gh"],
                "taxon": ["A", "A", "B"],
                "value": [5, None, 20],
                "family": ["familyA", "familyA", "familyB"],
                "genus": ["genusA", "genusA", "genusB"],
                "species": ["speciesA", "speciesA", "speciesB"],
            }
        )

        result = self.binding.process_analysis_per_taxon(
            primer_df, operation="coef_var", analysis_name="value"
        )

        expected = pd.DataFrame(
            {
                "family": [
                    "familyA",
                    "familyB",
                    "familyC",
                    "familyD",
                    "familyE",
                    "familyZ",
                ],
                "genus": ["genusA", "genusB", "genusC", "genusD", "genusE", "genusZ"],
                "species": [
                    "speciesA",
                    "speciesB",
                    "speciesC",
                    "speciesD",
                    "speciesE",
                    "speciesZ",
                ],
                "value": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            }
        )

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

    @patch("src.mozaiko.marker_scoring.scoring_utils.calculate_iupac_mismatches")
    @patch("Bio.SeqUtils.MeltingTemp.Tm_GC", return_value=60.0)
    @patch("Bio.SeqUtils.gc_fraction", return_value=0.5)
    @patch(
        "src.mozaiko.marker_scoring.metrics_system.Binding.get_primer_table",
        return_value=pd.DataFrame([...]),
    )
    @patch(
        "src.mozaiko.marker_scoring.metrics_system.Binding.parse_files_with_same_extension_in_folders",
        return_value=[...],
    )
    @patch(
        "src.mozaiko.marker_scoring.metrics_system.Binding.get_pbs_table",
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
            "mock_amplicon_folder", "mock_insert_folder", "mock_primer_table"
        )

        self.assertIsInstance(primer_pbs_df, dict)
        self.assertIsInstance(primer_gc_df, pd.DataFrame)

        self.assertEqual(primer_gc_df.index.name, "primer_name")
        self.assertEqual(primer_gc_df.index[0], "COI_TestAssay")

    @patch(
        "src.mozaiko.in_silico_analysis.amplification.InSilicoAmplification.validate_primer_table"
    )
    @patch("src.mozaiko.marker_scoring.scoring_utils.calculate_iupac_mismatches")
    @patch("Bio.SeqUtils.MeltingTemp.Tm_GC")
    @patch("Bio.SeqUtils.gc_fraction")
    @patch("src.mozaiko.marker_scoring.metrics_system.Binding.get_primer_table")
    @patch(
        "src.mozaiko.marker_scoring.metrics_system.Binding.parse_files_with_same_extension_in_folders"
    )
    @patch("src.mozaiko.marker_scoring.metrics_system.Binding.get_pbs_table")
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
                {
                    "header": ">seq1|taxid1|taxon1|rank1|lin1|lin2|lin3|lin4|family1|genus1|species1",
                    "fwd_seq": "AGCTT",
                    "rev_seq": "TCGAT",
                },
                {
                    "header": ">seq2|taxid2|taxon2|rank2|lin1|lin2|lin3|lin4|family2|genus2|species2",
                    "fwd_seq": "AGCTA",
                    "rev_seq": "TCGAC",
                },
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
            "mock_amplicon_folder", "mock_insert_folder", "mock_primer_table"
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
                ]
            )
        )

        self.created_files.extend(
            ["COI_TestAssay_comprehensive.csv", "primer_gc_fractions.csv"]
        )

    def test_get_priming_ratio(self):
        max_mismatch_full_len = pd.DataFrame(
            {"full_len_mismatch_sum": [10, 20]},
            index=pd.MultiIndex.from_tuples(
                [
                    ("familyA", "genusA", "speciesA"),
                    ("familyB", "genusB", "speciesB"),
                ],
                names=["family", "genus", "species"],
            ),
        )

        max_mismatch_three_end = pd.DataFrame(
            {"three_end_mismatch_sum": [5, 15]},
            index=pd.MultiIndex.from_tuples(
                [
                    ("familyA", "genusA", "speciesA"),
                    ("familyB", "genusB", "speciesB"),
                ],
                names=["family", "genus", "species"],
            ),
        )

        result = self.binding.get_priming_ratio(
            max_mismatch_full_len, max_mismatch_three_end
        )

        expected = pd.DataFrame(
            {"priming_ratio": [0.50, 0.75]},
            index=pd.MultiIndex.from_tuples(
                [
                    ("familyA", "genusA", "speciesA"),
                    ("familyB", "genusB", "speciesB"),
                ],
                names=["family", "genus", "species"],
            ),
        )

        self.assertTrue("priming_ratio" in result.columns)
        pd.testing.assert_frame_equal(result, expected)

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
        self.test_dir = Path(__file__).resolve().parent / "data/test_data"
        self.amplicon_dir = self.test_dir / "amplicon-test"
        self.insert_dir = self.test_dir / "insert-test"
        self.incomplete_pbs_dir = self.test_dir / "insert-test"
        self.traits = TraitsAndResolution(
            insert_folder_path=self.insert_dir,
            amplicon_folder_path=self.amplicon_dir,
            incomplete_pbs_folder_path=self.incomplete_pbs_dir,
            otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv",
        )
        self.created_files = []
        self.traits.multibarcode_output_folder = os.path.join(
            os.path.dirname(self.traits.insert_folder_path), "multibarcode"
        )
        self.primer_resolv_species = pd.DataFrame(
            {
                "primer": ["primerA", "primerB", "primerC"],
                "additional_resolved_species": [90, 8, 2],
                "cumulative_resolved_species": [90, 98, 100],
            }
        )
        self.traits.primer_resolv_species = self.primer_resolv_species

    def test_init_with_results_folder(self):
        results_folder = "path/to/results"
        expected_insert_path = os.path.join(results_folder, "insert/filtered")
        expected_amplicon_path = os.path.join(results_folder, "amplicon")
        expected_incomplete_pbs_path = os.path.join(
            results_folder, "incomplete_pbs/filtered"
        )

        traits = TraitsAndResolution(
            results_folder=results_folder, otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv"
        )

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
            incomplete_pbs_folder_path=incomplete_pbs_path,
            otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv",
        )

        self.assertEqual(traits.insert_folder_path, insert_path)
        self.assertEqual(traits.amplicon_folder_path, amplicon_path)
        self.assertEqual(traits.incomplete_pbs_path, incomplete_pbs_path)

    def test_init_with_missing_arguments(self):
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution(otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv")
        self.assertIn(
            "Either provide a path to the in-silico amplification results folder",
            str(context.exception),
        )

    def test_init_with_partial_arguments(self):
        # Test missing amplicon_folder_path
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution(
                insert_folder_path="path/to/insert",
                incomplete_pbs_folder_path="path/to/incomplete_pbs",
                otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv",
            )
        self.assertIn(
            "Either provide a path to the in-silico amplification results folder",
            str(context.exception),
        )

        # Test missing insert_folder_path
        with self.assertRaises(ValueError) as context:
            TraitsAndResolution(
                amplicon_folder_path="path/to/amplicon",
                incomplete_pbs_folder_path="path/to/incomplete_pbs",
                otl=Path(__file__).resolve().parent / "data/test_data/test_otl.tsv",
            )
        self.assertIn(
            "Either provide a path to the in-silico amplification results folder",
            str(context.exception),
        )

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
        self.mock_otl_handler.otl_taxa_set = set(["species1", "species2"])

        # Set up path patches
        self.mock_paths = {
            "insert": "/fake/path/results/insert/filtered",
            "amplicon": "/fake/path/results/amplicon/filtered",
            "incomplete_pbs": "/fake/path/results/incomplete_pbs/filtered/input_AC",
        }

    @patch("os.path.join")
    @patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler")
    def test_initialization(self, mock_otl_handler_class, mock_join):
        """
        Test the initialization of MetricsSystemExecutor.
        """
        mock_otl_handler_instance = Mock()
        mock_otl_handler_class.return_value = mock_otl_handler_instance
        mock_join.side_effect = lambda *args: "/".join(args)

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table,
        )

        self.assertEqual(executor.results_folder, self.results_folder)
        self.assertEqual(executor.otl, self.otl)
        self.assertEqual(executor.primer_table, self.primer_table)
        mock_otl_handler_instance.import_otl.assert_called_once()

    @patch("os.makedirs")
    @patch("tempfile.TemporaryDirectory")
    @patch("os.listdir")
    @patch("src.mozaiko.marker_scoring.metrics_system.ReferenceDatabaseQuality")
    @patch("os.path.exists")
    @patch("os.path.join")
    @patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler")
    def test_get_reference_database_quality(
        self,
        mock_otl_handler_class,
        mock_join,
        mock_exists,
        mock_ref_db_class,
        mock_listdir,
        mock_temp_dir,
        mock_makedirs,
    ):
        """
        Test reference database quality calculation.
        """
        # Mock the temporary directory to return a known path
        mock_temp_dir_instance = Mock()
        mock_temp_dir_instance.__enter__ = Mock(return_value="/fake/temp/dir")
        mock_temp_dir_instance.__exit__ = Mock(return_value=None)
        mock_temp_dir.return_value = mock_temp_dir_instance

        # Mock os.makedirs to avoid actual file system operations
        mock_makedirs.return_value = None

        mock_exists.return_value = True
        mock_join.side_effect = lambda *args: "/".join(args)
        # Mock os.listdir to return an empty list to avoid file operations
        mock_listdir.return_value = []

        # Mock OtlHandler
        mock_otl_instance = Mock()
        mock_otl_instance.total_taxa = 100
        mock_otl_instance.otl_taxa_set = set(["species1", "species2"])
        mock_otl_handler_class.return_value = mock_otl_instance

        # set ReferenceDatabaseQuality mock
        mock_ref_db_instance = Mock()
        mock_ref_db_class.return_value = mock_ref_db_instance
        expected_result = pd.DataFrame(
            {"barcoded_taxa": [80.0], "ratio_barcoded_taxa": [0.8]}
        )
        mock_ref_db_instance.barcoded_taxa_ratio.return_value = expected_result

        with patch("builtins.open", create=True):
            executor = MetricsSystemExecutor(
                results_folder=self.results_folder,
                otl=self.otl,
                primer_table=self.primer_table,
            )

            result = executor.get_reference_database_quality()

            # Verify results
            pd.testing.assert_frame_equal(result, expected_result)
            expected_calls = [
                call(
                    otl="/fake/path/otl.tsv",
                    all_inserts_path="/fake/path/results/all_complete_pbs/filtered",
                ),
                call().barcoded_taxa_ratio(total_taxa_count=100),
            ]
            mock_ref_db_class.assert_has_calls(expected_calls)

            # Verify ReferenceDatabaseQuality was initialized correctly
            mock_ref_db_instance.barcoded_taxa_ratio.assert_called_once_with(
                total_taxa_count=executor.total_otl_taxa_count
            )

    @patch("src.mozaiko.marker_scoring.metrics_system.Binding")
    @patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler")
    @patch("os.path.exists")
    @patch("os.path.join")
    def test_get_primer_pbs_analysis(
        self, mock_join, mock_exists, mock_otl_handler_class, mock_binding_class
    ):
        """
        Test primer PBS analysis.
        """
        # Configure mock
        mock_binding_instance = Mock()
        mock_binding_class.return_value = mock_binding_instance

        expected_primer_pbs = {
            "primer1": {"data": "value1"},
            "primer2": {"data": "value2"},
        }
        expected_gc_df = pd.DataFrame(
            {"gc_content": [0.5, 0.6]}, index=["primer1", "primer2"]
        )

        mock_binding_instance.primer_pbs_analysis.return_value = (
            expected_primer_pbs,
            expected_gc_df,
        )

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table,
        )

        primer_pbs_dict, gc_df = executor.get_primer_pbs_analysis()

        self.assertEqual(primer_pbs_dict, expected_primer_pbs)
        pd.testing.assert_frame_equal(gc_df, expected_gc_df)
        mock_binding_instance.primer_pbs_analysis.assert_called_once_with(
            insert_folder=executor.insert_folder_path,
            amplicon_folder=executor.amplicon_folder_path,
            primer_table=executor.primer_table,
        )

    @patch("src.mozaiko.marker_scoring.metrics_system.TraitsAndResolution")
    @patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler")
    @patch("os.path.exists")
    @patch("os.path.join")
    def test_get_traits_and_resolution(
        self, mock_join, mock_exists, mock_otl_handler_class, mock_traits_class
    ):
        """
        Test traits and resolution calculation.
        """
        mock_exists.return_value = True
        mock_join.side_effect = lambda *args: "/".join(args)

        # Configure mock
        mock_otl_handler_instance = Mock()
        mock_otl_handler_instance.get_taxa_count.return_value = 100
        mock_otl_handler_class.return_value = mock_otl_handler_instance

        mock_traits_instance = Mock()
        mock_traits_class.return_value = mock_traits_instance

        # Mock expected DataFrame with 'primer' index
        mock_traits_instance.get_taxonomic_resolution.return_value = pd.DataFrame(
            {"primer": ["primer1"], "taxonomic_resolution": [0.75],
             "ratio_taxonomic_resolution": [0.75]}
        ).set_index("primer")

        # Create instance with mocked dependencies
        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table,
        )
        executor.total_otl_taxa_count = 100

        # Run test
        result = executor.get_traits_and_resolution()

        # Verify results
        pd.testing.assert_frame_equal(
            result, mock_traits_instance.get_taxonomic_resolution.return_value
        )

    @patch("src.mozaiko.marker_scoring.metrics_system.OtlHandler")
    @patch("os.path.exists")
    @patch("os.path.join")
    @patch("os.makedirs")
    @patch("builtins.open", new_callable=mock_open)  # Mock the file open operation
    @patch("pandas.DataFrame.to_csv")  # Patch the DataFrame.to_csv method directly
    def test_rank_primers(
        self,
        mock_to_csv,
        mock_file_open,
        mock_makedirs,
        mock_join,
        mock_exists,
        mock_otl_handler_class,
    ):
        """
        Test primer ranking functionality.
        """
        mock_analysis_results = pd.DataFrame(
            {
                "primer": ["primer1", "primer2"],
                "barcoded_taxa": [90, 80],
                "ratio_barcoded_taxa": [0.9, 0.8],
                "normalized_mismatch_score": [2, 3],
                "normalized_priming_ratio_sum": [0.8, 0.7],
                "normalized_gc_matches_across_taxon": [15, 12],
                "min_tm_cv": [0.1, 0.2],
                "tm_score": [0.9, 0.8],
                "amplification_success_percent": [95, 85],
                "taxonomic_resolution": [0.2, 0.3],
                "ratio_taxonomic_resolution": [0.2, 0.3],
            },
            index=["primer1", "primer2"],
        )

        executor = MetricsSystemExecutor(
            results_folder=self.results_folder,
            otl=self.otl,
            primer_table=self.primer_table,
        )
        executor.join_analysis_results = Mock(return_value=mock_analysis_results)

        mock_join.return_value = "/mock/directory/output_file.tsv"
        mock_exists.return_value = True
        mock_makedirs.return_value = None
        mock_to_csv.return_value = None

        result = executor.rank_primers_flat(save_intermediate_ranks=False)
        result = result.set_index("primer")

        self.assertTrue("final_rank" in result.columns)
        self.assertEqual(len(result), 2)
        primer1_rank = result.loc["primer1", "final_rank"]
        primer2_rank = result.loc["primer2", "final_rank"]
        self.assertLess(primer1_rank, primer2_rank)

def _make_traits() -> TraitsAndResolution:
    """Return a TraitsAndResolution instance using the standard test OTL."""
    test_data = Path(__file__).resolve().parent / "data/test_data"
    return TraitsAndResolution(
        insert_folder_path=str(test_data / "insert-test"),
        amplicon_folder_path=str(test_data / "amplicon-test"),
        incomplete_pbs_folder_path=str(test_data / "insert-test"),
        otl=str(test_data / "test_otl.tsv"),
    )

class TestCleanAndSplit(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    # --- query prefix ---

    def test_query_splits_clean_tuple_string(self):
        """Standard tuple string '(familyA, genusA, speciesA)' is parsed correctly."""
        df = pd.DataFrame({"query_cat": ["('familyA', 'genusA', 'speciesA')"]})
        result = self.traits.clean_and_split(df, "query")
        self.assertEqual(result["query_family"].iloc[0], "familyA")
        self.assertEqual(result["query_genus"].iloc[0], "genusA")
        self.assertEqual(result["query_species"].iloc[0], "speciesA")

    def test_query_strips_whitespace(self):
        """Leading/trailing whitespace in each field is removed."""
        df = pd.DataFrame({"query_cat": ["(' familyA ', ' genusA ', ' speciesA ')"]})
        result = self.traits.clean_and_split(df, "query")
        self.assertEqual(result["query_family"].iloc[0], "familyA")
        self.assertEqual(result["query_genus"].iloc[0], "genusA")
        self.assertEqual(result["query_species"].iloc[0], "speciesA")

    def test_query_creates_expected_columns(self):
        """The three new columns are added and the original _cat column is preserved."""
        df = pd.DataFrame({"query_cat": ["('familyA', 'genusA', 'speciesA')"]})
        result = self.traits.clean_and_split(df, "query")
        for col in ["query_family", "query_genus", "query_species"]:
            self.assertIn(col, result.columns)
        self.assertIn("query_cat", result.columns)

    def test_query_multiple_rows(self):
        """Multiple rows are all parsed correctly."""
        df = pd.DataFrame({
            "query_cat": [
                "('familyA', 'genusA', 'speciesA')",
                "('familyB', 'genusB', 'speciesB')",
            ]
        })
        result = self.traits.clean_and_split(df, "query")
        self.assertEqual(result["query_family"].tolist(), ["familyA", "familyB"])
        self.assertEqual(result["query_genus"].tolist(), ["genusA", "genusB"])
        self.assertEqual(result["query_species"].tolist(), ["speciesA", "speciesB"])

    def test_query_nan_value_is_kept(self):
        """A NaN entry in the source column produces NaN split columns (not an error)."""
        df = pd.DataFrame({"query_cat": [np.nan]})
        result = self.traits.clean_and_split(df, "query")
        # All three split columns must exist; values may be NaN
        for col in ["query_family", "query_genus", "query_species"]:
            self.assertIn(col, result.columns)

    # --- target prefix ---

    def test_target_splits_correctly(self):
        """Same parsing logic works with the 'target' prefix."""
        df = pd.DataFrame({"target_cat": ["('familyX', 'genusX', 'speciesX')"]})
        result = self.traits.clean_and_split(df, "target")
        self.assertEqual(result["target_family"].iloc[0], "familyX")
        self.assertEqual(result["target_genus"].iloc[0], "genusX")
        self.assertEqual(result["target_species"].iloc[0], "speciesX")

    def test_target_creates_expected_columns(self):
        df = pd.DataFrame({"target_cat": ["('familyX', 'genusX', 'speciesX')"]})
        result = self.traits.clean_and_split(df, "target")
        for col in ["target_family", "target_genus", "target_species"]:
            self.assertIn(col, result.columns)

    def test_removes_brackets_and_quotes(self):
        """Parentheses and single-quotes are stripped before splitting."""
        df = pd.DataFrame({"query_cat": ["('familyA', 'genusA', 'speciesA')"]})
        result = self.traits.clean_and_split(df, "query")
        # No residual punctuation
        for col in ["query_family", "query_genus", "query_species"]:
            self.assertNotIn("(", result[col].iloc[0])
            self.assertNotIn(")", result[col].iloc[0])
            self.assertNotIn("'", result[col].iloc[0])

class TestMakeDataframeSymmetric(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()
        # Minimal valid dataframe with all required columns
        self.base_df = pd.DataFrame({
            "query":          ["seqA"],
            "query_family":   ["familyA"],
            "query_genus":    ["genusA"],
            "query_species":  ["speciesA"],
            "target":         ["seqB"],
            "target_family":  ["familyB"],
            "target_genus":   ["genusB"],
            "target_species": ["speciesB"],
            "divergence_prct": [5.0],
        })

    def test_output_has_double_the_rows(self):
        """Symmetric df must contain exactly 2× the original rows."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        self.assertEqual(len(result), 2 * len(self.base_df))

    def test_original_rows_preserved(self):
        """The original rows appear unchanged in the output."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        original_row = result.iloc[0]
        self.assertEqual(original_row["query"],         "seqA")
        self.assertEqual(original_row["query_family"],  "familyA")
        self.assertEqual(original_row["target"],        "seqB")
        self.assertEqual(original_row["target_family"], "familyB")

    def test_swapped_rows_present(self):
        """A swapped row (query ↔ target) must also appear in the output."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        swapped = result[result["query"] == "seqB"]
        self.assertFalse(swapped.empty, "Swapped row not found in symmetric df")
        self.assertEqual(swapped.iloc[0]["target"],        "seqA")
        self.assertEqual(swapped.iloc[0]["query_family"],  "familyB")
        self.assertEqual(swapped.iloc[0]["target_family"], "familyA")

    def test_divergence_prct_preserved_in_both_rows(self):
        """divergence_prct must be the same in both the original and swapped row."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        self.assertTrue((result["divergence_prct"] == 5.0).all())

    def test_column_order_starts_with_query_target_block(self):
        """Output columns must begin with the canonical query/target block."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        expected_prefix = [
            "query", "query_family", "query_genus", "query_species",
            "target", "target_family", "target_genus", "target_species",
        ]
        self.assertEqual(list(result.columns[:8]), expected_prefix)

    def test_multiple_rows_symmetric(self):
        """All rows are mirrored when the input has multiple rows."""
        df = pd.DataFrame({
            "query":          ["seqA", "seqC"],
            "query_family":   ["familyA", "familyC"],
            "query_genus":    ["genusA", "genusC"],
            "query_species":  ["speciesA", "speciesC"],
            "target":         ["seqB", "seqD"],
            "target_family":  ["familyB", "familyD"],
            "target_genus":   ["genusB", "genusD"],
            "target_species": ["speciesB", "speciesD"],
            "divergence_prct": [3.0, 7.0],
        })
        result = self.traits.make_dataframe_symmetric(df)
        self.assertEqual(len(result), 4)
        # Both original queries appear as targets somewhere
        self.assertIn("seqA", result["target"].values)
        self.assertIn("seqC", result["target"].values)

    def test_extra_columns_are_retained(self):
        """Extra columns beyond the standard block are retained in the output."""
        df = self.base_df.copy()
        df["extra_col"] = "extra_value"
        result = self.traits.make_dataframe_symmetric(df)
        self.assertIn("extra_col", result.columns)

    def test_index_is_reset(self):
        """The result index must be a default RangeIndex (no duplicates from concat)."""
        result = self.traits.make_dataframe_symmetric(self.base_df)
        self.assertEqual(list(result.index), list(range(len(result))))

    def test_original_dataframe_not_mutated(self):
        """The input dataframe must not be modified in place."""
        original_copy = self.base_df.copy()
        self.traits.make_dataframe_symmetric(self.base_df)
        pd.testing.assert_frame_equal(self.base_df, original_copy)

class TestAddTaxaFromMapping(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()
        # A minimal existing dataframe that already contains one taxon pair
        self.existing_df = pd.DataFrame({
            "query":          ["seq1"],
            "query_family":   ["familyA"],
            "query_genus":    ["genusA"],
            "query_species":  ["speciesA"],
            "target":         ["seq2"],
            "target_family":  ["familyB"],
            "target_genus":   ["genusB"],
            "target_species": ["speciesB"],
            "divergence_prct": [3.0],
        })

    def _write_mapping(self, tmp_path: Path, rows: list[dict]) -> Path:
        """Write a mapping TSV to tmp_path and return its path."""
        header = [
            "seq_id", "og_taxa", "scientificName", "rank",
            "kingdom", "phylum", "class", "order",
            "family", "genus", "species",
        ]
        df = pd.DataFrame(rows, columns=header)
        out = tmp_path / "mapping_primer1.tsv"
        df.to_csv(out, sep="\t", index=False)
        return tmp_path

    def test_no_mapping_file_returns_df_unchanged(self):
        """When no mapping file exists the original df is returned as-is."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        pd.testing.assert_frame_equal(
            result.reset_index(drop=True),
            self.existing_df.reset_index(drop=True),
        )

    def test_no_new_taxa_when_all_already_present(self):
        """Row count is unchanged when every mapping taxon is already in the df."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq1", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyA", "genus": "genusA", "species": "speciesA"},
                {"seq_id": "seq2", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyB", "genus": "genusB", "species": "speciesB"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        self.assertEqual(len(result), len(self.existing_df))

    def test_new_taxa_are_appended(self):
        """A taxon present in the mapping but absent from the df is appended."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq1", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyA", "genus": "genusA", "species": "speciesA"},
                # NEW taxon not in existing_df
                {"seq_id": "seq3", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyC", "genus": "genusC", "species": "speciesC"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        self.assertEqual(len(result), len(self.existing_df) + 1)

    def test_new_row_has_inf_divergence(self):
        """The divergence_prct for a newly added taxon must be np.inf."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq3", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyC", "genus": "genusC", "species": "speciesC"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        new_row = result[result["query_family"] == "familyC"].iloc[0]
        self.assertTrue(np.isinf(new_row["divergence_prct"]))

    def test_new_row_target_columns_are_nan_string(self):
        """Target columns for a newly added taxon are set to the string 'nan'."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq3", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyC", "genus": "genusC", "species": "speciesC"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        new_row = result[result["query_family"] == "familyC"].iloc[0]
        self.assertEqual(new_row["target_family"],  "nan")
        self.assertEqual(new_row["target_genus"],   "nan")
        self.assertEqual(new_row["target_species"], "nan")
        self.assertEqual(new_row["target"],         "nan")

    def test_new_row_query_id_matches_mapping(self):
        """The query seq_id for a new row is looked up from the mapping file."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq99", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyC", "genus": "genusC", "species": "speciesC"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        new_row = result[result["query_family"] == "familyC"].iloc[0]
        self.assertEqual(new_row["query"], "seq99")

    def test_multiple_new_taxa_all_appended(self):
        """Multiple new taxa are all appended, each with inf divergence."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seqX", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyX", "genus": "genusX", "species": "speciesX"},
                {"seq_id": "seqY", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyY", "genus": "genusY", "species": "speciesY"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        new_rows = result[result["query_family"].isin(["familyX", "familyY"])]
        self.assertEqual(len(new_rows), 2)
        self.assertTrue(new_rows["divergence_prct"].apply(np.isinf).all())

    def test_target_taxa_also_checked_for_duplicates(self):
        """A taxon that appears only as a target (not query) in the df is not re-added."""
        # familyB/genusB/speciesB is a target in existing_df — it must not be duplicated
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq2", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyB", "genus": "genusB", "species": "speciesB"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        self.assertEqual(len(result), len(self.existing_df))

    def test_warning_printed_when_no_mapping_file(self):
        """A warning message is printed when no mapping file is found."""
        with tempfile.TemporaryDirectory() as tmp, \
             patch("builtins.print") as mock_print:
            folder = Path(tmp)
            self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer_missing"
            )
        printed = " ".join(str(c) for c in mock_print.call_args_list)
        self.assertIn("WARNING", printed)

    def test_original_rows_retained_after_append(self):
        """Existing rows are not modified when new taxa are appended."""
        with tempfile.TemporaryDirectory() as tmp:
            folder = self._write_mapping(Path(tmp), [
                {"seq_id": "seq1", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyA", "genus": "genusA", "species": "speciesA"},
                {"seq_id": "seq3", "og_taxa": None, "scientificName": None,
                 "rank": None, "kingdom": None, "phylum": None,
                 "class": None, "order": None,
                 "family": "familyC", "genus": "genusC", "species": "speciesC"},
            ])
            result = self.traits.add_taxa_from_mapping(
                self.existing_df.copy(), folder, "primer1"
            )
        original_rows = result.iloc[:len(self.existing_df)]
        pd.testing.assert_frame_equal(
            original_rows.reset_index(drop=True),
            self.existing_df.reset_index(drop=True),
        )


def _make_traits() -> TraitsAndResolution:
    """Return a TraitsAndResolution instance using the standard test OTL."""
    test_data = Path(__file__).resolve().parent / "data/test_data"
    return TraitsAndResolution(
        insert_folder_path=str(test_data / "insert-test"),
        amplicon_folder_path=str(test_data / "amplicon-test"),
        incomplete_pbs_folder_path=str(test_data / "insert-test"),
        otl=str(test_data / "test_otl.tsv"),
    )


class TestRunCatnip(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    @patch("subprocess.run")
    @patch("shutil.copy2")
    @patch("os.makedirs")
    @patch("os.listdir")
    def test_run_catnip_processes_fasta_files(
        self,
        mock_listdir,
        mock_makedirs,
        mock_copy2,
        mock_subprocess,
    ):
        """FASTA files are processed and subprocess is called."""
        mock_listdir.return_value = ["primer1.fasta", "ignore.txt"]
        mock_subprocess.return_value = MagicMock(returncode=0)

        self.traits.run_catnip(10.0)

        self.assertEqual(mock_subprocess.call_count, 1)
        args = mock_subprocess.call_args[0][0]
        self.assertIn("primer1.fasta", args)

    def test_run_catnip_invalid_threshold_type_raises(self):
        """Invalid threshold type raises TypeError."""
        with self.assertRaises(TypeError):
            self.traits.run_catnip("invalid")

    def test_run_catnip_invalid_threshold_list_length_raises(self):
        """Threshold list must contain exactly 3 values."""
        with self.assertRaises(TypeError):
            self.traits.run_catnip([1.0, 2.0])


class TestProcessSinglePrimer(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    @patch.object(TraitsAndResolution, "add_taxa_from_mapping")
    @patch.object(TraitsAndResolution, "get_values_otl")
    @patch.object(TraitsAndResolution, "filter_results_by_otl")
    @patch.object(TraitsAndResolution, "filter_divergence_threshold")
    def test_process_single_primer_runs_pipeline(
        self,
        mock_filter_div,
        mock_filter_otl,
        mock_get_values,
        mock_add_mapping,
    ):
        """Pipeline methods are called and processed dataframes returned."""
        df = pd.DataFrame({
            "query": ["q1"],
            "query_cat": ["('famA', 'genA', 'spA')"],
            "target": ["t1"],
            "target_cat": ["('famB', 'genB', 'spB')"],
            "divergence_prct": [2.0],
        })

        processed = pd.DataFrame({
            "query": ["q1"],
            "query_family": ["famA"],
            "query_genus": ["genA"],
            "query_species": ["spA"],
            "target": ["t1"],
            "target_family": ["famB"],
            "target_genus": ["genB"],
            "target_species": ["spB"],
            "divergence_prct": [2.0],
        })

        mock_add_mapping.return_value = processed
        mock_get_values.return_value = processed
        mock_filter_otl.return_value = processed
        mock_filter_div.return_value = processed

        with tempfile.TemporaryDirectory() as tmpdir:
            result1, result2 = self.traits.process_single_primer(
                "primer1",
                df,
                Path(tmpdir),
                thresholds=[10.0, 5.0, 2.0],
            )

        self.assertTrue(mock_add_mapping.called)
        self.assertTrue(mock_get_values.called)
        self.assertTrue(mock_filter_otl.called)
        self.assertTrue(mock_filter_div.called)

        pd.testing.assert_frame_equal(result1, processed)
        pd.testing.assert_frame_equal(result2, processed)

    def test_process_single_primer_invalid_thresholds_raise(self):
        """Invalid threshold lists raise TypeError."""
        df = pd.DataFrame()

        with self.assertRaises(TypeError):
            self.traits.process_single_primer(
                "primer",
                df,
                Path("."),
                thresholds=[1.0, 2.0],
            )


class TestPostProcessCatnipPrimerResults(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    @patch.object(TraitsAndResolution, "process_single_primer")
    def test_post_process_calls_process_single_primer(self, mock_process):
        """Each valid primer folder is processed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            self.traits.catnip_dir = tmpdir

            primer_dir = Path(tmpdir) / "primer1"
            primer_dir.mkdir()

            df = pd.DataFrame({
                "query": ["q1"],
                "query_cat": ["('famA', 'genA', 'spA')"],
                "target": ["t1"],
                "target_cat": ["('famB', 'genB', 'spB')"],
                "divergence_prct": [2.0],
            })

            df.to_csv(
                primer_dir / "final_output_interclst.tsv",
                sep="\t",
                index=False,
            )

            self.traits.post_process_catnip_primer_results()

            self.assertEqual(mock_process.call_count, 1)


class TestGetValuesOtl(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    def test_get_values_otl_returns_minimum_divergence(self):
        """Lowest divergence value is selected."""
        self.traits.otl = pd.DataFrame({
            "family": ["famA"],
            "genus": ["genA"],
            "species": ["spA"],
            "rank": ["species"],
        })

        self.traits.df_inf = pd.DataFrame()

        df = pd.DataFrame({
            "query_family": ["famA", "famA"],
            "query_genus": ["genA", "genA"],
            "query_species": ["spA", "spA"],
            "target_family": ["famB", "famC"],
            "target_genus": ["genB", "genC"],
            "target_species": ["spB", "spC"],
            "divergence_prct": [5.0, 2.0],
        })

        result = self.traits.get_values_otl(df)

        self.assertEqual(result["divergence_prct"].iloc[0], 2.0)
        self.assertEqual(result["target_species"].iloc[0], "spC")

    def test_get_values_otl_multiple_equal_mins_returns_nan_targets(self):
        """Equal minimum values result in NaN targets."""
        self.traits.otl = pd.DataFrame({
            "family": ["famA"],
            "genus": ["genA"],
            "species": ["spA"],
            "rank": ["species"],
        })

        self.traits.df_inf = pd.DataFrame()

        df = pd.DataFrame({
            "query_family": ["famA", "famA"],
            "query_genus": ["genA", "genA"],
            "query_species": ["spA", "spA"],
            "target_family": ["famB", "famC"],
            "target_genus": ["genB", "genC"],
            "target_species": ["spB", "spC"],
            "divergence_prct": [2.0, 2.0],
        })

        result = self.traits.get_values_otl(df)

        self.assertTrue(pd.isna(result["target_species"].iloc[0]))
        self.assertEqual(result["divergence_prct"].iloc[0], 2.0)


class TestFilterResultsByOtl(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    def test_filter_results_by_otl_keeps_matching_targets(self):
        """Rows whose target taxa exist in OTL are retained."""
        self.traits.otl = pd.DataFrame({
            "family": ["famB"],
            "genus": ["genB"],
            "species": ["spB"],
            "rank": ["species"],
        })

        df = pd.DataFrame({
            "target_family": ["famB", "famX"],
            "target_genus": ["genB", "genX"],
            "target_species": ["spB", "spX"],
            "divergence_prct": [1.0, 2.0],
        })

        result = self.traits.filter_results_by_otl(df)

        self.assertEqual(len(result), 1)
        self.assertEqual(result["target_species"].iloc[0], "spB")

    def test_filter_results_by_otl_keeps_inf_rows(self):
        """Rows with infinite divergence are preserved."""
        self.traits.otl = pd.DataFrame({
            "family": ["famB"],
            "genus": ["genB"],
            "species": ["spB"],
            "rank": ["species"],
        })

        df = pd.DataFrame({
            "target_family": ["famX"],
            "target_genus": ["genX"],
            "target_species": ["spX"],
            "divergence_prct": [np.inf],
        })

        result = self.traits.filter_results_by_otl(df)

        self.assertEqual(len(result), 1)
        self.assertTrue(np.isinf(result["divergence_prct"].iloc[0]))


class TestExcludeCommonAncestry(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    def test_excludes_same_species(self):
        """Rows with identical query/target species are removed."""
        df = pd.DataFrame({
            "query_family": ["famA"],
            "query_genus": ["genA"],
            "query_species": ["spA"],
            "target_family": ["famA"],
            "target_genus": ["genA"],
            "target_species": ["spA"],
            "divergence_prct": [1.0],
        })

        result = self.traits.exclude_common_ancestry(df)

        self.assertEqual(len(result), 0)

    def test_keeps_different_species(self):
        """Rows with different species are retained."""
        df = pd.DataFrame({
            "query_family": ["famA"],
            "query_genus": ["genA"],
            "query_species": ["spA"],
            "target_family": ["famA"],
            "target_genus": ["genA"],
            "target_species": ["spB"],
            "divergence_prct": [1.0],
        })

        result = self.traits.exclude_common_ancestry(df)

        self.assertEqual(len(result), 1)


class TestFilterDivergenceThreshold(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    def test_filter_divergence_threshold_single_float(self):
        """Rows above single threshold are retained."""
        df = pd.DataFrame({
            "query_family": ["famA"],
            "query_genus": [np.nan],
            "query_species": [np.nan],
            "target_family": ["famB"],
            "target_genus": [np.nan],
            "target_species": [np.nan],
            "divergence_prct": [15.0],
        })

        result = self.traits.filter_divergence_threshold(df, 10.0)

        self.assertEqual(len(result), 1)

    def test_filter_divergence_threshold_list(self):
        """Rank-specific thresholds are applied."""
        df = pd.DataFrame({
            "query_family": ["famA"],
            "query_genus": ["genA"],
            "query_species": [np.nan],
            "target_family": ["famB"],
            "target_genus": ["genB"],
            "target_species": [np.nan],
            "divergence_prct": [6.0],
        })

        result = self.traits.filter_divergence_threshold(
            df,
            [10.0, 5.0, 2.0],
        )

        self.assertEqual(len(result), 1)

    def test_filter_divergence_threshold_list(self):
        """Rank-specific thresholds are applied."""
        df = pd.DataFrame({
            "query_family": ["famA"],
            "query_genus": ["genA"],
            "query_species": [np.nan],
            "target_family": ["famB"],
            "target_genus": ["genB"],
            "target_species": [np.nan],
            "divergence_prct": [4.0],
        })

        result = self.traits.filter_divergence_threshold(
            df,
            [10.0, 5.0, 2.0],
        )

        self.assertEqual(len(result), 0)

    def test_filter_divergence_threshold_invalid_type_raises(self):
        """Invalid threshold types raise TypeError."""
        df = pd.DataFrame()

        with self.assertRaises(TypeError):
            self.traits.filter_divergence_threshold(df, "bad")


class TestGetTaxonomicResolution(unittest.TestCase):

    def setUp(self):
        self.traits = _make_traits()

    def test_get_taxonomic_resolution_returns_dataframe(self):
        """Taxonomic resolution metrics are calculated."""
        self.traits.otl_handler.total_taxa = 10

        with tempfile.TemporaryDirectory() as tmpdir:
            self.traits.catnip_dir = tmpdir
            self.traits.country_name = "test"

            primer_dir = Path(tmpdir) / "primer1"
            primer_dir.mkdir()

            catnip_df = pd.DataFrame({
                "divergence_prct": [1.0, 2.0],
            })

            otl_df = pd.DataFrame({
                "divergence_prct": [1.0],
            })

            catnip_df.to_csv(
                primer_dir / "catnip_target_primer1_test.tsv",
                sep="\t",
                index=False,
            )

            otl_df.to_csv(
                primer_dir / "otl_target_primer1_test.tsv",
                sep="\t",
                index=False,
            )

            result = self.traits.get_taxonomic_resolution()

            print(result)

            self.assertEqual(len(result), 1)
            self.assertEqual(result["primer"].iloc[0], "primer1")
            self.assertEqual(result["taxonomic_resolution"].iloc[0], 0.1)
            self.assertEqual(result["ratio_taxonomic_resolution"].iloc[0], 2.0)