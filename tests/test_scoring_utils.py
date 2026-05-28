"""
Unit tests for scoring_utils.py and aux_analysis.py.
"""

import shutil
import tempfile
import unittest
from io import StringIO
from pathlib import Path
import os
import shutil
from unittest.mock import MagicMock, patch

import pandas as pd

from src.mozaiko.marker_scoring.scoring_utils import *
from src.mozaiko.marker_scoring.aux_analysis import *

from src.mozaiko.marker_scoring.aux_analysis import (
    _calculate_tax_coverage,
    _process_fasta_file,
)


class TestScoringUtils(unittest.TestCase):
    def setUp(self):
        self.test_directory = Path(__file__).resolve().parent / "data/test_data"

    def test_calculate_iupac_mismatches(self):
        """
        Function to test the calculate_iupac_mismatches function.
        """
        mismatches_cases = self.test_directory / "iupac_mismatches_cases.txt"

        with open(mismatches_cases, "r") as file:
            for line in file:
                seq1, seq2, expected_mismatches, gc_matches = line.strip().split(" ")
                mismatches = calculate_iupac_mismatches(seq1, seq2)
                self.assertEqual(mismatches, int(expected_mismatches))

    def test_calculate_gc_matches(self):
        """
        Function to test the search_gc_clamp functonality, within the calculate_iupac_mismatches function.
        """
        mismatches_cases = self.test_directory / "iupac_mismatches_cases.txt"

        with open(mismatches_cases, "r") as file:
            for line in file:
                print(line)
                seq1, seq2, expected_mismatches, expected_gc_matches = (
                    line.strip().split(" ")
                )
                mismatches, real_gc_matches = calculate_iupac_mismatches(
                    seq1, seq2, search_gc_clamp=True
                )
                self.assertEqual(real_gc_matches, int(expected_gc_matches))

    def test_calculate_ambiguous_percentage(self):
        """
        Function to test the calculate_ambiguous_percentage.
        """
        sequences = ["ACGTACGTAG", "ACGTANGATC", "RWNBVDNWSN"]  # 0%  # 10%  # 100%
        ambiguity = []

        for sequence in sequences:
            result = calculate_ambiguous_percentage(sequence)
            ambiguity.append(result)

        self.assertEqual(ambiguity[0], 0.0)

        self.assertEqual(ambiguity[1], 0.1)

        self.assertEqual(ambiguity[2], 1.0)

    def test_filter_sequences_by_ambiguity_not_real_path(self):
        """
        Test if filter_sequences_by_ambiguity detects a non-existent path.
        """
        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity("not_real_path")
        self.assertIn("does not exist", str(context.exception))

    def test_filter_sequences_by_ambiguity_empty_dir(self):
        """
        Test if filter_sequences_by_ambiguity detects a empty folder.
        """
        empty_dir = self.test_directory / "empty"
        empty_dir.mkdir()

        with self.assertRaises(ValueError) as context:
            filter_sequences_by_ambiguity(empty_dir)
        self.assertIn("No FASTA files found in", str(context.exception))

        shutil.rmtree(empty_dir)

    def test_filter_sequences_by_ambiguity_file(self):
        """
        Test if filter_sequences_by_ambiguity correctly reads a input file.
        """
        test_file = self.test_directory / "ambiguous_sequences.fasta"

        results = filter_sequences_by_ambiguity(
            test_file, max_ambiguous_percentage=0.25
        )

        output_dir = self.test_directory / "filtered"
        self.assertTrue(output_dir.exists())

        output_file = self.test_directory / "filtered" / "ambiguous_sequences.fasta"
        self.assertTrue(output_file.exists())

        with open(output_file) as f:
            contents = f.read()

            self.assertIn(">seq1", contents)
            self.assertIn("ATGCATGC", contents)
            self.assertIn(">seq2", contents)
            self.assertIn("ATGNNTGC", contents)
            self.assertNotIn(">seq3", contents)
            self.assertNotIn("NNNNNNNN", contents)
            self.assertNotIn(">seq4", contents)
            self.assertNotIn(
                "ACGTACGACGAGCATTCGAAGGTCAGTCGNNNNATTAGCTACTGATCGATCGACTAGCTCCGCATCGATGATGCATGCTAGTCGATGCATGCATCG",
                contents,
            )

    def test_extract_primer_binding_sites(self):
        amplicon_file = self.test_directory / "amplicon-test/primerC.fasta"
        insert_file = self.test_directory / "insert-test/primerA.fasta"

        result = extract_primer_binding_sites(amplicon_file, insert_file)

        self.assertIsInstance(result, pd.DataFrame)

        self.assertEqual(
            list(result.columns),
            ["header", "fwd_seq", "rev_seq", "fwd_seq_len", "rev_seq_len"],
        )

        self.assertEqual(result["header"].iloc[0], ">abc")
        self.assertEqual(result["header"].iloc[1], ">def")
        self.assertEqual(result["fwd_seq"].iloc[0], "ACGTAGCA")
        self.assertEqual(result["rev_seq"].iloc[0], "ACCATCA")
        self.assertEqual(result["fwd_seq"].iloc[1], "AGTGA")
        self.assertEqual(result["rev_seq"].iloc[1], "GTCGATAGCAT")
        self.assertEqual(result["fwd_seq_len"].iloc[0], 8)
        self.assertEqual(result["rev_seq_len"].iloc[0], 7)
        self.assertEqual(result["fwd_seq_len"].iloc[1], 5)
        self.assertEqual(result["rev_seq_len"].iloc[1], 11)

    def tearDown(self):
        """
        Clean up any files created during tests.
        """
        filtered_dir = self.test_directory / "filtered"
        if filtered_dir.exists():
            shutil.rmtree(filtered_dir)


class TestSequenceCountTracking(unittest.TestCase):
    def setUp(self):
        self.analysis_folder = tempfile.TemporaryDirectory()
        self.original_database = tempfile.NamedTemporaryFile(
            delete=False, suffix=".fasta"
        )

        with open(self.original_database.name, "w") as f:
            f.write(">seq\nAGTGCA\n>seq2\nGTCAGCGA\n>seq3\nGGGGCA")

        required_dirs = [
            "all_complete_pbs/input_ABC",
            "insert/input_A",
            "input_B",
        ]

        for d in required_dirs:
            path = os.path.join(self.analysis_folder.name, d)
            os.makedirs(path, exist_ok=True)

            fasta = os.path.join(path, "primerA.fasta")

            with open(fasta, "w") as f:
                f.write(">seq1\nATGC\n>seq2\nATGC\n")

    def tearDown(self):
        self.analysis_folder.cleanup()
        if os.path.exists(self.original_database.name):
            os.remove(self.original_database.name)

    def test_sequence_count_tracking(self):
        df_pivoted = sequence_count_tracking(
            self.original_database.name, self.analysis_folder.name
        )

        self.assertIsInstance(df_pivoted, pd.DataFrame)

        run_name = os.path.basename(self.analysis_folder.name)

        output_path = os.path.join(
            self.analysis_folder.name, f"{run_name}-sequence_count_track.tsv"
        )

        self.assertTrue(os.path.exists(output_path))

        result_df = pd.read_csv(output_path, sep="\t")

        self.assertIn("primer_name", result_df.columns)

        # Check at least one analysis-step column exists (amplicon)
        self.assertTrue(len(result_df.columns) > 1)

        self.assertIn(
            "percentage_of_sequences_with_PBS",
            df_pivoted.columns
        )

        self.assertIn(
            "percentage_of_amplified_sequences_with_PBS",
            df_pivoted.columns
        )

        self.assertFalse(df_pivoted.empty)

        print(result_df)

class TestMultiBarcodeToolsInput(unittest.TestCase):
    def setUp(self):
        """
        Set up a temporary directory and create sample FASTA files for testing.
        """
        self.test_dir = tempfile.mkdtemp()
        self.create_sample_fasta_files()

    def tearDown(self):
        """
        Clean up the temporary directory after tests.
        """
        shutil.rmtree(self.test_dir)

    def create_sample_fasta_files(self):
        """
        Method to create sample FASTA files for testing different scenarios.
        """
        # FASTA with standard header
        with open(os.path.join(self.test_dir, "primer1.fasta"), "w") as f:
            f.write(">seq1|Species Name 1|harmonized_species|rank|kingdom|phylum|class|order|family|genus|species\n")
            f.write("ATCGATCGATCG\n")
            f.write(">seq2|Species Name 2|harmonized_species|rank|kingdom|phylum|class|order|family|genus|species\n")
            f.write("GCTAGCTAGCTA\n")

        # FASTA with multiple sequences
        with open(os.path.join(self.test_dir, "primer2.fasta"), "w") as f:
            f.write(">seq3|Species Name 3|harmonized_species|rank|kingdom|phylum|class|order|family|genus|species\n")
            f.write("TAGCTAGCTAGC\n")
            f.write(">seq4|Species Name 4|harmonized_species|rank|kingdom|phylum|class|order|family|genus|species\n")
            f.write("CGATCGATCGAT\n")

        # FASTA with problematic header
        with open(os.path.join(self.test_dir, "primer3.fasta"), "w") as f:
            f.write(">seq5\n")  # Missing species name
            f.write("ATCGATCG\n")

def _write_fasta(path, records):
    """Write a minimal FASTA file. records = [(header, seq), ...]"""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f"{header}\n{seq}\n")


def _write_tsv(path, rows):
    """Write a TSV with a taxa_info column. rows = [taxa_info_string, ...]"""
    df = pd.DataFrame({"taxa_info": rows})
    df.to_csv(path, sep="\t", index=False)


def _make_otl_handler(families, genera, species_list, total):
    """
    Return a mock OtlHandler whose .otl DataFrame has the expected shape.
    Mirrors the columns accessed in _calculate_tax_coverage:
        family, genus, species, rank
    """
    records = []
    for fam, gen, sp in zip(families, genera, species_list):
        records.append({"family": fam, "genus": gen, "species": sp, "rank": "species"})
    df = pd.DataFrame(records)

    handler = MagicMock()
    handler.total_taxa = total
    handler.otl = df
    return handler

class TestProcessFastaFile(unittest.TestCase):
    """Tests for _process_fasta_file."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # Directory structure mirrors what the function expects:
        #   input_file = <output_dir>/<step_name>/<file_name>.fasta
        self.step_dir = os.path.join(self.tmp, "step_name")
        os.makedirs(self.step_dir)
        self.output_dir = os.path.join(self.tmp, "output")
        os.makedirs(self.output_dir)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_returns_dataframe_when_tsv_cached(self):
        """Reads the pre-existing TSV and returns family/genus/species columns."""
        taxa_rows = [
            "seq1|Species A|harmonized|rank|k|p|c|o|Familyaceae|GenusA|SpeciesA",
            "seq2|Species B|harmonized|rank|k|p|c|o|Familyaceae|GenusB|SpeciesB",
        ]
        tsv_path = os.path.join(self.output_dir, "primerA_step_name.tsv")
        _write_tsv(tsv_path, taxa_rows)

        input_file = os.path.join(self.step_dir, "primerA.fasta")
        _write_fasta(input_file, [(">seq1", "ATGC")])

        result = _process_fasta_file(input_file, self.output_dir)

        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(list(result.columns), ["family", "genus", "species"])
        self.assertEqual(len(result), 2)
        self.assertEqual(result["family"].iloc[0], "Familyaceae")
        self.assertEqual(result["genus"].iloc[0], "GenusA")
        self.assertEqual(result["species"].iloc[0], "SpeciesA")

    def test_tsv_cached_last_three_pipe_fields_parsed(self):
        """Verifies that only the last three pipe-separated fields are used."""
        # Many extra fields before family|genus|species
        taxa_rows = [
            "seq1|Common Name|harmonized|rank|Animalia|Chordata|Mammalia|Carnivora|Felidae|Felis|catus",
        ]
        tsv_path = os.path.join(self.output_dir, "primerA_step_name.tsv")
        _write_tsv(tsv_path, taxa_rows)

        input_file = os.path.join(self.step_dir, "primerA.fasta")
        _write_fasta(input_file, [(">seq1", "ATGC")])

        result = _process_fasta_file(input_file, self.output_dir)

        self.assertEqual(result["family"].iloc[0], "Felidae")
        self.assertEqual(result["genus"].iloc[0], "Felis")
        self.assertEqual(result["species"].iloc[0], "catus")

    def test_creates_tsv_and_returns_dataframe_when_not_cached(self):
        """
        When no TSV is present the function calls CustomFastaImport and
        writes a TSV before returning the DataFrame.
        """
        taxa_rows = [
            "seq1|Species A|harmonized|rank|k|p|c|o|Familyaceae|GenusA|SpeciesA",
        ]

        input_file = os.path.join(self.step_dir, "primerA.fasta")
        _write_fasta(input_file, [(">seq1", "ATGC")])

        expected_tsv = os.path.join(self.output_dir, "primerA_step_name.tsv")

        mock_importer = MagicMock()

        def fake_df2csv(path):
            _write_tsv(path, taxa_rows)

        mock_importer.df2csv.side_effect = fake_df2csv

        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.CustomFastaImport",
            return_value=mock_importer,
        ):
            result = _process_fasta_file(input_file, self.output_dir)

        self.assertTrue(os.path.exists(expected_tsv))
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(list(result.columns), ["family", "genus", "species"])

    def test_output_dir_created_when_missing(self):
        """The output directory is created if it does not already exist."""
        new_output = os.path.join(self.tmp, "brand_new_output")
        self.assertFalse(os.path.exists(new_output))

        input_file = os.path.join(self.step_dir, "primerA.fasta")
        _write_fasta(input_file, [(">seq1", "ATGC")])

        taxa_rows = [
            "seq1|Species A|harmonized|rank|k|p|c|o|Familyaceae|GenusA|SpeciesA",
        ]
        mock_importer = MagicMock()
        mock_importer.df2csv.side_effect = lambda p: _write_tsv(p, taxa_rows)

        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.CustomFastaImport",
            return_value=mock_importer,
        ):
            _process_fasta_file(input_file, new_output)

        self.assertTrue(os.path.exists(new_output))


class TestCalculateTaxCoverage(unittest.TestCase):
    """Tests for _calculate_tax_coverage."""

    def _input_df(self, families, genera, species_list, rank="species"):
        return pd.DataFrame(
            {"family": families, "genus": genera, "species": species_list, "rank": rank}
        )

    def _mock_otl_path(self):
        return "/fake/otl.tsv"

    def test_full_coverage_returns_100(self):
        """All OTL taxa appear at least once: 100%."""
        input_df = self._input_df(
            ["Fam1", "Fam2"],
            ["Gen1", "Gen2"],
            ["sp1", "sp2"],
        )
        handler = _make_otl_handler(
            ["Fam1", "Fam2"], ["Gen1", "Gen2"], ["sp1", "sp2"], total=2
        )
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(input_df, self._mock_otl_path())

        self.assertAlmostEqual(result, 100.0)

    def test_partial_coverage(self):
        """Half the OTL taxa are covered: 50%."""
        input_df = self._input_df(["Fam1"], ["Gen1"], ["sp1"])
        handler = _make_otl_handler(
            ["Fam1", "Fam2"],
            ["Gen1", "Gen2"],
            ["sp1", "sp2"],
            total=2,
        )
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(input_df, self._mock_otl_path())

        self.assertAlmostEqual(result, 50.0)

    def test_zero_coverage(self):
        """No OTL taxa appear in the input: 0%."""
        input_df = self._input_df(["Unknown"], ["Unknown"], ["unknown_sp"])
        handler = _make_otl_handler(
            ["Fam1", "Fam2"],
            ["Gen1", "Gen2"],
            ["sp1", "sp2"],
            total=2,
        )
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(input_df, self._mock_otl_path())

        self.assertAlmostEqual(result, 0.0)

    def test_return_type_is_float(self):
        """Coverage is always returned as a plain Python float."""
        input_df = self._input_df(["Fam1"], ["Gen1"], ["sp1"])
        handler = _make_otl_handler(["Fam1"], ["Gen1"], ["sp1"], total=1)
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(input_df, self._mock_otl_path())

        self.assertIsInstance(result, float)

    def test_cutoff_val_filters_low_count_taxa(self):
        """
        cutff_val=2 means a taxon must appear at least twice.
        With the taxon appearing once it should not count.
        """
        input_df = self._input_df(["Fam1"], ["Gen1"], ["sp1"])
        handler = _make_otl_handler(["Fam1"], ["Gen1"], ["sp1"], total=1)
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(
                input_df, self._mock_otl_path(), cutff_val=2
            )

        self.assertAlmostEqual(result, 0.0)

    def test_cutoff_val_1_is_default_behaviour(self):
        """Default cutff_val=1: a single occurrence qualifies the taxon."""
        input_df = self._input_df(["Fam1"], ["Gen1"], ["sp1"])
        handler = _make_otl_handler(["Fam1"], ["Gen1"], ["sp1"], total=1)
        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ):
            result = _calculate_tax_coverage(input_df, self._mock_otl_path())

        self.assertAlmostEqual(result, 100.0)

    def test_otl_handler_called_with_correct_path(self):
        """OtlHandler receives the path passed to _calculate_tax_coverage."""
        otl_path = "/some/path/otl.tsv"
        input_df = self._input_df([], [], [], [])
        handler = _make_otl_handler(["Fam1"], ["Gen1"], ["sp1"], total=1)


        with patch(
            "src.mozaiko.marker_scoring.aux_analysis.OtlHandler",
            return_value=handler,
        ) as MockOtl:
            _calculate_tax_coverage(input_df, otl_path)

        MockOtl.assert_called_once_with(otl_path)
        handler.import_otl.assert_called_once()


class TestBarcodedTaxa(unittest.TestCase):
    """Integration-style tests for barcoded_taxa (mocks its two helpers)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.tmp, "output")
        os.makedirs(self.output_dir)
        # Minimal FASTA file — content doesn't matter here because
        # _process_fasta_file is mocked.
        self.fasta = os.path.join(self.tmp, "primerA.fasta")
        _write_fasta(self.fasta, [(">seq1", "ATGC")])
        self.otl = "/fake/otl.tsv"

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _sample_df(self):
        return pd.DataFrame(
            {
                "family": ["Familyaceae"],
                "genus": ["GenusA"],
                "species": ["SpeciesA"],
            }
        )

    def test_returns_coverage_float(self):
        """barcoded_taxa returns the float produced by _calculate_tax_coverage."""
        with (
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._process_fasta_file",
                return_value=self._sample_df(),
            ),
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._calculate_tax_coverage",
                return_value=75.0,
            ),
        ):
            result = barcoded_taxa(self.fasta, self.otl, self.output_dir)

        self.assertIsInstance(result, float)
        self.assertAlmostEqual(result, 75.0)

    def test_process_fasta_called_with_correct_args(self):
        """_process_fasta_file is called with input_ABC_file and output_dir."""
        with (
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._process_fasta_file",
                return_value=self._sample_df(),
            ) as mock_proc,
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._calculate_tax_coverage",
                return_value=0.0,
            ),
        ):
            barcoded_taxa(self.fasta, self.otl, self.output_dir)

        mock_proc.assert_called_once_with(self.fasta, self.output_dir)

    def test_calculate_tax_coverage_called_with_correct_args(self):
        """_calculate_tax_coverage receives the DataFrame and OTL path."""
        df = self._sample_df()
        with (
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._process_fasta_file",
                return_value=df,
            ),
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._calculate_tax_coverage",
                return_value=50.0,
            ) as mock_cov,
        ):
            barcoded_taxa(self.fasta, self.otl, self.output_dir)

        args, kwargs = mock_cov.call_args
        pd.testing.assert_frame_equal(args[0], df)
        self.assertEqual(args[1], self.otl)

    def test_zero_coverage_propagated(self):
        """A coverage of 0.0 is passed through without modification."""
        with (
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._process_fasta_file",
                return_value=self._sample_df(),
            ),
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._calculate_tax_coverage",
                return_value=0.0,
            ),
        ):
            result = barcoded_taxa(self.fasta, self.otl, self.output_dir)

        self.assertAlmostEqual(result, 0.0)

    def test_full_coverage_propagated(self):
        """A coverage of 100.0 is passed through without modification."""
        with (
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._process_fasta_file",
                return_value=self._sample_df(),
            ),
            patch(
                "src.mozaiko.marker_scoring.aux_analysis._calculate_tax_coverage",
                return_value=100.0,
            ),
        ):
            result = barcoded_taxa(self.fasta, self.otl, self.output_dir)

        self.assertAlmostEqual(result, 100.0)
