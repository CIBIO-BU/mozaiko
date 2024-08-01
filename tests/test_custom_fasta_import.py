"""
Unit tests for the CustomFastaImport class.
"""
import os
import unittest
from unittest.mock import patch, MagicMock
from src.sequence_import import CustomFastaImport, LineageFileLoader


class TestCustomFastaImport(unittest.TestCase):
    """
    Class to test the CustomFastaImport class.
    """

    def setUp(self):
        """
        Set up the test class and data.
        """
        self.data_dir = "data/test_data"
        self.fasta_file = os.path.join(self.data_dir, "fasta_example_file.fasta")
        self.fasta_taxid_file = os.path.join(self.data_dir, "fasta_example_file_taxid.fasta")
        self.fasta_import = CustomFastaImport(None)
        self.lineage_loader = MagicMock()
        self.lineage_loader.load_lineage_file.return_value = 'dummy_lineage.tsv'

    def test_empty_file(self):
        """
        Test if test_empty_file raises a ValueError when the input file is empty.
        """
        empty_file = os.path.join(self.data_dir, "empty_fasta.fasta")
        with self.assertRaises(ValueError) as context:
            self.fasta_import.read_fasta(empty_file)
        self.assertEqual(str(context.exception), 'Input file is empty.')

    def test_txt_file(self):
        """
        Test if test_txt_file raises a ValueError when the input file is not a fasta file.
        """
        txt_file = os.path.join(self.data_dir, "not_fasta_example.txt")
        with self.assertRaises(ValueError) as context:
            self.fasta_import.read_fasta(txt_file)
        self.assertEqual(str(context.exception), 'Input file must be a fasta file.')

    def test_missing_file(self):
        """
        Test if test_missing_file raises a FileNotFoundError when the input file does not exist.
        """
        missing_file = os.path.join(self.data_dir, "no_file.fasta")
        with self.assertRaises(FileNotFoundError) as context:
            self.fasta_import.read_fasta(missing_file)
        self.assertEqual(str(context.exception), 'Input file does not exist in the directory.')

    def test_wrong_input(self):
        """
        Test if test_wrong_input raises a ValueError when the input is not a string.
        """
        not_string = 1
        with self.assertRaises(ValueError) as context:
            self.fasta_import.read_fasta(not_string)
        self.assertEqual(str(context.exception), 'Directory must be a string.')

    def test_read_fasta(self):
        """
        Test if read_fasta reads the fasta file correctly.
        """
        data = self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(data.shape, (3, 4))

    def test_get_number_of_sequences(self):
        """
        Test if get_number_of_sequences returns the correct number of sequences.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.fasta_import.get_number_of_sequences(), 3)

    def test_get_sequence_lengths(self):
        """
        Test if get_sequence_lengths returns the correct lengths of the sequences.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.fasta_import.get_sequence_lengths(), [16, 19, 8])

    def test_get_sequence_ids(self):
        """
        Test if get_sequence_ids returns the correct sequence IDs.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.fasta_import.get_sequence_ids(),
                         ['CM074756.1', 'NC_088426.1', 'PP475397.1'])

    def test_get_sequences(self):
        """
        Test if get_sequences returns the correct sequences.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.fasta_import.get_sequences(),
                         ['GTTATTGTAGCTTATC', 'GCATAAAGCATGGCACTGA', 'GTTATTGA'])

    def test_df2fasta_structure(self):
        """
        Test if df2fasta correctly writes the data to a fasta file.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.fasta_import.df2csv()
        example_file = 'processed_input_fasta.csv'
        with open(example_file, 'r', encoding='UTF-8') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], 'seq_id,sequence,length,taxid\n')
            self.assertEqual(lines[1], 'CM074756.1,GTTATTGTAGCTTATC,16,8481\n')
            self.assertEqual(lines[2], 'NC_088426.1,GCATAAAGCATGGCACTGA,19,12345\n')
            self.assertEqual(lines[3], 'PP475397.1,GTTATTGA,8,106731\n')
        os.remove(example_file)

    def test_get_taxids(self):
        """
        Test if get_taxids returns the correct taxids.
        """
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.fasta_import.get_taxids(), ['8481', '12345', '106731'])

    def test_check_for_taxids_with_taxid(self):
        """
        Test if _check_for_taxids requests lineage file if taxids are found.
        """
        self.fasta_import.add_taxids = MagicMock()
        self.fasta_import.read_fasta(self.fasta_taxid_file)
        self.fasta_import.add_taxids.assert_called_once_with(self.fasta_taxid_file)

    @patch('builtins.print')
    @patch('Bio.SeqIO.parse')
    def test_check_for_taxids_without_taxid(self, mock_seqio_parse, mock_print):
        """
        Test if _check_for_taxids requests lineage file if taxids are not found.
        """
        # mock the SeqIO.parse function to return a list of sequences without taxids
        mock_seqio_parse.return_value = [
            MagicMock(description='example 1 without taxonomic id'),
            MagicMock(description='example 2 without taxonomic id')
        ]

        with patch.object(self.fasta_import.lineage_file_loader,
                          'load_lineage_file',
                          return_value='dummy_lineage.tsv') as mock_load_lineage_file:
            with patch.object(self.fasta_import, 'add_taxids') as mock_add_taxids:
                self.fasta_import.read_fasta(self.fasta_file) # Use real fasta file so
                #FileNotFoundError isn't raised

                mock_add_taxids.assert_not_called()

        mock_print.assert_called_once_with("No TaxIDs found in the fasta file. " +
                                           "Starting lineage file upload process.")

        mock_load_lineage_file.assert_called_once()

        self.assertEqual(self.fasta_import.lineage_file, 'dummy_lineage.tsv')

    def tearDown(self):
        """
        Tear down the test classes.
        """
        del self.fasta_import

class TestLinageFileLoader(unittest.TestCase):
    """
    Class to test the LineageFileLoader class.
    """
    def setUp(self):
        """
        Initialize the test class.
        """
        self.lineage_loader = LineageFileLoader()
        self.fasta_file_no_taxid = "data/test_data/fasta_example_file.fasta"
        self.fasta_import = CustomFastaImport(None)
        self.lineage_file = 'dummy_lineage_file.tsv'


    def test_validate_file_nofile(self):
        """
        Test if validate_file returns output when the input file does not exist.
        """

    @patch('builtins.print')
    def test_help_message(self, mock_print):
        """
        Test if the help message is displayed correctly.
        """
        self.lineage_loader._print_help_message()

        expected_message = self.lineage_loader.help_message_template.format(
            columns=self.lineage_loader.str_requirements)
        mock_print.assert_called_once_with(expected_message)

if __name__ == '__main__':
    unittest.main()
