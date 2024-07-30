"""
Unit tests for the CustomFastaImport class.
"""
import unittest
import os
from src.reference_database.sequence_import import CustomFastaImport


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
        self.data_transformer = CustomFastaImport(None)

    def test_empty_file(self):
        """
        Test if test_empty_file raises a ValueError when the input file is empty.
        """
        empty_file = os.path.join(self.data_dir, "empty_fasta.fasta")
        with self.assertRaises(ValueError) as context:
            self.data_transformer.read_fasta(empty_file)
        self.assertEqual(str(context.exception), 'Input file is empty.')

    def test_txt_file(self):
        """
        Test if test_txt_file raises a ValueError when the input file is not a fasta file.
        """
        txt_file = os.path.join(self.data_dir, "not_fasta_example.txt")
        with self.assertRaises(ValueError) as context:
            self.data_transformer.read_fasta(txt_file)
        self.assertEqual(str(context.exception), 'Input file must be a fasta file.')
        print('Test #1 passed.')

    def test_missing_file(self):
        """
        Test if test_missing_file raises a FileNotFoundError when the input file does not exist.
        """
        missing_file = os.path.join(self.data_dir, "no_file.fasta")
        with self.assertRaises(FileNotFoundError) as context:
            self.data_transformer.read_fasta(missing_file)
        self.assertEqual(str(context.exception), 'Input file does not exist in the directory.')
        print('Test #2 passed.')

    def test_wrong_input(self):
        """
        Test if test_wrong_input raises a ValueError when the input is not a string.
        """
        not_string = 1
        with self.assertRaises(ValueError) as context:
            self.data_transformer.read_fasta(not_string)
        self.assertEqual(str(context.exception), 'Directory must be a string.')
        print('Test #3 passed.')

    def test_read_fasta(self):
        """
        Test if read_fasta reads the fasta file correctly.
        """
        data = self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.assertEqual(data.shape, (3, 4))
        print('Test #4 passed.')

    def test_get_number_of_sequences(self):
        """
        Test if get_number_of_sequences returns the correct number of sequences.
        """
        self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.data_transformer.get_number_of_sequences(), 3)
        print('Test #5 passed.')

    def test_get_sequence_lengths(self):
        """
        Test if get_sequence_lengths returns the correct lengths of the sequences.
        """
        self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.data_transformer.get_sequence_lengths(), [16, 19, 8])
        print('Test #6 passed.')

    def test_get_sequence_ids(self):
        """
        Test if get_sequence_ids returns the correct sequence IDs.
        """
        self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.data_transformer.get_sequence_ids(),
                         ['CM074756.1', 'NC_088426.1', 'PP475397.1'])
        print('Test #7 passed.')

    def test_get_sequences(self):
        """
        Test if get_sequences returns the correct sequences.
        """
        self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.assertEqual(self.data_transformer.get_sequences(),
                         ['GTTATTGTAGCTTATC', 'GCATAAAGCATGGCACTGA', 'GTTATTGA'])
        print('Test #8 passed.')

    def test_df2fasta_structure(self):
        """
        Test if df2fasta correctly writes the data to a fasta file.
        """
        self.data_transformer.read_fasta(self.fasta_taxid_file)
        self.data_transformer.df2csv()
        example_file = 'processed_input_fasta.csv'
        with open(example_file, 'r', encoding='UTF-8') as f:
            lines = f.readlines()
            self.assertEqual(lines[0], 'seq_id,sequence,length,taxid\n')
            self.assertEqual(lines[1], 'CM074756.1,GTTATTGTAGCTTATC,16,8481\n')
            self.assertEqual(lines[2], 'NC_088426.1,GCATAAAGCATGGCACTGA,19,12345\n')
            self.assertEqual(lines[3], 'PP475397.1,GTTATTGA,8,106731\n')
        os.remove(example_file)

        print('Test #10 passed.')

    def tearDown(self):
        """
        Tear down the test class.
        """
        del self.data_transformer


if __name__ == '__main__':
    unittest.main()
