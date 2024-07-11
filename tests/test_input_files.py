import unittest as unittest
import os
import sys

# Add root directory to the path in order to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.data_transformer import DataTransformer

class TestDataTransformer(unittest.TestCase):
    
    def setUp(self):
        self.DataTransformer = DataTransformer
        self.data_dir = "tests/test_data"
        self.txt_file = os.path.join(self.data_dir, "not_fasta_example.txt")
        self.empty_file = os.path.join(self.data_dir, "empty_fasta.fasta")
        self.fasta_file = os.path.join(self.data_dir, "fasta_example_file.fasta")
        self.missing_file = os.path.join(self.data_dir, "no_file.fasta")
        self.not_string = 1
        self.open_files = []

    def test_empty_file(self):
        data_transformer = DataTransformer(self.empty_file)
        with self.assertRaises(ValueError) as context:
            data_transformer.read_fasta(self.empty_file)
        self.assertEqual(str(context.exception), 'Input file is empty.')

    def test_txt_file(self):
        data_transformer = DataTransformer(self.txt_file)

        with self.assertRaises(ValueError) as context:
            data_transformer.read_fasta(self.txt_file)
        self.assertEqual(str(context.exception), 'Input file must be a fasta file.')
        print('Test #1 passed.')

    def test_missing_file(self):
        data_transformer = DataTransformer(self.missing_file)

        with self.assertRaises(FileNotFoundError) as context:
            data_transformer.read_fasta(self.missing_file)
        self.assertEqual(str(context.exception), 'Input file does not exist in the directory.')
        print('Test #2 passed.')
    
    def test_wrong_input(self):
        data_transformer = DataTransformer(self.not_string)

        with self.assertRaises(ValueError) as context:
            data_transformer.read_fasta(self.not_string)
        self.assertEqual(str(context.exception), 'Directory must be a string.')
        print('Test #3 passed.')
    
    def test_read_fasta(self):
        data_transformer = DataTransformer(self.fasta_file)
        data = data_transformer.read_fasta(self.fasta_file)
        self.assertEqual(data.shape, (3, 3))
        print('Test #4 passed.')

    def test_get_number_of_sequences(self):
        data_transformer = DataTransformer(self.fasta_file)
        data = data_transformer.read_fasta(self.fasta_file)
        self.assertEqual(data_transformer.get_number_of_sequences(), 3)
        print('Test #5 passed.')

    
    def test_get_sequence_lengths(self):
        data_transformer = DataTransformer(self.fasta_file)
        data = data_transformer.read_fasta(self.fasta_file)
        self.assertEqual(data_transformer.get_sequence_lengths(), [16, 19, 8])
        print('Test #6 passed.')

    def test_get_sequence_ids(self):
        data_transformer = DataTransformer(self.fasta_file)
        data = data_transformer.read_fasta(self.fasta_file)
        self.assertEqual(data_transformer.get_sequence_ids(), ['CM074756.1', 'NC_088426.1', 'PP475397.1'])
        print('Test #7 passed.')

    def test_get_sequences(self):
        data_transformer = DataTransformer(self.fasta_file)
        data = data_transformer.read_fasta(self.fasta_file)
        self.assertEqual(data_transformer.get_sequences(), ['GTTATTGTAGCTTATC', 'GCATAAAGCATGGCACTGA', 'GTTATTGA'])
        print('Test #8 passed.')
    
if __name__ == '__main__':
    unittest.main()
        
        