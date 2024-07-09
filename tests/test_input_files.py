import unittest as unittest
from src.data_transformer import DataTransformer

class TestInputFiles(unittest.TestCase):
    
    def setUp(self):
        self.data_dir = "../data"
        self.txt_file = "not_fasta_example.txt"
        self.empty_file = "empty_fasta.fasta"
        self.fasta_file = "fasta_example_file.fasta"
        self.no_file = "no_file.fasta"
        self.not_string = 1

    def test_empty_file(self):
        with self.assertRaises(FileNotFoundError):
            self.DataTransformer.read_fasta(f"{self.data_dir}/{self.no_file}")

if __name__ == '__main__':
    unittest.main()
        
        