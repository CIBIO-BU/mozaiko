"""
This module contains the DataTransformer class which is responsible for handling and 
transforming data into desired outputs.

The DataTransformer class contains the following methods:
- read_fasta: Reads a fasta file.
- get_number_of_sequences: Returns the number of sequences in the data.
- get_sequence_lengths: Returns the lengths of the sequences in the data.
- get_sequence_ids: Returns the sequence IDs in the data.
- get_sequences: Returns the sequences in the data.
- df2csv: Write the data frame to a csv file.
"""
import os
from Bio import SeqIO
import pandas as pd

class DataTransformer:
    """
    Handles and transforms data into desired outputs.
    """

    def __init__(self, data):
        self.data = data

    def read_fasta(self, input_file):
        """
        Reads a fasta file.
        
        Parameters
        input_file (str): Path to the fasta file.
            
        Returns
        pd.DataFrame
        """

        if not isinstance(input_file, str):
            raise ValueError('Directory must be a string.')

        if not os.path.exists(input_file):
            raise FileNotFoundError('Input file does not exist in the directory.')

        if not input_file.endswith('.fasta'):
            raise ValueError('Input file must be a fasta file.')

        if os.path.getsize(input_file) == 0:
            raise ValueError('Input file is empty.')

        with open(input_file, 'r', encoding='UTF-8'):
            fasta_sequences = SeqIO.parse(input_file, 'fasta')

            data = []

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                seq_len = len(sequence)

                data.append([name, sequence, seq_len])

            self.data = pd.DataFrame(
                data, columns=['seq_id', 'sequence', 'lenght']
                )

            return self.data

    def get_number_of_sequences(self):
        """
        Returns the number of sequences in the data.
        
        Returns
        int
        """

        return len(self.data)

    def get_sequence_lengths(self):
        """
        Returns the lengths of the sequences in the data.
        
        Returns
        pd.Series
        """

        seq_lengths = self.data['lenght']
        seq_lengths = seq_lengths.to_list()

        return seq_lengths

    def get_sequence_ids(self):
        """
        Returns the sequence IDs in the data.
        
        Returns
        pd.Series
        """

        seq_ids = self.data['seq_id']
        seq_ids = seq_ids.to_list()

        return seq_ids

    def get_sequences(self):
        """
        Returns the sequences in the data.
        """

        seq_list = self.data['sequence']
        seq_list = seq_list.to_list()

        return seq_list

    def df2csv(self, output_name: str = 'processed_input_fasta.csv'):
        """
        Write the data frame to a csv file.
        """

        self.data.to_csv(output_name, index=False)
