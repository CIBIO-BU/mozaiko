from Bio import SeqIO
import pandas as pd
import re
import os

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
        
        
        with open(input_file, 'r') as input_file:
            fasta_sequences = SeqIO.parse(input_file, 'fasta')

            data = []

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                seq_len = len(sequence)

                data.append([name, sequence, seq_len])

            self.data = pd.DataFrame(
                data, columns=['SeqID', 'Sequence', 'Lenght']
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
        
        seq_lengths = self.data['Lenght']
        seq_lengths = seq_lengths.to_list()

        return seq_lengths
    
    def get_sequence_ids(self):
        """
        Returns the sequence IDs in the data.
        
        Returns
        pd.Series
        """
        
        return self.data['SeqID']
    
    def get_sequences(self):
        """
        Returns the sequences in the data.
        """

        return self.data['Sequence']
