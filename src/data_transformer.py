from Bio import SeqIO
import pandas as pd
import re
import os

class DataTransformer(object):
    """
    Transforms data into desired outputs.
    """

    def __init__(self):
        self.data = None

    def read_fasta(self, input_file):
        """
        Reads a fasta file.
        
        Parameters
        input_file (str): Path to the fasta file.
            
        Returns
        pd.DataFrame
        """

        try:
            
            if not isinstance(input_file, str):
                raise ValueError('Input file must be a string.')
            
            if not os.path.exists(input_file):
                raise FileNotFoundError('Input file does not exist.')
            
            if not input_file.endswith('.fasta'):
                raise ValueError('Input file must be a fasta file.')
            
            if os.path.getsize(input_file) == 0:
                raise ValueError('Input file is empty.')
            
            fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

            data = []

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                seq_len = len(sequence)

                # TODO: Would we include a search for TaxID in the fasta header? Or require a TSV?
                taxid_match = re.search(r'(?<=taxid=)([0-9]+)(?=;)', fasta.description) # retrieved from Joana's code, doesn't extract anything

                if taxid_match:
                    taxid = taxid_match.group(1)
                else:
                    taxid = None

                data.append([name, sequence, seq_len, taxid])

            self.data = pd.DataFrame(data, columns=['SeqID', 'Sequence', 'Lenght', 'TaxID'])

            return self.data

        except Exception as e:
            print(e)