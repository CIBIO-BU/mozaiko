from utils.general_utils import calculate_ambiguous_percentage, write_filtered_sequence
from utils.general_utils import calculate_ambiguous_percentage, write_filtered_sequence
from collections import defaultdict
from Bio import SeqIO

class PreProcessor:
    """
    Class to pre-process DNA sequences.
    """

    def __init__(self, max_ambiguous_percentage=0.05, extract_primers=False): # TODO: check extract_primers usage. Don't we always want to extract primers?
        """
        Initialize the PreProcessor object.

        Parameters:
        - max_ambiguous_percentage (float): Maximum allowed percentage of ambiguous bases in a sequence.
        - extract_primers (bool): Flag indicating whether to extract primers from the sequences.

        Returns:
        - None
        """

        self.max_ambiguous_percentage = max_ambiguous_percentage
        self.extract_primers = extract_primers
        self.primers = defaultdict(lambda: [None, None, None]) if extract_primers else None

    def detect_fwd_rev_primer_len(self, sequence):
        """
        Calculate the length of the forward and reverse primer.

        Parameters:
        - sequence (str): DNA sequence.

        Returns:
        - int: Length of the forward and reverse primer.
        """

        fwd_len = len(sequence) // 2
        rev_len = len(sequence) - fwd_len
        
        return fwd_len, rev_len
    
    def filter_sequences(self, input_file, out_file):
        """
        Filter DNA sequences based on the maximum allowed percentage of ambiguous bases.

        Parameters:
        - input_file (str): Path to the input file containing DNA sequences in FASTA format.
        - out_file (str): Path to the output file to write the filtered sequences.

        Returns:
        - dict or None: Dictionary of extracted primers if `extract_primers` is True, otherwise None.
        """
        
        with open(out_file, 'w') as output_handle:
            
            for record in SeqIO.parse(input_file, 'fasta'):
                sequence = str(record.seq)
                ambiguous_percentage = calculate_ambiguous_percentage(sequence)
                
                if ambiguous_percentage <= self.max_ambiguous_percentage:
                    
                    if self.extract_primers:
                        fwd_len, rev_len = self.detect_fwd_rev_primer_len(sequence)
                        accession_number = record.id
                        self.primers[accession_number][0] = sequence[:fwd_len] # stores forward primer
                        self.primers[accession_number][1] = sequence[-rev_len:] # stores reverse primer
                        self.primers[accession_number][2] = len(sequence) - fwd_len - rev_len # stores sequence lenght (total lenght minus lenght of both primers)
                    
                    write_filtered_sequence(output_handle, record)
        
        return self.primers if self.extract_primers else None