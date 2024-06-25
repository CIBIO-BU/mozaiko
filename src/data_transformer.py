from Bio import SeqIO
import pandas as pd
import re

class DataTransformer:
    """
    Transforms data into desired outputs.
    """

    ## TODO: Review function. This method was transformed to a class without revision.

    def fasta_2_df(input_file):
        seq_object = SeqIO.parse(input_file, "fasta")
        sequences = []
        
        for seq in seq_object:
            sequences.append(seq)
        
        samples = []
        tax = []
        seq_lenghts = []
        sequences2 = []
        
        for record in sequences:
            seq = str(record.seq)
            length = len(seq)
            sample = record.id
            taxid = re.search(r'(?<=taxid=)([0-9]+)(?=;)',record.description).group(1)
            samples.append(sample)
            tax.append(taxid)
            seq_lenghts.append(length)
            sequences2.append(seq)
        
        data =  pd.DataFrame(list(zip(samples, tax, seq_lenghts, sequences2)), columns=["AC", "TaxID", "Length", "Sequence"])
        
        return(data)