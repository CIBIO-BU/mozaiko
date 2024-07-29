"""
This module contains the CustomFastaImport class which is responsible for handling and
transforming custom data into desired inputs/outputs.

The CustomFastaImport class contains the following methods:
- _validate_input: Validates the input provided by the user.
- add_taxids: Joins the TaxIDs of sequences in the data.
- _check_for_taxids: Checks if fasta file contains TaxIDs.
- _request_lineage_file: Requests users to upload lineage file if no taxids are found in fasta file.
- read_fasta: Reads a fasta file.
- print_data: Print the DataFrame.
- get_taxids: Returns the TaxIDs of the sequences in the data.
- get_number_of_sequences: Returns the number of sequences in the data.
- get_sequence_lengths: Returns the lengths of the sequences in the data.
- get_sequence_ids: Returns the sequence IDs in the data.
- get_sequences: Returns the sequences in the data.
- df2csv: Write the data frame to a csv file.
"""
import os
from Bio import SeqIO
import pandas as pd
import tkinter as tk
import re
from tkinter.filedialog import askopenfilename


class CustomFastaImport:
    """
    Handles and transforms fasta into desired outputs.
    """

    def __init__(self, data):
        self.data = data

    def _validate_input(self, input_file):
        """
        Validades the input provided by the user.

        Parameters
        input_file (str): Path to fasta file.

        Raises
        ValueError: If the input is not a string.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the input file is not a fasta file.  
        ValueError: If the input file is empty.
        """

        if not isinstance(input_file, str):
            raise ValueError('Directory must be a string.')

        if not os.path.exists(input_file):
            raise FileNotFoundError('Input file does not exist in the directory.')

        if not input_file.endswith('.fasta'):
            raise ValueError('Input file must be a fasta file.')

        if os.path.getsize(input_file) == 0:
            raise ValueError('Input file is empty.')
        
    def add_taxids(self, input_file):
        """
        Joins the TaxIDs of sequences in the data.

        Returns
        pd.DataFrame
        """
            
        with open(input_file, 'r', encoding='UTF-8'):
            records = SeqIO.parse(input_file, 'fasta')

            taxids = []

            for seq in records:
                taxid = re.search(r'(?<=taxid=)([0-9]+)', seq.description).group(1)
                taxids.append(taxid)
                
            self.data['taxid'] = taxids

        return self.data
        
    def _check_for_taxids(self, input_file):
        """
        Checks if fasta file contains TaxIDs.
        """
        with open(input_file, 'r', encoding='UTF-8'):
            for record in SeqIO.parse(input_file, 'fasta'):

                if 'taxid' in record.description.lower():
                    
                    return True
                
                else:

                    self._request_lineage_file(input_file)
                
    def _request_lineage_file(self, input_file):
        """
        Requests users to upload lineage files whem taxids are not present in fasta file.

        Returns
        pd.DataFrame
        """

        header_requirements = ['seq_ID',
                                'Species',
                                'Genus',
                                'Family',
                                'Order',
                                'Class',
                                'Phylum',
                                'Subkingdom',
                                'Kingdom',
                                'Empire']
        
        str_requirements = ', '.join(header_requirements)

        print("File given as input does not contain Taxonomic IDs. \n" +
              "Please upload a TSV file containing sequence IDs and the respective lineage. \n" + 
              "The file should contain the following columns: [" + str_requirements + "]. \n" +
              "Please make sure columns follow the same format and order before uploading. \n" +
              "For the taxonomic levels where no assignment is available, "
              "please leave the cells blank. \n" +
              "\n" +
              "------ \n" +
              "\n" +
              "If your file contains Taxonomic IDs, please make sure these are present " +
              "in the header of each record of fasta file. Additionally, these must be "
              "identified with 'taxid=' beforehand. For example: 'CM074756.1;taxid=8481'.")
        
        root = tk.Tk()
        root.withdraw()
        
        while True:
            print("Please select the TSV file to upload.")

            file = askopenfilename(filetypes=[('TSV files', '*.tsv')])
            
            if not file:
                raise ValueError("No file selected. Please try again.")
            
            if not file.endswith('.tsv'):
                raise ValueError("File must be a TSV file. Please try again.")

            try:
                lineage_file = pd.read_csv(file, header=0, sep='\t')
            
            except Exception as e:
                raise ValueError("Error reading the file: {e}. Please try.")

            lineage_header = [column.lower() for column in lineage_file.columns]
            
            if lineage_header == header_requirements:
                print("File was successfully loaded.")
                root.destroy()
                
                self.lineage_file = lineage_file

                return lineage_file
                
            else:
                print("Columns in TSV file do not match the requirements: [" +
                             {header_requirements} + "]. Please try again.")

    def read_fasta(self, input_file):
        """
        Reads a fasta file.

        Parameters
        input_file (str): Path to the fasta file.

        Returnscheck_for_taxids(input_file)

        self.
        pd.DataFrame
        """

        self._validate_input(input_file)

        self._check_for_taxids(input_file)
        
        with open(input_file, 'r', encoding='UTF-8'):
            records = SeqIO.parse(input_file, 'fasta')

            data = []

            for seq in records:
                name, sequence = seq.id, str(seq.seq)
                seq_len = len(sequence)

                data.append([name, sequence, seq_len])

            self.data = pd.DataFrame(
                data, columns=['seq_id', 'sequence', 'lenght']
                )

        self.add_taxids(input_file)

        return self.data
        
    def print_data(self):
        """
        Print the DataFrame
        """
        return self.data

    def get_taxids(self):
        """
        Returns the TaxIDs of the sequences in the data.

        Returns
        pd.Series
        """

        taxids = self.data['taxid']
        taxids = taxids.to_list()

        return taxids

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
        
