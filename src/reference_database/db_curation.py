"""
This modules generates scripts to borrow methods from CRBAS
(https://github.com/gjeunen/reference_database_creator/tree/main) with the aim to retrieve taxonomic
information, assing taxonomic IDs, generate lineage files and dereplicate sequences.
"""

import json


class CrabsScriptGenerator:
    def __init__(self, json_file):
        self.json_file = json_file

    def load_parameters(self):
        with open(self.json_file) as file:
            self.params = json.load(file)

    def generate_assign_tax_script(self, sh_filename):
        self.load_parameters()

        pass

    def generate_dereplicate_script(self, sh_filename):
        pass


class DatabaseCuration:
    def __init__(self, fasta, lineage):
        self.fasta = fasta
        self.lineage_file = lineage

    pass
