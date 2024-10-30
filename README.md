
# DNAquaIMG

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Conda Build](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml)
[![Lint Status](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/DNAquaIMG/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/DNAquaIMG)

![alt text](data/images/mosaiko-logo.png)
mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring

## Installation instructions

### Prerequisites

- Python 3.x
- Conda (Miniconda or Anaconda)
- Git

### Installation

1. Clone the repository:

   ```bash
   git clone git@github.com:CIBIO-BU/DNAquaIMG.git
   ```
      ```bash
   cd DNAquaIMG
   ```

2. Run the installation script:

   ```bash
   chmod +x conda_env_setup.sh
   ```
      ```bash
   ./conda_env_setup.sh
   ```

This script will:

- Check if Conda is installed;
- Create a new Conda environment named "dnaquaimg", if it does not yet exist;
- Activate the Conda environment;
- Clone the DNAquaIMG repository, if not already cloned;
- Install the DNAquaIMG package.

## Running mosaiko

To upload custom FASTA sequences:

   ```bash
   mosaiko --load_custom_fasta --input path/to/file
   ```

To assign taxonomic information to sequences:

   ```bash
   mosaiko --assign_tax --json_file path/to/json/file
   ```
The best practice is to include a JSON file that
specifies all the correct parameters.
A template for the JSON file can be found in
[here](https://github.com/CIBIO-BU/DNAquaIMG/blob/main/src/reference_database/assign_tax_parameters.json). This command contemplates eight parameters, five are required ('input', 'output', 'acc2tax', 'taxid' and 'name') and three are optional ('web', 'ranks' and 'missing'). 'axx2tax', 'taxid' and 'name' are related to necessary taxonomic information to complete the task, Despite being required, these fields are already filled in and are available when the taxonomic files are downloaded upon calling the command (no user action required). The 'input' related to the input FASTA file and the 'output' to the output TSV file. For the optional ones, 'web' allows for a web-search to query NCBI's EFetch for missing taxonomic information; 'ranks' allows the user to choose which ranks should represent the organism lineage; 'missing' allows the user to provide a file path to write sequences where taxonomic information was not possible to be retrieved.


To dereplicate sequences in the reference database, the --dereplicate command can be used.

 ```bash
   mosaiko --dereplicate --json_file path/to/json/file
 ```

Similarly to --assign_tax, the best practice is to provide a JSON file that specifies all the correct parameters. A template for the dereplication JSON file can be found in [here](https://github.com/CIBIO-BU/DNAquaIMG/blob/main/src/reference_database/dereplicate_parameters.json). This command contemplates four parameters, two being required ('input' and 'output') and two being optional ('method', 'ranks'). Both the 'input' and 'output' should be a TSV file. For the optional ones, 'ranks' allows the user to choose which ranks should represent the organism lineage and 'method' allows the user to choose which method should be used for the dereplication. Please refer to [CRABS' original documentation](https://github.com/gjeunen/reference_database_creator/tree/main?tab=readme-ov-file#6-dereplicate) for further details.

To run the in-silico amplification analysis, the --in_silico_analysis command can be used.

 ```bash
   mosaiko --in_silico_analysis --input path/to/fasta/file
 ```

The in-silico analysis command will run the amplification process. It requires a FASTA file to be provided as input. This FASTA file should consider headers and sequences in different lines, with the header formatted as ">AB123 | taxid=1234". After running the command, the user will be prompted to provide a primer table and a name for the folder where results will be outputted. The primer table must be provided as a TSV file with following fields: ["target_group", "barcode_region" "assay_name", "fw_seq", "rev_seq"]. After all inputs are validated, the workflow will perform in-silico amplification, retriving amplicons and inserts where amplification was successful (< 3 mismatches) and inserts where amplification was unsuccessful due to too many mismatches (> 3 mismatches) or due to incomplete forward or reverse primer-binding sites.


## Contacts

In case of enquiry, please reach out to <bu@cibioup.pt>.

## Issues/Problems

   ```markdown
   Where the issue tracking is i.e. github issues
   API documentation (f applicable)
   ```

## Dev/test/staging/production environments (when applicable)

institution-name/repository-name
