
# mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lint Status](https://github.com/CIBIO-BU/mozaiko/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/mozaico/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/mozaiko/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/mozaico/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/mozaiko/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/mozaico)

![alt text](data/images/mosaiko-logo.png)

mozaiko is a bioinformatics tool designed to help researchers select optimized sets of primers for complete coverage in biomonitoring studies. Taking inspiration from mosaic art—where small pieces fit together to form a whole—mozaiko supports comprehensive genetic marker analysis by suggesting a fitting combination of primers.

The name comes from the Esperanto word 'Mozaiko', reflecting the idea of bringing different elements together. With mozaiko, researchers can efficiently select primer sets for a range of applications, making biomonitoring and ecological studies more reliable and comparable.

## Installation instructions

### Prerequisites

- Python 3.x
- Conda (Miniconda or Anaconda)
- Git

### Installation

1. Clone the repository:

   ```bash
   git clone git@github.com:CIBIO-BU/mozaico.git
   ```
      ```bash
   cd mozaico
   ```

2. Run the installation script:

   ```bash
   chmod +x conda_env_setup.sh
   ```
      ```bash
   ./conda_env_setup.sh
   ```

3. Activate the environment:

      ```bash
   conda activate mozaiko
   ```

Th installation script will:

- Check if Conda is installed;
- Create a new Conda environment named "mozaiko", if it does not yet exist;
- Activate the Conda environment;
- Clone the mozaico repository, if not already cloned;
- Install the mozaiko package;
- Install required dependencies and tools.

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
[here](https://github.com/CIBIO-BU/mozaico/blob/main/src/reference_database/assign_tax_parameters.json). This command contemplates eight parameters, five are required ('input', 'output', 'acc2tax', 'taxid' and 'name') and three are optional ('web', 'ranks' and 'missing'). 'axx2tax', 'taxid' and 'name' are related to necessary taxonomic information to complete the task, Despite being required, these fields are already filled in and are available when the taxonomic files are downloaded upon calling the command (no user action required). The 'input' related to the input FASTA file and the 'output' to the output TSV file. For the optional ones, 'web' allows for a web-search to query NCBI's EFetch for missing taxonomic information; 'ranks' allows the user to choose which ranks should represent the organism lineage; 'missing' allows the user to provide a file path to write sequences where taxonomic information was not possible to be retrieved.


To dereplicate sequences in the reference database, the --dereplicate command can be used.

 ```bash
   mosaiko --dereplicate --json_file path/to/json/file
 ```

Similarly to --assign_tax, the best practice is to provide a JSON file that specifies all the correct parameters. A template for the dereplication JSON file can be found in [here](https://github.com/CIBIO-BU/mozaico/blob/main/src/reference_database/dereplicate_parameters.json). This command contemplates four parameters, two being required ('input' and 'output') and two being optional ('method', 'ranks'). Both the 'input' and 'output' should be a TSV file. For the optional ones, 'ranks' allows the user to choose which ranks should represent the organism lineage and 'method' allows the user to choose which method should be used for the dereplication. Please refer to [CRABS' original documentation](https://github.com/gjeunen/reference_database_creator/tree/main?tab=readme-ov-file#6-dereplicate) for further details.

To run the in-silico amplification analysis, the --in_silico_analysis command can be used.

 ```bash
   mosaiko --in_silico_analysis --input path/to/fasta/file
 ```

The in-silico analysis command will run the amplification process. It requires a FASTA file to be provided as input. This FASTA file should consider headers and sequences in different lines, with the header formatted as ">AB123 | taxid=1234". After running the command, the user will be prompted to provide a primer table and a name for the folder where results will be outputted. The primer table must be provided as a TSV file with following fields: ["target_group", "barcode_region" "assay_name", "fw_seq", "rev_seq", "min_len", "max_len"]. After all inputs are validated, the workflow will perform in-silico amplification, retriving amplicons and inserts where amplification was successful (< 3 mismatches) and inserts where amplification was unsuccessful due to too many mismatches (> 3 mismatches) or due to incomplete forward or reverse primer-binding sites.


## Contacts

In case of enquiry, please reach out to <bu@biopolis.up.pt>.
