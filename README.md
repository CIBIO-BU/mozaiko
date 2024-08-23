
# DNAquaIMG

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Conda Build](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml)
[![Lint Status](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/DNAquaIMG/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/DNAquaIMG)

![alt text](data/images/mosaiko-logo.png)
mosaiko: Piecing Together Complete Genetic Coverage for Biomonitoring

## Installation instructions

### Prerequisites

- Python 3.x
- Conda (Miniconda or Anaconda)
- Git

### Installation

1. Clone the repository:

   ```bash
   git clone git@github.com:CIBIO-BU/DNAquaIMG.git
   cd DNAquaIMG
   ```

2. Run the installation script:

   ```bash
   chmod +x conda_env_setup.sh
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

To assign taxonomic information to sequences,
the best practice is to include a JSON file that
specifies all the correct parameters.
A template for the JSON file can be found in
[here](https://github.com/CIBIO-BU/DNAquaIMG/blob/main/src/reference_database/assign_tax_parameters.json).

   ```bash
   mosaiko --assign_tax --json_file path/to/json/file
   ```

To dereplicate sequences in the reference database, the --dereplicate command can be used.
Similarly to --assign_tax, the best practice is to provide a JSON file that specifies all the
correct parameters. The input file for this task needs to be a .tsv file.
A template for the JSON file can be found in [here](https://github.com/CIBIO-BU/DNAquaIMG/blob/main/src/reference_database/dereplicate_parameters.json).

 ```bash
   mosaiko --dereplicate --json_file path/to/json/file
 ```

## Contacts

In case of enquiry, please reach out to <bu@cibioup.pt>.

## Issues/Problems

   ```markdown
   Where the issue tracking is i.e. github issues
   API documentation (f applicable)
   ```

## Dev/test/staging/production environments (when applicable)

institution-name/repository-name
