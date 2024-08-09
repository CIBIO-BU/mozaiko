
# DNAquaIMG

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Conda Build](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/test-conda-setup.yml)
[![Lint Status](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/DNAquaIMG/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/DNAquaIMG/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/DNAquaIMG)

![alt text](<mosaiko-logo.png>)
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

# Contacts

In case of enquiry, please reach out to <bu@cibioup.pt>.

## Issues/Problems

   ```markdown
   Where the issue tracking is i.e. github issues
   API documentation (f applicable)
   ```

## Dev/test/staging/production environments (when applicable)

institution-name/repository-name
