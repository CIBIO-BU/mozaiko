
# mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lint Status](https://github.com/CIBIO-BU/mozaiko/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/mozaiko/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/mozaiko/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/mozaiko/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/mozaiko/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/mozaiko)

![alt text](data/images/mozaiko-logo.png)

mozaiko is a bioinformatics tool designed to help researchers select optimized sets of primers for complete coverage in biomonitoring studies. Taking inspiration from mosaics, where small pieces fit together to form a whole, mozaiko supports comprehensive genetic marker analysis by ranking primers' fitness.

The name comes from the Esperanto word 'Mozaiko', reflecting the idea of bringing different elements together. With mozaiko, researchers can efficiently select primer sets for a range of applications, making biomonitoring and ecological studies more reliable and comparable.

## Installation instructions

### Prerequisites

- Python 3.x
- Conda (Miniconda or Anaconda)
- Git

### Installation

1. Clone the repository:

   ```bash
   git clone git@github.com:CIBIO-BU/mozaiko.git
   ```
      ```bash
   cd mozaiko
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

3. Run mozaiko:

      ```bash
   mozaiko --help
   ```

Th installation script will:

- Check if Conda is installed;
- Create a new Conda environment named "mozaiko", if it does not yet exist;
- Activate the Conda environment;
- Install the mozaiko package;
- Install required dependencies and tools.

## mozaiko Metrics' System

mozaiko contains three main categories to evaluate and rank primer sets:

### **Module 1:** Reference Database Quality

- **_barcoded_taxa_one_plus_**: percentage of taxa in OTL with more than one barcode. A barcode must include the target insert to be considered.
- **_ratio_barcoded_taxa_**: proportion of taxa in OTL with high barcode coverage (more than five barcodes) relative to taxa with minimal barcode coverage (at least one barcode). The value ranges from 0 to 1, 1 representing the optimal scenario.

### **Module 2:** Binding

- **_mismatch_score_**: the maximum number of mismatches between the forward primer and its binding site and the reverse primer and its binding site is recorded for each taxon. The maximum mismatch values are then summed to provide the score for the OTL list. The lowest values indicate lower mismatches between primer and primer-binding sites, facilitating amplification.
- **_priming_ratio_sum_**: sum of the priming ratio across taxon. The priming ratio is computed as the ratio of the maximum number of mismatches at the 3’ end of the primer binding site to the maximum number of mismatches across the entire primer binding site. The lowest values indicate fewer mismatches at the 3’ end of the primer binding site, hence higher binding strength.
- **_gc_matches_across_taxon_**:  sum of G-C matches at the 3’ end across all taxa present in the OTL. Higher values are preferable, as a content of 40-60% of GC matches promotes binding.
- **_min_tm_cv_**: The minimum melting temperature (Tm) between each pair of forward and reverse primers is calculated for each taxon. The coefficient of variation across taxa is then determined. Lower values indicate a more consistent thermal performance and are preferable.
- **_tm_score_**: proportion of taxa with a lower or equal variation of Tm below 2ºC.  Higher values are preferable as they indicate a better thermal performance across taxa in the OTL.
- **_amplification_success_percent_**: the ratio of taxa that amplify to the total number of taxa with sequences containing primer binding sites, expressed as a percentage. Higher values represent higher amplification success across taxa.

### **Module 3:** Traits and Resolution

- **_taxonomic_resolution_**: percentage of taxa whose genetic divergence is higher than 2%. Higher values are preferable as they indicate an increased possibility of distinguishing between closely related taxa.



The final ranking position is determined based on the individual ranking scores for each metric, presented in the output file intermediate_ranks, with all metrics weighted equally. Each metric is ranked based on whether higher or lower values are more desirable:
   - Descending (higher is better):
      - barcoded_taxa_one_plus
      - ratio_barcoded_taxa
      - gc_matches_across_taxon
      - tm_score
      - amplification_success_percent
   - Ascending (lower is better):
      - mismatch_score
      - priming_ratio_sum
      - min_tm_cv
      - taxonomic_resolution

For metrics ranked ascending, primers with lower values are preferred. For example, a lower ‘mismatch_score’ is better because it means fewer mismatches. For metrics ranked descending, primers with higher values are preferred.

## mozaiko Workflow

Primer rankings are always relative to a specific run, if different primers are given the results will vary.

## Contacts

In case of enquiry, please reach out to <bu@biopolis.up.pt>.
