
# mozaiko: Piecing Together Complete Genetic Coverage for Biomonitoring

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lint Status](https://github.com/CIBIO-BU/mozaiko/actions/workflows/super-linter.yml/badge.svg)](https://github.com/CIBIO-BU/mozaiko/actions/workflows/super-linter.yml)
[![Packge Tests](https://github.com/CIBIO-BU/mozaiko/actions/workflows/python-test-check.yml/badge.svg)](https://github.com/CIBIO-BU/mozaiko/actions/workflows/python-test-check.yml)
[![codecov](https://codecov.io/gh/CIBIO-BU/mozaiko/graph/badge.svg?token=21eBYKePwR)](https://codecov.io/gh/CIBIO-BU/mozaiko)

![alt text](data/images/mozaiko-logo.png)

Mozaiko is a bioinformatics tool for ranking primer sets for biomonitoring studies. It evaluates primer performance based on reference database representation, primer-binding efficiency, and taxonomic coverage for a user-defined list of taxa, helping researchers identify the most suitable primer sets for reliable and comparable genetic marker analyses.

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

### **Category 1:** Reference Database Quality

- **_barcoded_taxa_one_plus_**: percentage of taxa in OTL with more than one barcode. A barcode must include the target insert to be considered.
- **_ratio_barcoded_taxa_**: proportion of taxa in OTL with high barcode coverage (more than five barcodes) relative to taxa with minimal barcode coverage (at least one barcode). The value ranges from 0 to 1, 1 representing the optimal scenario.

### **Category 2:** Binding

- **_mismatch_score_**: the maximum number of mismatches between the forward primer and its binding site and the reverse primer and its binding site is recorded for each taxon. The maximum mismatch values are then summed to provide the score for the OTL list. The lowest values indicate lower mismatches between primer and primer-binding sites, facilitating amplification.
- **_priming_ratio_sum_**: sum of the priming ratio across taxon. The priming ratio is computed as the ratio of the maximum number of mismatches at the 3’ end of the primer binding site to the maximum number of mismatches across the entire primer binding site. The lowest values indicate fewer mismatches at the 3’ end of the primer binding site, hence higher binding strength.
- **_gc_clamp_score_**:  sum of GC matches at the 3’ end (GC Clamp) across all taxa present in the OTL. Higher values are preferable, as a content of 40-60% of GC matches promotes binding.
- **_min_tm_cv_**: The minimum melting temperature (Tm) between each pair of forward and reverse primers is calculated for each taxon. The coefficient of variation across taxa is then determined. Lower values indicate a more consistent thermal performance and are preferable.

### **Category 3:** Traits and Resolution

- **_taxonomic_resolution_**: percentage of taxa whose genetic divergence is higher than 2%. Higher values are preferable as they indicate an increased possibility of distinguishing between closely related taxa.
- **_resolution_ratio_**: percentage of taxa with genetic divergence higher than a cutoff (default cutoffs are 10%, 5%, and 2% for families, genus, and species, respectively), divided by the total number of taxa considered. This metric indicates the primer's ability to distinguish the target taxonomy from non-target taxa.



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

- The tool assumes that the Reference Database used as input, as well as the taxonomic indetification declared are correct.
- Primer rankings are always relative to a specific run, if different primers are given the results will vary.


```mermaid
flowchart TD

    %% ─────────────────────────────────────────────
    %% INPUTS+OUTPUTS
    %% ─────────────────────────────────────────────
    FASTA@{ shape: lean-r, label: "Harmonised Reference Database\n(FASTA)"}
    PRIMERS@{ shape: lean-r, label: "Primer set table\n(TSV)"}
    OTL@{ shape: lean-r, label: "Operational Taxonomic List (OTL)\n(TSV)"}
    RNAME{{"Run Name\n (String)"}}
    HRMRFDB@{ shape: lean-r, label: "Processed Reference Database\n (Dataframe)"}
    NMM{{"Max. Number of Mismatches (int)"}}
    MPI{{"Min. Percentage Identity\n(float, 0-1)"}}
    MAC{{"Min. Alignment Coverage\n(int, 0-1010)"}}
    MAP{{"Max. % of Ambigous Bases\n(int, 0-1)"}}
    MAXILUM{{"Max. Read Length According to Illumina Equipment\n(bool)"}}
    AMP@{ shape: lean-r, label: "Amplicons"}
    INST@{ shape: lean-r, label: "Inserts"}
    INCMP@{ shape: lean-r, label: "Amplicons w/ Incomplete PBS"}
    COMP@{ shape: lean-r, label: "Amplicons w/ Complete PBS"}

    %% ─────────────────────────────────────────────
    %% STAGE 1 — DATA IMPORT
    %% ─────────────────────────────────────────────
    subgraph S1["Module 1: Reference Database Pre-Processing"]
        LD["Load Reference Database"]
        VAL["Validate Input File"]
        CLNHD["Clean FASTA Headers"]
        EXTTX["Extract Sequence ID, Sequence, Lenght and Taxonomy"]
        HARMONIZE["Pre-process Harmonised Taxonomic Information"]
    end

    FASTA --> S1
    LD --> VAL --> CLNHD --> EXTTX --> HARMONIZE --> HRMRFDB

    %% ─────────────────────────────────────────────
    %% STAGE 2 — IN SILICO AMPLIFICATION
    %% ─────────────────────────────────────────────
    subgraph S2["Module 2: in-silico PCR Amplification"]
        LDPT["Load and Validate Primer Table"]
        CRTADAPT["Create Adapter Sequence"]
        SETUP["Setup Output Directories\n(Amplicon, Insert, Complete PCR, Incomplete PCR)"]
        INSIL["Run in-silico PCR\n (cutadapt)"]
        LOOP@{ shape: dbl-circ, label: "For each\nprimer pair …"}
        PGA["PGA Alignment\n(Retrieve Amplicons w/ Complete and Incomplete PBS)"]
        FILTER_PBS["Filter Amplicons by % of Ambiguous Bases"]
        AMP@{ shape: lean-r, label: "Amplicons"}
        PBS@{ shape: lean-r, label: "PBS"}
    end

    HRMRFDB --> S2
    PRIMERS --> LDPT
    LDPT --> CRTADAPT
    CRTADAPT --> SETUP
    SETUP --> LOOP
    LOOP --> INSIL
    RNAME --> INSIL
    NMM --> INSIL
    MPI --> PGA
    MAC --> PGA
    MAXILUM --> LDPT
    INSIL --> AMP
    INSIL --> INST
    INST --> PGA
    INSIL --> PBS
    PGA --> FILTER_PBS
    MAP --> FILTER_PBS
    FILTER_PBS --> INCMP
    FILTER_PBS --> COMP


    %% ─────────────────────────────────────────────
    %% STAGE 3 — METRICS EVALUATION
    %% ─────────────────────────────────────────────
    subgraph S3["Module 3: Metrics System"]

        LOAD_OTL["Read and Validate OTL"]
        OTL_PROC["Proccess OTL\n• Remove non-harmonised entries\n•Remove non-ASCII characters\n•Only keep entries at family, genus, species ranks"]
        OTL_IN@{ shape: lean-r, label: "Proccessed OTL"}

        subgraph CAT1["Category 1"]
            RDQ["Reference Database Quality"]
            M1A["Barcoded Taxa\n(% taxa with ≥1 barcode)"]
            M1B["Ratio of Barcoded Taxa\n(high-coverage / low-coverage)"]
            RDQ --> M1A & M1B
        end

        subgraph CAT2["Category 2"]
            BIND["Binding Efficency"]
            M2A["Mismatch Score\n(sum of max. mismatches)"]
            M2B["Priming Ratio Sum\n(3′ mismatches / total)"]
            M2C["GC Clamp Score\n(GC matches at 3′ end)"]
            M2D["Coefficient of Variation of the Minimum Melting Temperature"]
            BIND --> M2A & M2B & M2C & M2D
        end

        subgraph CAT3["Category 3"]
            CATNIP["Taxonomic Resolution"]
            M3A["Taxonomic Resolution\n(% taxa divergence > cutoff)"]
            M3B["Resolution Ratio\n(divergence > cutoff / total taxa)"]
            CATNIP --> M3A & M3B
        end

        RANK["Rank Primers"]
        RNKMD{{"Ranking Mode (&quot;flat&quot; or &quot;category&quot;, string)"}}
        FLAT["Flat\nAll metrics weighted equally"]
        CATW["Categorical\nWeighted by category"]
        SCINTM{{"Save Intermediate Ranks (bool)"}}
        INTER["Save Metric Ranks Used to Compute Primer Ranks"]
        SORT["Aggregate metric values per OTL taxonomy"]

        OTL_IN --> RDQ & BIND & CATNIP

    end

    OTL --> LOAD_OTL
    LOAD_OTL --> OTL_PROC
    OTL_PROC --> OTL_IN
    INST & COMP & INCMP --> RDQ
    INST --> BIND
    INST & INCMP --> CATNIP

    %% ─────────────────────────────────────────────
    %% STAGE 4 — RANKING
    %% ─────────────────────────────────────────────


    M1A & M1B --> RANK
    M2A & M2B & M2C & M2D --> RANK
    M3A & M3B --> RANK

    RANK --> SCINTM
    SCINTM --> INTER
    RANK --> SORT
    RANK --> RNKMD
    RNKMD --> FLAT
    RNKMD --> CATW


    %% ─────────────────────────────────────────────
    %% OUTPUT
    %% ─────────────────────────────────────────────
    FLAT & CATW --> RESULT["Ranked primer set per OTL\n{country}_ranked_primers.tsv"]

    %% ─────────────────────────────────────────────
    %% STYLING
    %% ─────────────────────────────────────────────
    classDef inout     fill:#34eb8f,stroke:#000000,color:#085041,font-weight:500
    classDef param  fill:#FFAB03,stroke:#000000,color:#085041
    classDef stage2    fill:#fff8e6,stroke:#000000,color:#633806
    classDef stage3    fill:#f5f0ff,stroke:#000000,color:#26215c
    classDef stage4    fill:#fef6f0,stroke:#000000,color:#4a1b0c
    classDef metric    fill:#eef9f4,stroke:#000000,color:#085041,font-size:12px
    classDef decision  fill:#faeeda,stroke:#000000,color:#412402
    classDef explain  fill:#34eb8f,stroke:#000000,color:#085041
    classDef metpy  fill:#e7f268,stroke:#000000,color:#085041

    class FASTA,PRIMERS,OTL,RNAME,RESULT,HRMRFDB,NMM,AMP,INST,INCMP,COMP,OTL_IN inout
    class NMM,RNAME,MPI,MAC,MAP,MAXILUM,RNKMD,SCINTM param
    class RDQ,BIND,CATNIP stage3
    class M1A,M1B,M2A,M2B,M2C,M2D,M3A,M3B metric
```

## Contacts

In case of enquiry, please reach out to <bu@biopolis.up.pt>.
