#!/bin/bash

# Get input args
INPUT_FILE=$1
MAPPING_FILE=$2
MAPPING_COLS=$3
THRESHOLD=$4

source ~/miniconda3/etc/profile.d/conda.sh
conda activate catnip

catnip -i "$INPUT_FILE" -f "$MAPPING_FILE" -c "$MAPPING_COLS" -p "$THRESHOLD"

conda deactivate