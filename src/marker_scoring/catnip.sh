#!/bin/bash
set -euo pipefail

OUTPUT_DIR=$1      # primer output directory
INPUT_FILE=$2
MAPPING_NAME=$3
MAPPING_COLS=$4
THRESHOLD=$5

conda activate catnip

# move into primer directory
cd "$OUTPUT_DIR" || exit 1

catnip -i "$INPUT_FILE" -M -o "$MAPPING_NAME"
catnip -i "$INPUT_FILE" -f "$MAPPING_NAME" -c "$MAPPING_COLS" -p "$THRESHOLD"

conda deactivate