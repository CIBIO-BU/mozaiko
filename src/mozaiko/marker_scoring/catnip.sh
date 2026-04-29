#!/usr/bin/env bash
set -euo pipefail

INPUT_FILE=$1
MAPPING_NAME=$2
MAPPING_COLS=$3
THRESHOLD=$4

catnip -i "$INPUT_FILE" -M -o "$MAPPING_NAME"
catnip -i "$INPUT_FILE" -f "$MAPPING_NAME" -c "$MAPPING_COLS" -p "$THRESHOLD"