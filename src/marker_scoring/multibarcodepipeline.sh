#!/bin/bash

# Get input args
OUTPUT_FOLDER=$1
INPUT_FILE=$2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate multibarcode

echo $OUTPUT_FOLDER
echo $INPUT_FILE

multi-barcode -o "$OUTPUT_FOLDER" "$INPUT_FILE"

conda deactivate