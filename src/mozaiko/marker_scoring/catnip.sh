#!/usr/bin/env bash
set -euo pipefail

INPUT_FILE=$1
MAPPING_NAME=$2
MAPPING_COLS=$3
THRESHOLD=$4
# Optional resource settings; left empty, catnip's own defaults are used.
THREADS=${5:-}
MEMORY=${6:-}          # MB
PARALLEL_JOBS=${7:-}   # requires catnip-seq >= 0.2.0

EXTRA_ARGS=()
[ -n "$THREADS" ] && EXTRA_ARGS+=(-t "$THREADS")
[ -n "$MEMORY" ] && EXTRA_ARGS+=(-m "$MEMORY")
[ -n "$PARALLEL_JOBS" ] && EXTRA_ARGS+=(-j "$PARALLEL_JOBS")

catnip -i "$INPUT_FILE" -M -o "$MAPPING_NAME"
catnip -i "$INPUT_FILE" -f "$MAPPING_NAME" -c "$MAPPING_COLS" -p "$THRESHOLD" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
