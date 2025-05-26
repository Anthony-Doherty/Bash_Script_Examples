#!/bin/bash

#Load module
module load igmm/apps/python/3.12.3

# Notes:
# Use this script to summarize signal over regions using multiBigwigSummary
# Good for PCA plots or correlation across multiple samples

# === CONFIGURATION ===
TARGET="P300"
BED_FILE="complete_SS18_Mega_bed_merged.bed"
DMSO_BW="DMS0_${TARGET}_Rep123.av.bw"
ACBI1_BW="ACBI1_${TARGET}_Rep123.av.bw"
LABELS="DMSO_FP ACBI1_FP"
OUTPUT_MATRIX="scores_per_transcript_mega_${TARGET}.npz"
OUTPUT_TABLE="scores_per_transcript_mega_${TARGET}.tab"

echo "Running multiBigwigSummary for target: $TARGET"
echo "Using BED file: $BED_FILE"
echo "BigWig files: $DMSO_BW, $ACBI1_BW"
echo "Labels: $LABELS"

# Validate inputs
if [[ ! -f "$BED_FILE" ]]; then
  echo "ERROR: BED file not found: $BED_FILE"
  exit 1
fi

if [[ ! -f "$DMSO_BW" || ! -f "$ACBI1_BW" ]]; then
  echo "ERROR: One or more BigWig files are missing."
  exit 1
fi

# Run multiBigwigSummary
multiBigwigSummary BED-file \
  --bwfiles "$DMSO_BW" "$ACBI1_BW" \
  --BED "$BED_FILE" \
  --labels $LABELS \
  -out "$OUTPUT_MATRIX" \
  --outRawCounts "$OUTPUT_TABLE"

echo "multiBigwigSummary complete. Outputs:"
echo "  Matrix: $OUTPUT_MATRIX"
echo "  Raw table: $OUTPUT_TABLE"
