#!/bin/bash

#Load required module
module load igmm/apps/python/3.12.3

#Notes:
#663399 Purple SMARCA4
#ED1C24 Red SS18SSX
#00A651 Green P300


# === CONFIGURATION ===
TARGET="P300"
BED_FILE="ConsensusPeaks.bed"
BED_LABEL="ConsensusPeaks"
COLOR_HEX="#00A651"  # Example: forest green; you can change to any valid hex or named color

# Input files
DMSO_BW="DMS0_${TARGET}_Rep123.av.bw"
ACBI1_BW="ACBI1_${TARGET}_Rep123.av.bw"
MATRIX_OUTPUT="matrix_${TARGET}_${BED_LABEL}_test.tab.gz"
HEATMAP_OUTPUT="heatmap_${TARGET}_${BED_LABEL}.svg"

echo "Running test for target: $TARGET with BED: $BED_FILE and color: $COLOR_HEX"

# Run computeMatrix
computeMatrix reference-point \
  -S "$DMSO_BW" "$ACBI1_BW" \
  -R "$BED_FILE" \
  --referencePoint center \
  --samplesLabel "DMSO" "ACBI1" \
  --missingDataAsZero \
  -a 25000 -b 25000 \
  --numberOfProcessors 20 \
  -out "$MATRIX_OUTPUT"

# Run plotHeatmap
plotHeatmap -m "$MATRIX_OUTPUT" \
  --regionsLabel "$BED_LABEL" \
  --colorList "w,${COLOR_HEX}" \
  --legendLocation none \
  --plotTitle "${TARGET} Signal - ${BED_LABEL}" \
  --yAxisLabel 'r/CPM' \
  -out "$HEATMAP_OUTPUT"
