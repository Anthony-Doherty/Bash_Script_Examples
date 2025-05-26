#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N Consensus_peaks.sh  
#$ -m eas
#$ -M anthony.doherty@ed.ac.uk            
#$ -cwd                  
#$ -l h_rt=24:00:00 
#$ -l h_vmem=24G
#$ -t 1

#Works really well. 

# Aim: Run the MSPC software for multiple prefixes automatically.

# Load required Modules
. /etc/profile.d/modules.sh
module load igmm/apps/mspc/6.0.0

# # Define an array of prefixes
# PREFIXES=(
#   "dCBP_FP"
#   "dCBP_P300_AM"
#   "dCBP_P300_CST"
#   "dCBP_CBP"
#   "dCBP_H3K27Ac"
#   "dCBP_H2BK20Ac"
# )

# Define an array of prefixes
PREFIXES=(
"ACBI1_P300"
"ACBI1_SMARCA4"
"ACBI1_SS18SSX"
"DMS0_P300"
"DMS0_SMARCA4"
"DMS0_SS18SSX"
)

# Define the input directory and output directory
INPUT_DIR="$(pwd)"  # Current directory
OUTPUT_DIR="ReplicatePeaks"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define MSPC parameters
RELATIONSHIP="Bio"
WEIGHT="1e-4"
SIGNIFICANCE="1e-8"
CUTOFF="2"
PARSER_CONFIG="parser_config.json"

# Loop through each prefix and process the files
for PREFIX in "${PREFIXES[@]}"; do
  echo "Processing prefix: $PREFIX"
  
  # Find replicate files matching the pattern
  INPUT_FILES=$(ls ${INPUT_DIR}/${PREFIX}_Rep_[123]*.broadPeak 2>/dev/null | xargs)

  # Check if the input files were found
  if [ -z "$INPUT_FILES" ]; then
    echo "No matching files found for prefix: $PREFIX with Rep_1, Rep_2, Rep_3 patterns."
    continue
  fi

  # Define output file name
  OUTPUT_FILE="${OUTPUT_DIR}/${PREFIX}.bed"

  # Build the MSPC command
  COMMAND="mspc -i ${INPUT_FILES} \
  -r ${RELATIONSHIP} \
  -w ${WEIGHT} \
  -s ${SIGNIFICANCE} \
  -o ${OUTPUT_FILE} \
  -c ${CUTOFF} \
  --parser ${PARSER_CONFIG}"

  # Execute the command
  echo "Running MSPC with the following command:"
  echo "$COMMAND"
  eval "$COMMAND"

  # Confirm completion for the prefix
  echo "MSPC processing completed for $PREFIX. Output saved to $OUTPUT_FILE."
done

# Final message
echo "All prefixes have been processed."
