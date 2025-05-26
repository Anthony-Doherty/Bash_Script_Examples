#!/bin/bash

# Name of the log file
LOG_FILE="peak_counts.log"

# Initialize or clear the log file
echo "Sample Name, Number of Peaks" > "$LOG_FILE"

# Loop through each .broadPeakForHomerFiltered.bed file
for file in *.broadPeakForHomerFiltered.bed; do
    if [[ -f "$file" ]]; then
        # Count the number of lines (peaks) in the file
        peak_count=$(wc -l < "$file")
        # Append the file name and peak count to the log file
        echo "$(basename -- "$file"), $peak_count" >> "$LOG_FILE"
    fi
done

echo "Peak counts logged to $LOG_FILE."
