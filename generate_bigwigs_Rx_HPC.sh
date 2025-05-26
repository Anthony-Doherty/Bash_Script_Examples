#!/bin/bash

# Grid Engine options (lines prefixed with #$)
#$ -N generate_bigwigs_Rx_RL9.sh  
#$ -m beas
#$ -M anthony.doherty@ed.ac.uk            
#$ -cwd                  
#$ -l h_rt=24:00:00 
#$ -l h_vmem=16G
#$ -pe sharedmem 24
#$ -t 1-9

#Set the THREADS variable
THREADS="24"

#Script name: generate_bigwigs_Rx_RL9.sh

#Aim: Generate alignment stats for dedup.bam files > Calculate scale factors > generate log of scale factors > index dedup.bam files > Use values to make Rx scaled bigwig

#Usage example: 

##make output directories

if [ ! -d analysis/bigwig ]; then
    mkdir -p analysis/bigwig
fi

#Stop on errors.

set -uex

# Target directory

TARGET_DIR="."

# Check if SampleIDs_Minimal.txt exists

if [ ! -f "$TARGET_DIR/SampleIDs_No_IgG.txt" ]; then
    echo "SampleIDs_No_IgG.txt not found in $TARGET_DIR"
    exit 1
fi

# Get list of files in target directory

files=($(cat "${TARGET_DIR}/SampleIDs_No_IgG.txt"))

# Get file to be processed by *this* task 
# extract the Nth file in the list of files, $files, where N == $SGE_TASK_ID

this_file="${files[$SGE_TASK_ID - 1]}" # Adjusted to use zero-based indexing
echo "Processing file: ${this_file} on $HOSTNAME"

# Initialise the environment modules

. /etc/profile.d/modules.sh

#Load required Modules

module load igmm/apps/python
module load igmm/apps/samtools 

#Use samtools flagstat to calculate alignment stats

REFCOUNTS=$(samtools flagstat -@ $THREADS analysis/reference/${this_file}.dedup.bam | head -n1| cut -f1 -d ' ')
SPIKECOUNTS=$(samtools flagstat -@ $THREADS analysis/spikein/${this_file}.dedup.bam | head -n1| cut -f1 -d ' ')

#Use bc with matlib for floating points to do math in bash
SCALEFACTOR=$( echo "1/($SPIKECOUNTS/1000000)" |bc -l)
PERCSPIKE=$( echo "(1/$REFCOUNTS)*$SPIKECOUNTS" | bc -l)


# Define the log file
LOGFILE="analysis/bigwig/scaleFactorLog.txt"

# If the log file does not exist, add a header
if [ ! -f "$LOGFILE" ]; then
    echo -e "Sample\tRefCounts\tSpikeCounts\tPercSpike\tScaleFactor" > "$LOGFILE"
fi

# Generate a formatted log entry
printf "%s\t%s\t%s\t%.6f\t%.6f\n" "${this_file}" "${REFCOUNTS}" "${SPIKECOUNTS}" "${PERCSPIKE}" "${SCALEFACTOR}" >> "$LOGFILE"

#Index dedup.bam files

samtools index -@ $THREADS analysis/reference/${this_file}.dedup.bam

#Use scalefactors to make bigwigs

bamCoverage -of bigwig -b analysis/reference/${this_file}.dedup.bam -bs 10 --smoothLength 30 -e 200 -ignore ChrM --scaleFactor $SCALEFACTOR -o analysis/bigwig/${this_file}.Rx.bw -p $THREADS
