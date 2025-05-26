#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N Peak_Call_Array_RL9  
#$ -m beas
#$ -M anthony.doherty@ed.ac.uk            
#$ -cwd
#$ -l h_rt=24:00:00 
#$ -l rl9=false
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -t 1

# Stop on errors
set -uex

# Directory where the sample files are located (not the CWD)
ORIGINAL_DIR="./analysis/reference/"

# Current directory for the sample ID file
TARGET_DIR="."

# Check if SampleIDs_No_IgG.txt exists in the CWD
if [ ! -f "$TARGET_DIR/SampleIDs_No_IgG.txt" ]; then
    echo "SampleIDs_No_IgG.txt not found in $TARGET_DIR"
    exit 1
fi

# Read sample names from the SampleIDs_No_IgG.txt file in the CWD
files=($(cat "${TARGET_DIR}/SampleIDs_No_IgG.txt"))

# Get sample for this task
this_file="${files[$SGE_TASK_ID - 1]}"
echo "Processing sample: ${this_file} on $HOSTNAME"

# Load environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/MACS2/2.1.1  
module load igmm/apps/homer/4.10    
module load igmm/apps/python 
module load roslin/bedtools/2.29.2

# Check that Input_BAM is provided
if [ -z "${Input_BAM:-}" ]; then
    echo "Error: Input_BAM variable not set. Use qsub -v Input_BAM=/path/to/input.bam"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p Macs2output

# Define the path for the current sample's BAM file in the original directory
this_bam="${ORIGINAL_DIR}/${this_file}.dedup.bam"

# Check if the BAM file exists in the original directory
if [ ! -f "$this_bam" ]; then
    echo "Error: BAM file ${this_bam} not found in ${ORIGINAL_DIR}"
    exit 1
fi

# Call peaks with MACS2
macs2 callpeak -t "$this_bam" -c "${Input_BAM}" -f BAMPE -g hs -n "$this_file" --broad --broad-cutoff 0.05 --outdir Macs2output

# Create a basic BED file for HOMER
cut -f 1,2,3,4,6 Macs2output/${this_file}_peaks.broadPeak > Macs2output/${this_file}.broadPeakForHomer.bed

# Format for HOMER
awk '{print $1"\t"$2"\t"$3"\t"$4"\t\t"$5}' Macs2output/${this_file}.broadPeakForHomer.bed > Macs2output/${this_file}_correctFormat.broadPeakForHomer.bed

# Report overlap with blacklist regions
bedtools intersect -a Macs2output/${this_file}_correctFormat.broadPeakForHomer.bed -b ~/Blacklists/custom_hg38.blacklist.merge.bed -wa | wc -l > Macs2output/${this_file}.blacklistreport.txt

# Filter out blacklist regions
bedtools intersect -a Macs2output/${this_file}_correctFormat.broadPeakForHomer.bed -b ~/Blacklists/custom_hg38.blacklist.merge.bed -v > Macs2output/${this_file}.broadPeakForHomerFiltered.bed

# Annotate peaks with HOMER
annotatePeaks.pl Macs2output/${this_file}.broadPeakForHomerFiltered.bed hg38 > Macs2output/${this_file}.broadPeakForHomerFiltered.homerannotate
