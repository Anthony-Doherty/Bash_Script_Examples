#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N align_CUT_And_RUN_Rx_RL9.sh  
#$ -m beas
#$ -M anthony.doherty@ed.ac.uk            
#$ -cwd                  
#$ -l h_rt=24:00:00 
#$ -l h_vmem=19G
#$ -pe sharedmem 16
#$ -t 1-14

#Set the THREADS variable
THREADS="16"

#Script name: align_CUT_And_RUN_Rx_RL9.sh   

#Aim: Align ChIPRX fastq files to hg_38_mm10 index > Separate reference and spike in reads > remove duplicates from paired end alignments > tidy up redundant intermediate files > generate alignment report > index the bam files

#Usage example: qsub ~/Central_Script_Repo/align_CUT_And_RUN_Rx_RL9.sh   

#Note: In order to use 1 input variable for paired-end reads specify the common root of the fastq file.
#This is enough as the bowtie command includes _1 _2 .

#make output directories

if [ ! -d analysis/ ]; then
    mkdir -p analysis/reference analysis/spikein
fi

#Stop on errors.

set -uex

#The reference genome indexes

REF38="/exports/igmm/eddie/gbrien-lab/index/hg38_mm10/Bowtie2/hg38_mm10"

# Target directory

TARGET_DIR="."

# Check if SampleIDs_Minimal.txt exists

if [ ! -f "$TARGET_DIR/SampleIDs.txt" ]; then
    echo "SampleIDs.txt not found in $TARGET_-14DIR"
    exit 1
fi

# Get list of files in target directory

files=($(cat "${TARGET_DIR}/SampleIDs.txt"))

# Get file to be processed by *this* task 
# extract the Nth file in the list of files, $files, where N == $SGE_TASK_ID

this_file="${files[$SGE_TASK_ID - 1]}" # Adjusted to use zero-based indexing
echo "Processing file: ${this_file} on $HOSTNAME"

# Initialise the environment modules

. /etc/profile.d/modules.sh

#Load required Modules

module load igmm/apps/python/3.12.3
module load igmm/apps/samtools/1.20
module load igmm/apps/bowtie/2.5.3

#Align paired-end reads to hg38 (and mm10) and output a sorted BAM file.

bowtie2 --very-sensitive -x $REF38 -p $THREADS -1 ${this_file}_1.trimmed.fq.gz -2 ${this_file}_2.trimmed.fq.gz | samtools view -@ $THREADS -b -S | samtools sort -@ $THREADS > ${this_file}.bam

##Separate mouse (spikein) reads

samtools view ${this_file}.bam -@ $THREADS -q 2 -h | grep mm10 | samtools view -@ $THREADS -b | samtools view -@ $THREADS -h | sed 's/mm10_//g' | samtools view -@ $THREADS -b > analysis/spikein/${this_file}_mm10.bam

#Separate human (reference) reads

samtools view ${this_file}.bam -@ $THREADS -q 2 -h | grep -v mm10 | samtools view -@ $THREADS -b > analysis/reference/${this_file}_hg38.bam

##Reference genome remove duplicates
#Sort the sub-sampled bam file by name.

samtools sort -@ $THREADS -n -o analysis/reference/${this_file}_hg38.bam.namesort.bam analysis/reference/${this_file}_hg38.bam

#Add ms and MC tags for markdup to use later.

samtools fixmate -m -@ $THREADS analysis/reference/${this_file}_hg38.bam.namesort.bam analysis/reference/${this_file}.fixmate.bam

#Markdup needs position order.

samtools sort -@ $THREADS -o analysis/reference/${this_file}.positionsort.bam analysis/reference/${this_file}.fixmate.bam

#Remove duplicates.

samtools markdup -r -@ $THREADS analysis/reference/${this_file}.positionsort.bam analysis/reference/${this_file}.dedup.bam

##Spikein genome remove duplicates
#Sort the sub-sampled bam file by name.

samtools sort -@ $THREADS -n -o analysis/spikein/${this_file}_mm10.bam.namesort.bam analysis/spikein/${this_file}_mm10.bam

#Add ms and MC tags for markdup to use later.

samtools fixmate -m -@ $THREADS analysis/spikein/${this_file}_mm10.bam.namesort.bam analysis/spikein/${this_file}.fixmate.bam

#Markdup needs position order.

samtools sort -@ $THREADS -o analysis/spikein/${this_file}.positionsort.bam analysis/spikein/${this_file}.fixmate.bam

#Remove duplicates.

samtools markdup -r -@ $THREADS analysis/spikein/${this_file}.positionsort.bam analysis/spikein/${this_file}.dedup.bam

#Tidy up redundant intermediate files.

rm -f analysis/reference/*namesort*
rm -f analysis/reference/*fixmate*
rm -f analysis/reference/*positionsort*
rm -f analysis/spikein/*namesort*
rm -f analysis/spikein/*fixmate*
rm -f analysis/spikein/*positionsort*
