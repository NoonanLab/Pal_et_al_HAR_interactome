#!/bin/bash

#SBATCH --job-name=Hicup_preprocessing
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=7G

# Load required modules
module load R
module load Bowtie2
module load SAMtools

# Full path to the conf file as defined in the HiCUP vignette
conf_file="$1"

# Create the "results" directory if it doesn't exist
mkdir -p results

# Change to the "results" directory
cd results

# Read the conf file and print its parameters
while IFS= read -r line; do
    echo "$line"
done < "$conf_file"

# Link to the folder where hicup package is downloaded
ln -s ~/bin/HiCUP-0.8.2

# Run all preprocessing steps of HiCUP with the specified config file
./hicup --config "$conf_file"

# Link to the folder where umi-tools package is downloaded
ln -s ~/bin/umi_tools

# Run deduplication on the filtered reads (here named R1_2.filt.sort.bam) to remove PCR duplicates based on UMI
umi_tools dedup -I R1_2.filt.sort.bam --paired --umi-separator=":" -S R1_2.dedup.bam

# Convert to the CHiCAGO input format (after downloading the chicagoTools package)
bam2chicago.sh R1_2.dedup.bam designDir/HARs_HGEs.baitmap designDir/GRCh38_DpnII.rmap HAR_HGE_chinput nodelete