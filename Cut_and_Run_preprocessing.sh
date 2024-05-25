#!/bin/bash

ml Trimmomatic
ml Bowtie2
ml SAMtools
ml Picard
ml BEDTools
ml deepTools


# Function to print usage
usage() {
  echo "Usage: $0 -s sample -r reference -t cutoff -o output_dir -f output_file -g genome_file"
  echo "  -s sample         Sample name"
  echo "  -r reference      Path to Bowtie2 index"
  echo "  -t cutoff         Cutoff value for SEACR"
  echo "  -o output_dir     Output directory"
  echo "  -f output_file    Output file for SEACR results"
  echo "  -g genome_file    Genome file for bamCoverage"
  exit 1
}

# Parse arguments
while getopts "s:r:t:o:f:g:" opt; do
  case $opt in
    s) sample="$OPTARG" ;;
    r) ref="$OPTARG" ;;
    t) cutoff="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    f) output_file="$OPTARG" ;;
    g) genome_file="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check if all arguments are provided
if [ -z "$sample" ] || [ -z "$ref" ] || [ -z "$cutoff" ] || [ -z "$output_dir" ] || [ -z "$output_file" ] || [ -z "$genome_file" ]; then
  usage
fi

# Define file paths
trimmed_r1="${output_dir}/${sample}_R1_trim.fastq.gz"
trimmed_r2="${output_dir}/${sample}_R2_trim.fastq.gz"
orphan_r1="${output_dir}/${sample}_R1_orphan.fastq.gz"
orphan_r2="${output_dir}/${sample}_R2_orphan.fastq.gz"
bowtie2_log="${output_dir}/${sample}.bowtie2.log"
sorted_bam="${output_dir}/${sample}.sorted.bam"
marked_bam="${output_dir}/${sample}.marked.bam"
metrics_file="${output_dir}/${sample}.sorted.metrics"
rmdup_bam="${output_dir}/${sample}.rmDup.bam"
filtered_bam="${output_dir}/${sample}.filtered.bam"
bed_file="${output_dir}/${sample}.bed"
clean_bed="${output_dir}/${sample}.clean.bed"
fragments_bed="${output_dir}/${sample}.fragments.bed"
bedgraph="${output_dir}/${sample}.fragments.bedgraph"
seq_depth_norm_bw="${output_dir}/${sample}.SeqDepthNorm.bw"

# Create output directory if it does not exist
mkdir -p $output_dir

# Step 1: Trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
  ${sample}/C4_s3_125_068_S15_L002_R1_001.fastq.gz \
  ${sample}/C4_s3_125_068_S15_L002_R2_001.fastq.gz \
  $trimmed_r1 $orphan_r1 $trimmed_r2 $orphan_r2 \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# Step 2: Alignment with Bowtie2 and sorting with Samtools
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
  -x $ref -p 18 -1 $trimmed_r1 -2 $trimmed_r2 2> $bowtie2_log | \
  samtools sort -@ 18 -O bam -o $sorted_bam

# Step 3: Mark duplicates with Picard
java -XX:ParallelGCThreads=10 -Djava.io.tmpdir=/tmp -jar picard.jar MarkDuplicates \
  QUIET=true INPUT=$sorted_bam OUTPUT=$marked_bam METRICS_FILE=$metrics_file \
  REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

# Step 4: Remove duplicates with Samtools
samtools view -h -b -F 1024 $sorted_bam > $rmdup_bam

# Step 5: Filter alignments with Samtools
samtools view -h -b -F 1804 -f 2 $sorted_bam > $filtered_bam

# Step 6: Convert BAM to BED, clean, and generate fragments
bedtools bamtobed -bedpe -i $sorted_bam > $bed_file
awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_file > $clean_bed
cut -f 1,2,6 $clean_bed | sort -k1,1 -k2,2n -k3,3n > $fragments_bed

# Step 7: Generate bedGraph with bamCoverage
bamCoverage --bam $filtered_bam -of "bedgraph" -o $seq_depth_norm_bw \
  --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads

# Step 8: Run SEACR; SEACR_1.3.sh is a wrapper script for SEACR
bash SEACR_1.3.sh $seq_depth_norm_bw $cutoff non stringent $output_file
