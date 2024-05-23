#!/bin/bash

data_directory=$1
output_directory=$2
reference_genome_dir=$3
rsem_dir=$4
annotation_file_path=$5
threads=$6
adapter_file=$7

mkdir -p "$output_directory"

mkdir -p "$output_directory/fastqc"
fastqc -t 5 -o "$output_directory/fastqc" "$data_directory"/*.gz
mkdir -p "$output_directory/trimmed_files"

fastq_files=( "$data_directory"/*.fastq.gz "$data_directory"/*.fq.gz )
true_names_and_file=()
for f in "${fastq_files[@]}"; do
    if [[ $f == *r1.fastq.gz || $f == *r1.fq.gz ]]; then
        true_names_and_file+=( "${f%r1*}" "$f" )
    fi
done
names=()
for f in "${true_names_and_file[@]}"; do
    names+=( "${f%%.*}" )
done

for name in "${names[@]}"; do
    echo "Starting new analysis for $name"

    unzip "$output_directory/fastqc/${name}_r1_fastqc.zip"
    unzip "$output_directory/fastqc/${name}_r2_fastqc.zip"

    R1_F=$(grep -P -c "Adapter Content\tfail" "$output_directory/fastqc/${name}_r1_fastqc/fastqc_data.txt")
    R2_F=$(grep -P -c "Adapter Content\tfail" "$output_directory/fastqc/${name}_r2_fastqc/fastqc_data.txt")

    if [[ $R1_F -eq 0 && $R2_F -eq 0 ]]; then

        A="$data_directory/${name}_r1.fastq.gz"
        B="$data_directory/${name}_r2.fastq.gz"
        C="$output_directory/trimmed_files/${name}_r1_paired.fq.gz"
        D="$output_directory/trimmed_files/${name}_r1_unpaired.fq.gz"
        E="$output_directory/trimmed_files/${name}_r2_paired.fq.gz"
        F="$output_directory/trimmed_files/${name}_r2_unpaired.fq.gz"

        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads "$threads" "$A" "$B" "$C" "$D" "$E" "$F" ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    else
        mv "$data_directory/${name}_r1.fastq.gz" "$output_directory/trimmed_files/${name}_r1_paired.fq.gz"
        mv "$data_directory/${name}_r2.fastq.gz" "$output_directory/trimmed_files/${name}_r2_paired.fq.gz"
    fi

    mkdir -p "$output_directory/trimmed_files/fastqc"
    fastqc -t 5 -o "$output_directory/trimmed_files/fastqc" "$output_directory/trimmed_files/${name}_r1_paired.fq.gz"
    fastqc -t 5 -o "$output_directory/trimmed_files/fastqc" "$output_directory/trimmed_files/${name}_r2_paired.fq.gz"

    unzip "$output_directory/trimmed_files/fastqc/${name}_r1_paired_fastqc.zip"
    unzip "$output_directory/trimmed_files/fastqc/${name}_r2_paired_fastqc.zip"

    R1_F=$(grep -P -c "Adapter Content\tfail" "$output_directory/trimmed_files/fastqc/${name}_r1_paired_fastqc/fastqc_data.txt")
    R2_F=$(grep -P -c "Adapter Content\tfail" "$output_directory/trimmed_files/fastqc/${name}_r2_paired_fastqc/fastqc_data.txt")

    if [[ $R1_F -eq 0 && $R2_F -eq 0 ]]; then

        echo "Adapter still present after trimming. Do you want to continue? (y/n)"
        read answer

        if [[ $answer == "n" ]]; then
            echo "Check adapter file and rerun the analysis"
            exit 1
        elif [[ $answer == "y" ]]; then
            echo "Continuing with the analysis"
        fi

    fi

    mkdir -p "$output_directory/${name}_STAR"

    STAR --genomeDir "$reference_genome_dir" --readFilesCommand zcat --readFilesIn "$output_directory/trimmed_files/${name}_r1_paired.fq.gz" "$output_directory/trimmed_files/${name}_r2_paired.fq.gz" --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN "$threads" --outFileNamePrefix "$output_directory/${name}_STAR/${name}"

    mkdir -p "$output_directory/RSEM_counts"

    rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 "$output_directory/${name}_STAR/${name}Aligned.toTranscriptome.out.bam" "$rsem_dir" "$output_directory/RSEM_counts/${name}" >& "$output_directory/RSEM_counts/${name}.log"

done