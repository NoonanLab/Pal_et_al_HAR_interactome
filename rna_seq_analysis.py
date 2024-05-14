import os
import subprocess
import argparse
import gzip
import tarfile

code = []
code.append("#!/bin/bash")
code.append("\n")
'''
ml Trimmomatic
ml STAR
ml RSEM
ml SAMtools
ml FastQC'''

code.append("ml Trimmomatic")
code.append("ml STAR")
code.append("ml RSEM")
code.append("ml SAMtools")
code.append("ml FastQC")
code.append("\n")
# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description='RNA-Seq analysis pipeline using STAR for alignment.')
parser.add_argument('--data_dir', required=True, help='Directory containing your FASTQ files')
parser.add_argument('--output_dir', required=True, help='Output directory for results')
parser.add_argument('--ref_genome_dir', help='Directory of the STAR reference genome indices')
parser.add_argument('--rsem', help='Path to the RSEM reference directory')
parser.add_argument('--ann', help='Path to the annotation GTF file')
parser.add_argument('--th', type=int, default=2, help='Number of threads to use for alignment')
parser.add_argument('--adap', help='Path to the adapter file for Trimmomatic')
# Parse the command-line arguments
args = parser.parse_args()

data_directory = args.data_dir
output_directory = args.output_dir
reference_genome_dir = args.ref_genome_dir
rsem_dir = args.rsem
annotation_file_path = args.ann
threads = args.th
adapter_file = args.adap

# mkdir -p $OUTPUT_DIR
code.append("mkdir -p " + os.path.join(output_directory))

'''mkdir -p $OUTPUT_DIR/fastqc
fastqc -t 5 -o $OUTPUT_DIR/fastqc/ $INPUT_DIR/*.gz
mkdir -p $OUTPUT_DIR/trimmed_files
'''

code.append("mkdir -p " + os.path.join(output_directory, "fastqc"))
code.append("fastqc -t 5 -o " + os.path.join(output_directory, "fastqc") + " " + os.path.join(data_directory, "*.gz"))
code.append("mkdir -p " + os.path.join(output_directory, "trimmed_files"))


fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
true_names_and_file = [(f.split('r1')[0], f) for f in fastq_files if f.endswith('r1.fastq.gz') or f.endswith('r1.fq.gz')]
names = [f[0].rstrip("._") for f in true_names_and_file]

# add a for loop in the code list over the list names using bash script

code.append("for name in " + " ".join(names) + "; do")
code.append("echo starting new analysis for ${name}")

# unzip $OUTPUT_DIR/fastqc/${name}_r1_fastqc.zip
# unzip $OUTPUT_DIR/fastqc/${name}_r2_fastqc.zip
# R1_F=grep -q "Adapter Content\tfail" os.path.join(output_directory, "fastqc", "${name}" + "_r1_fastqc", "fastqc_data.txt")
# R2_F=grep -q "Adapter Content\tfail" os.path.join(output_directory, "fastqc", "${name}" + "_r2_fastqc", "fastqc_data.txt")

# if [ $R1_F -eq 0 ] && [ $R2_F -eq 0 ]
# then

code.append("unzip " + os.path.join(output_directory, "fastqc", "${name}" + "_r1_fastqc.zip"))
code.append("unzip " + os.path.join(output_directory, "fastqc", "${name}" + "_r2_fastqc.zip"))

code.append("R1_F=grep -P -c \"Adapter Content\tfail\" " + os.path.join(output_directory, "fastqc", "${name}" + "_r1_fastqc", "fastqc_data.txt"))
code.append("R2_F=grep -P -c \"Adapter Content\tfail\" " + os.path.join(output_directory, "fastqc", "${name}" + "_r2_fastqc", "fastqc_data.txt"))

code.append("if [ $R1_F -eq 0 ] && [ $R2_F -eq 0 ]")
code.append("then")
code.append("\n")

'''
A=$INPUT_DIR/hcc1395_normal_rep1_r1.fastq.gz  # Replace with the actual file names .fastq.gz or .fq.gz
B=$INPUT_DIR/hcc1395_normal_rep1_r2.fastq.gz 
C=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_paired.fq.gz # STAR uses paired files
D=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_unpaired.fq.gz 
E=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_paired.fq.gz # STAR uses paired files
F=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_unpaired.fq.gz 
#G=ILLUMINACLIP:$INPUT_DIR/TruSeq3-PE.fa:2:30:10  # TruSeq3-PE.fa is the adapter file in form .fa or .fasta
#H=LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $THREADS $A $B $C $D $E $F ILLUMINACLIP:/gpfs/gibbs/pi/noonan/ap2549/RNA-seq_NSC/new_analysis_20230101/trimmomatic_files/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
'''
code.append("A=" + os.path.join(data_directory, "${name}" + "_r1.fastq.gz"))
code.append("B=" + os.path.join(data_directory, "${name}" + "_r2.fastq.gz"))
code.append("C=" + os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_paired.fq.gz"))
code.append("D=" + os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_unpaired.fq.gz"))
code.append("E=" + os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_paired.fq.gz"))
code.append("F=" + os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_unpaired.fq.gz"))
code.append("java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads "+ str(threads) +" $A $B $C $D $E $F ILLUMINACLIP:" + adapter_file + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
code.append("\n")

# else
# mv os.path.join(data_directory, "${name}" + "_r1.fastq.gz") os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_paired.fq.gz")
# mv os.path.join(data_directory, "${name}" + "_r2.fastq.gz") os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_paired.fq.gz")

code.append("else")
code.append("mv " + os.path.join(data_directory, "${name}" + "_r1.fastq.gz") + " " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_paired.fq.gz"))
code.append("mv " + os.path.join(data_directory, "${name}" + "_r2.fastq.gz") + " " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_paired.fq.gz"))
code.append("fi")
code.append("\n")

#mkdir -p $OUTPUT_DIR/trimmed_files/fastqc
#fastqc -t 5 -o $OUTPUT_DIR/trimmed_files/fastqc $OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_paired.fq.gz
#fastqc -t 5 -o $OUTPUT_DIR/trimmed_files/fastqc $OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_paired.fq.gz

code.append("mkdir -p " + os.path.join(output_directory, "trimmed_files", "fastqc"))
code.append("fastqc -t 5 -o " + os.path.join(output_directory, "trimmed_files", "fastqc") + " " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_paired.fq.gz"))
code.append("fastqc -t 5 -o " + os.path.join(output_directory, "trimmed_files", "fastqc") + " " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_paired.fq.gz"))


# unzip $OUTPUT_DIR/trimmed_files/fastqc/${name}_r1_paired_fastqc.zip
# unzip $OUTPUT_DIR/trimmed_files/fastqc/${name}_r2_paired_fastqc.zip

code.append("unzip " + os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r1_paired_fastqc.zip"))
code.append("unzip " + os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r2_paired_fastqc.zip"))

# R1_F=grep -q "Adapter Content\tfail" os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r1_paired_fastqc", "fastqc_data.txt")
# R2_F=grep -q "Adapter Content\tfail" os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r2_paired_fastqc", "fastqc_data.txt")

# if [ $R1_F -eq 0 ] && [ $R2_F -eq 0 ]
# then

code.append("R1_F=grep -P -c \"Adapter Content\tfail\" " + os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r1_paired_fastqc", "fastqc_data.txt"))
code.append("R2_F=grep -P -c \"Adapter Content\tfail\" " + os.path.join(output_directory, "trimmed_files", "fastqc", "${name}" + "_r2_paired_fastqc", "fastqc_data.txt"))

code.append("if [ $R1_F -eq 0 ] && [ $R2_F -eq 0 ]")
code.append("then")
code.append("\n")

# echo "Adapter still present after trimmming. Do you want to continue? (y/n)"
# read answer
# if [ $answer == "n" ]
# then
# echo "Check adapter file and rerun the analysis"
# exit 1
# elif [ $answer == "y" ]
# then
# echo "Continuing with the analysis"
# fi
code.append("echo \"Adapter still present after trimmming. Do you want to continue? (y/n)\"")
code.append("read answer")
code.append("if [ $answer == \"n\" ]")
code.append("then")
code.append("echo \"Check adapter file and rerun the analysis\"")
code.append("exit 1")
code.append("elif [ $answer == \"y\" ]")
code.append("then")
code.append("echo \"Continuing with the analysis\"")
code.append("\n")
code.append("fi")

# fi
code.append("fi")





# mkdir -p $OUTPUT_DIR/STAR_alignment

code.append("mkdir -p " + os.path.join(output_directory, "${name}" + "_STAR"))

# STAR --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn $OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_paired.fq.gz $OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_paired.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN 10 --outFileNamePrefix $OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1

code.append("STAR --genomeDir " + reference_genome_dir + " --readFilesCommand zcat --readFilesIn " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r1_paired.fq.gz") + " " + os.path.join(output_directory, "trimmed_files", "${name}" + "_r2_paired.fq.gz") + " --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN " + str(threads) + " --outFileNamePrefix " + os.path.join(output_directory, "${name}" + "_STAR", "${name}"))

# mkdir -p $OUTPUT_DIR/RSEM_counts

code.append("mkdir -p " + os.path.join(output_directory,"RSEM_counts"))

# rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 $OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1Aligned.toTranscriptome.out.bam $RSEM_DIR $OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1 >& $OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1.log

code.append("rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 " + os.path.join(output_directory, "${name}" + "_STAR", "${name}" + "Aligned.toTranscriptome.out.bam") + " " + rsem_dir + " " + os.path.join(output_directory, "RSEM_counts", "${name}") + " >& " + os.path.join(output_directory, "RSEM_counts", "${name}") + ".log")

code.append("done")

# Write the code to a bash script
with open('gen_bash_trial.sh', 'w') as f:
    f.write('\n'.join(code))

# subprocess.run(["chmod", "+x", "gen_bash.py"], shell=True)
# # Make the bash script executable
# subprocess.run(["chmod", "+x", "gen_bash_trial.sh"], shell=True)

# Run the bash script
subprocess.run(["./gen_bash_trial.sh"], shell=True)