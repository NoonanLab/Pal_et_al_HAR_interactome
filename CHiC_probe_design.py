import argparse
import subprocess
import pybedtools

# Setup argument parser
parser = argparse.ArgumentParser(description="Process BED files for genomic data.")
parser.add_argument("input1", help="HAR coordinates BED file, should contain 'hg19' or 'panTro6' in its name")
parser.add_argument("input2", help="DA fragments in HGEs- coordinates in BED file, should contain 'hg19' or 'panTro6' in its name")
args = parser.parse_args()

# run script first on human (hg19) to get hg38 equivalents, and then run on chimp (panTro6) files
# Determine genome version and set filenames
if 'hg19' in args.input1 and 'hg19' in args.input2:
    chain_file = 'hg19ToHg38.over.chain'
    fasta_file = 'hg38.fa'
    suffix = '_hg38'
elif 'panTro6' in args.input1 and 'panTro6' in args.input2:
    chain_file = 'hg38TopanTro6.over.chain'
    fasta_file = 'panTro6.fa'
    suffix = '_panTro6'
else:
    raise ValueError("Input files must be consistent and contain either 'hg19' or 'panTro6'.")

# Function to run bash commands
def run_command(command):
    subprocess.run(command, shell=True, check=True)

# Liftover elements to genome version accepted by Capsequm2
run_command(f'liftOver {args.input1} {chain_file} {args.input1.split(".")[0]}{suffix}.bed unMapped')
run_command(f'liftOver {args.input2} {chain_file} {args.input2.split(".")[0]}{suffix}.bed unMapped')

# Find DpnII cut sites and generate DpnII-digested genome
run_command(f'restriction-finder.py --fasta {fasta_file} --enzyme DpnII --bed DpnII-sites{suffix}.bed')
run_command(f'sort -k1,1 -k2,2n DpnII-sites{suffix}.bed > DpnII-sites_sort{suffix}.bed')
run_command(f'bedtools complement -i DpnII-sites{suffix}.bed -g {fasta_file.split(".")[0]}.chrom.sizes > DpnII-fragments{suffix}.bed')

# Functions to process BED files using pybedtools
def process_bed_files(file_a, file_b, output_file, d_value, condition):
    a = pybedtools.BedTool(file_a)
    b = pybedtools.BedTool(file_b)
    result = a.intersect(b, u=True).merge(d=d_value).filter(lambda x: condition(int(x.end) - int(x.start)))
    result.saveas(output_file)

process_bed_files(f'DpnII-fragments{suffix}.bed', f'{args.input2.split(".")[0]}{suffix}.bed', f'da_intersect_DpnIIfrags_gt1k{suffix}.bed', 9, lambda d: d > 1000)

# Extend intervals less than 1000 bp
def extend_intervals(input_file, chrom_sizes, extend_by):
    intervals = pybedtools.BedTool(input_file)
    extended = intervals.slop(g=chrom_sizes, b=extend_by)
    intersected = extended.intersect(pybedtools.BedTool(f'DpnII-fragments{suffix}.bed'), u=True)
    intersected.saveas(input_file.replace('lt1k', 'lt1k_extend'))

extend_intervals(f'da_intersect_DpnIIfrags_lt1k{suffix}.bed', f'{fasta_file.split(".")[0]}.chrom.sizes', 109)

# Generate batches for submission
def generate_batches(input_file, batch_size):
    with open(input_file) as file:
        lines = file.readlines()
        for i in range(0, len(lines), batch_size):
            with open(f'submission_batch{int(i/batch_size) + 1}{suffix}.bed', 'w') as batch_file:
                batch_file.writelines(lines[i:i+batch_size])

generate_batches(f'target_intersect_DpnIIfrags_center_forSubmission{suffix}.bed', 1000)

# Collect results after CapSequm2 runs are complete
def collect_results(files, output_file):
    with open(output_file, 'w') as outfile:
        for f in files:
            with open(f) as infile:
                outfile.write(infile.read())
                outfile.write('\n')

collect_results([f'submission_batch{i}_100_GATC_oligonucleotides{suffix}.txt' for i in range(1, 6)], f'submissions_100_GATC_oligonucleotides{suffix}.txt')
