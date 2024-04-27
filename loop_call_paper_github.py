#!/usr/bin/env python3

import argparse
import pybedtools
import pandas as pd
from gff2bed import gff_to_bed  # Import function from gff2bed.py

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calling active and repressed loops.")
    parser.add_argument("--chromatin_mark", required=True, help="Consensus BED file for K27me3 or K27ac")
    parser.add_argument("--chrom_sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--gff3", required=True, help="GFF3 annotation file")
    parser.add_argument("--promoters_bed", required=True, help="Path to promoters BED file")
    parser.add_argument("--mapping", required=True, help="Mapping file with HARs/HGEs anf gene targets")
    parser.add_argument("--all_elements", required=True, help="Combined coordinates of all HARs and HGEs")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    return parser.parse_args()

def run_gff2bed(gff3_path):
    gff3_lines = open(gff3_path).readlines()
    bed_lines = gff_to_bed(gff3_lines)
    gff_bed = pybedtools.BedTool('\n'.join(bed_lines), from_string=True)
    return gff_bed

def process_files(args):
    # Processing the chromatin_mark file by adding a 500bp window
    chromatin_mark = pybedtools.BedTool(args.chromatin_mark)
    mark_slopped = chromatin_mark.slop(g=args.chrom_sizes, b=500)
    slopped_bed = mark_slopped.cut(x=[0, 1, 2])

    # Convert GFF3 to BED
    gff_bed = run_gff2bed(args.gff3)

    # Intersect GFF and slopped BED
    intersection = gff_bed.intersect(slopped_bed, wa=True, wb=True)
    intersect_df = intersection.to_dataframe(names=["chr1", "start1", "end1", "chr2", "start2", "end2", "gene_region", "score", "strand", "attributes"])

    # Truncating and filtering out unnecessary columns
    intersect_df['attributes'] = intersect_df['attributes'].apply(lambda x: dict(item.split("=") for item in x.split(";") if "=" in item))
    intersect_df['gene_name'] = intersect_df['attributes'].apply(lambda x: x.get('gene_name', ''))
    intersect_df['gene_type'] = intersect_df['attributes'].apply(lambda x: x.get('gene_type', ''))

    # Filter rows based on feature type
    relevant_features = ["transcript", "five_prime_UTR", "three_prime_UTR"]
    filtered_df = intersect_df[intersect_df['gene_region'].isin(relevant_features)]

    # Ensure unique gene names and sort
    unique_genes_df = filtered_df[['gene_name']].drop_duplicates().sort_values(by='gene_name')

    final_df = intersect_df[(intersect_df['gene_name'] != '') & (intersect_df['gene_type'] != '')]
    final_df = final_df[['gene_name']].drop_duplicates()
    
    # Subset out gene names with intersections to promoter regions (+/- 5kb)
    final_bed = pybedtools.BedTool.from_dataframe(final_df)

    # Load promoters BED file
    promoters_bed = pybedtools.BedTool(args.promoters_bed)

    # Perform intersection with 'u' option
    intersected_bed = final_bed.intersect(promoters_bed, u=True)

    # Convert to DataFrame and rename 
    final_df = intersected_bed.to_dataframe()

    # Intersection process for all elements. 
    # For just finding active/repressed loops I use the regular intersection of any element within the 500bp window of chromatin mark
    # For finding switch from repressed (NSC) to active (neuron), I use a slightly more relaxed criteria of intersection of within 500bp window of the element or the chromatin mark or both
    all_elements = pybedtools.BedTool(args.all_elements)
    all_elements_intersected = all_elements.intersect(slopped_bed, u=True)
    
    # Read mapping data
    df_map = pd.read_table(args.mapping, sep='\t')
    mapping = {row['har']: row['gene'] for index, row in df_map.iterrows()}
    
    # Perform the analysis on intersected data
    all_elements_df = all_elements_intersected.to_dataframe(names=["chr", "start", "end", "name"])
    results = []

    # Find loops
    for index, row in all_elements_df.iterrows():
        har = row['name']
        if har in mapping:
            genes = mapping[har].split(',')
            for gene in genes:
                if gene in final_df['gene_name'].values:
                    results.append((har, gene))

    # Save the results
    with open(f"{args.output_dir}/hNeuron_active_interactions.txt", "w") as file:
        file.write("Number of interactions: " + str(len(results)) + "\n")
        for har, gene in results:
            file.write(f"{har}\t{gene}\n")
    # note, for switch cases, intersect the K27me3_hNSC_loops with K27ac_neuron_loops (can do in bash)

def main():
    args = parse_arguments()
    process_files(args)

if __name__ == "__main__":
    main()
