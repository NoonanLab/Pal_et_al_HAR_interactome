def modify_gtf(input_file, output_file, upstream_distance, downstream_distance):
    """
    Modify the start and end positions of gene and transcript features in a GTF file to refine set of human and chimpanzee-specific gene targets.

    Refined species specific targets are those present in one species for the HAR/HGE but there is no interaction of the HAR/HGE-ortholog within 50kb of the target in the other genome

    Parameters:
    input_file (str): Path to the input GTF file.
    output_file (str): Path to the output GTF file.
    upstream_distance (int): Distance to extend upstream of the start position.
    downstream_distance (int): Distance to extend downstream of the end position.
    """
    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Iterate through each line in the input file
        for line in infile:
            if line.startswith('#'):
                # Preserve comment lines as they are
                outfile.write(line)
            else:
                # Split the line into fields
                fields = line.strip().split('\t')
                feature_type = fields[7]  # The feature type is in the 8th column (index 7)

                if feature_type == 'gene' or feature_type == 'transcript':
                    # Modify the start and end positions for genes and transcripts
                    start = max(1, int(fields[1]) - upstream_distance)  # Ensure start is not less than 1
                    end = int(fields[2]) + downstream_distance
                    fields[1] = str(start)
                    fields[2] = str(end)

                # Write the modified line to the output file
                outfile.write('\t'.join(fields) + '\n')

# Example usage:
input_gtf = 'gencode.v43.annotation.bed'
output_gtf = 'gencode.v43.annotation_slop50k.bed'
upstream_distance = 50000  # 50kb
downstream_distance = 50000  # 50kb

# Call the function to modify the GTF file
modify_gtf(input_gtf, output_gtf, upstream_distance, downstream_distance)
