library(data.table)

# Function to generate promoter regions from GTF file
generate_promoter_regions <- function(gtf_file, promoter_distance, output_dir) {
  promoters_for_bed <- file.path(output_dir, "promoters.for1.bed")
  promoters_rev_bed <- file.path(output_dir, "promoters.rev1.bed")
  promoters_bed <- file.path(output_dir, "promoters_5kb.bed")
  
  system(paste(
    "awk -vFS='\t' -vOFS='\t' '($7 == \"+\"){ print $1, ($4 - 1), $4, $9, $5, $7, $9, $3; }'",
    gtf_file,
    "| bedops --range -", promoter_distance, ":", promoter_distance, "--everything - >",
    promoters_for_bed
  ))

  system(paste(
    "awk -vFS='\t' -vOFS='\t' '($7 == \"-\"){ print $1, ($5), ($5 + 1), $9, $5, $7, $9, $3; }'",
    gtf_file,
    "| bedops --range -", promoter_distance, ":", promoter_distance, "--everything - >",
    promoters_rev_bed
  ))

  system(paste("bedops --everything", promoters_for_bed, promoters_rev_bed, ">", promoters_bed))

  return(promoters_bed)
}

# Function to process the GTF or promoter regions and intersect with targets
process_intersection <- function(targets_bed, intersect_bed, output_dir, result_prefix) {
  intersect_output <- file.path(output_dir, paste0(result_prefix, "_intersect.bed"))
  trunc_output1 <- file.path(output_dir, paste0(result_prefix, "_trunc1.bed"))
  trunc_output2 <- file.path(output_dir, paste0(result_prefix, "_trunc2.bed"))
  trunc_output3 <- file.path(output_dir, paste0(result_prefix, "_trunc3.bed"))
  trunc_output4 <- file.path(output_dir, paste0(result_prefix, "_trunc4.bed"))
  trunc_output5 <- file.path(output_dir, paste0(result_prefix, "_trunc5.bed"))
  trunc_output6 <- file.path(output_dir, paste0(result_prefix, "_trunc6.bed"))
  final_output <- file.path(output_dir, paste0(result_prefix, "_alt_promoters.bed"))
  gene_list_output <- file.path(output_dir, paste0(result_prefix, "_gene_list.txt"))
  
  system(paste(
    "bedtools intersect -wa -wb -a", targets_bed, "-b", intersect_bed, ">",
    intersect_output
  ))

  system(paste("cut -f1,2,3,4,5,6,7,8,11,12,14", intersect_output, ">", trunc_output1))
  system(paste("cut -d';' -f1,5,6,7", trunc_output1, ">", trunc_output2))
  system(paste("cat", trunc_output2, "| tr ';' '\\t' >", trunc_output3))
  system(paste("cut -f1,2,3,4,5,6,7,8,9,12", trunc_output3, ">", trunc_output4))
  system(paste("awk -F'\\t' '$10 ~ /^gene_name=/ { print $0 }'", trunc_output4, ">", trunc_output5))
  system(paste("awk -v OFS='\\t' '{$10=$10;sub(/gene_name=/, \"\", $10); print}'", trunc_output5, ">", trunc_output6))

  system(paste(
    "awk 'BEGIN{FS=OFS=\"\\t\"} {if($9==\"transcript\" || $9==\"five_prime_UTR\" || $9==\"three_prime_UTR\"){print $1,$2,$3,$4,$5,$6,$7,$10}}'",
    trunc_output6, ">", final_output
  ))

  system(paste(
    "awk '{print $8}'", final_output, "| sort | uniq >", gene_list_output
  ))

  return(gene_list_output)
}

# Main function to run the full process
process_gtf_with_gene_list <- function(
  gtf_file,
  targets_bed,
  output_dir = "output",
  promoter_distance = 5000
) {
  # Create output directory if it does not exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Generate promoter regions
  promoters_bed <- generate_promoter_regions(gtf_file, promoter_distance, output_dir)
  
  # Process promoter regions and get gene list
  promoter_gene_list_file <- process_intersection(targets_bed, promoters_bed, output_dir, "promoters")

  # Process the entire GTF file and get gene list
  gtf_gene_list_file <- process_intersection(targets_bed, gtf_file, output_dir, "gtf")
  
  # Read gene lists
  promoter_gene_list <- fread(promoter_gene_list_file, header = FALSE)$V1
  gtf_gene_list <- fread(gtf_gene_list_file, header = FALSE)$V1
  
  # Combine and deduplicate gene lists
  final_gene_list <- unique(c(promoter_gene_list, gtf_gene_list))
  
  # Write final gene list to output
  final_gene_list_output <- file.path(output_dir, "final_gene_list.txt")
  fwrite(data.table(final_gene_list), final_gene_list_output, col.names = FALSE)
  
  return(final_gene_list_output)
}

# Example of how to call the function for processing hNSC interaction data to GENCODEv43 annotation
process_gtf_with_gene_list(
  gtf_file = "gencode_v43_annotation.bed",
  targets_bed = "HAR_HGE_interactions.bed",
  output_dir = "output",
  promoter_distance = 5000
)
