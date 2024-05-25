# Does permutation tests to determine if the number of elements (HARs, HGEs, VISTA elements) that overlap with desired regions (TADs with hsSVs, introns, functional annotation etc) is significantly different from what would be expected by chance

library(data.table) # for fread and fwrite

# Define the function
shuffle_and_intersect <- function(
  input_bed, 
  incl_bed,
  genome_sizes,
  intersect_bed,
  output_file,
  iterations = 20000 # Number of trials
) {
  # Create output directory if it does not exist
  output_dir <- "intermediateFiles"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize result storage
  results <- data.table(withCGI = integer(), denominator = integer())

  # Loop for specified iterations
  for (i in 1:iterations) {
    # Shuffle BED file
    shuffle_cmd <- paste(
      "bedtools shuffle -chrom -incl", incl_bed, 
      "-i", input_bed, 
      "-g", genome_sizes, 
      ">", file.path(output_dir, "hg38_HAR_hsSV_shuffled.bed")
    )
    system(shuffle_cmd)
    
    # Count total number of lines (denominator)
    denominator <- as.integer(system(paste("wc -l < ", file.path(output_dir, "hg38_HAR_hsSV_shuffled.bed")), intern = TRUE))
    
    # Count lines with CGI intersection
    withCGI <- as.integer(system(paste(
      "bedtools intersect -wa -u -a", file.path(output_dir, "hg38_HAR_hsSV_shuffled.bed"), 
      "-b", intersect_bed, "| wc -l"
    ), intern = TRUE))
    
    # Append results to the data table
    results <- rbind(results, data.table(withCGI = withCGI, denominator = denominator))
    
    # Print progress
    if (i %% 100 == 0) {
      cat("Iteration:", i, "\n")
    }
  }

  # Write results to file
  fwrite(results, output_file, sep = "\t", col.names = FALSE)
}

# Example of how to call the function with an example of testing if HARs are enriched in TADs (called in GZ of human brain) containing hsSVs
shuffle_and_intersect(
  input_bed = "HAR_list_hg38.bed",
  incl_bed = "Phastcons_background.bed",
  genome_sizes = "hg38.chrom.sizes",
  intersect_bed = "TADs_GZ_hsSVs_all.bed",
  output_file = "Enrichment_results_HARs_in TAD_with_hsSV.txt",
  iterations = 20000
)