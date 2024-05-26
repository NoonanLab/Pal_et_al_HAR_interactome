# Load necessary libraries
library(Rsubread)
library(DESeq2)
library(apeglm)

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Please provide the hg38 and panTro6 pseudo peak BED files (H3K27ac and H3K27me3) and metadata files as arguments.")
}

# Pseudo peak files are the union of consensus peaks from human and chimpanzee NSCs
# Pseudo peak files have the following format: {mark}_{species}_consensus_peaks.bed

hg38_H3K27ac_bed_file <- args[1] # BED file for H3K27ac pseudo peaks in hg38 coordinates
hg38_H3K27me3_bed_file <- args[2] # BED file for H3K27me3 pseudo peaks in hg38 coordinates
panTro6_H3K27ac_bed_file <- args[3] # BED file for H3K27ac pseudo peaks in panTro6 coordinates
panTro6_H3K27me3_bed_file <- args[4] # BED file for H3K27me3 pseudo peaks in panTro6 coordinates
metadata_H3K27ac_file <- args[5] # Metadata file for H3K27ac samples
metadata_H3K27me3_file <- args[6] # Metadata file for H3K27me3 samples

# Define conditions and replicates
marks <- c("H3K27ac", "H3K27me3")
species <- c("human", "chimpanzee")
reps <- c("rep1", "rep2")

# Function to create SAF file
create_saf_file <- function(bed_file, saf_file) {
  bed_data <- read.table(bed_file, header = FALSE)
  saf_data <- data.frame(
    GeneID = bed_data$V4,
    Chr = bed_data$V1,
    Start = bed_data$V2 + 1,
    End = bed_data$V3,
    Strand = "."
  )
  saf_data
}

# Create SAF files from BED files for each mark and replicate
for (mark in marks) {
  hg38_bed_file <- if (mark == "H3K27ac") hg38_H3K27ac_bed_file else hg38_H3K27me3_bed_file
  panTro6_bed_file <- if (mark == "H3K27ac") panTro6_H3K27ac_bed_file else panTro6_H3K27me3_bed_file
  
  for (sp in species) {
    for (rep in reps) {
      sample <- paste(sp, mark, rep, sep = "_")
      bed_file <- if (sp == "human") hg38_bed_file else panTro6_bed_file
      saf_file <- paste0(sample, "_consensus_peaks.saf")
      
      # Check if the BED file exists
      if (!file.exists(bed_file)) {
        stop(paste("BED file not found:", bed_file))
      }
      
      # Create SAF file
      saf_data <- create_saf_file(bed_file, saf_file)
      assign(saf_file, saf_data)
    }
  }
}

# Define directory for mapped files
mapped_dir <- "mapped/"

# Differential enrichment analysis for each mark
for (mark in marks) {
  output_files <- list()
  
  for (sp in species) {
    for (rep in reps) {
      sample <- paste(sp, mark, rep, sep = "_")
      bam_file <- file.path(mapped_dir, paste0(sample, ".bam"))
      saf_file <- paste0(sample, "_consensus_peaks.saf")
      
      # Check if the BAM file exists
      if (!file.exists(bam_file)) {
        stop(paste("BAM file not found:", bam_file))
      }
      
      # Check if the SAF file exists
      if (!exists(saf_file)) {
        stop(paste("SAF file not found:", saf_file))
      }
      
      saf_data <- get(saf_file)
      
      # Run featureCounts
      fc <- featureCounts(files = bam_file, annot.ext = saf_data, isPairedEnd = TRUE, 
                          GTF.featureType = "SAF", useMetaFeatures = FALSE, 
                          nthreads = 4)
      
      output_files[[paste(sp, rep, sep = "_")]] <- fc$counts
    }
  }
  
  # Create counts matrix
  counts_matrix <- do.call(cbind, output_files)
  colnames(counts_matrix) <- paste(rep(species, each = length(reps)), mark, reps, sep = "_")
  
  # Differential analysis using DESeq2
  countdata <- round(counts_matrix)
  metadata_file <- if (mark == "H3K27ac") metadata_H3K27ac_file else metadata_H3K27me3_file
  
  metadata <- read.delim(metadata_file, row.names = 1)
  metadata$sampleid <- row.names(metadata)
  metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]
  
  ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                   colData = metadata,
                                   design = ~ Group)
  
  ddsMat <- DESeq(ddsMat)
  
  res <- results(ddsMat, name = "Group_Human_vs_Chimp", pAdjustMethod = "BH", alpha = 0.01)
  
  res_apeglm <- lfcShrink(ddsMat, coef = "Group_Human_vs_Chimp", res = res, type = "apeglm")
  output_res <- as.data.frame(res_apeglm)
  resordered <- data.frame(output_res[order(output_res$padj, na.last = NA), ])
  de_peaks <- resordered[resordered$padj < 0.01 & abs(resordered$log2FoldChange) >= 1, ]
  
  # Output the differentially enriched peaks
  print(de_peaks)
}
