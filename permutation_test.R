set.seed(1)

library(dplyr)

# Get a list of all files in the directory
file_list <- list.files(path = "~/Downloads/Enrichment_psychencode/Macaque/", full.names = TRUE)

# Filter the list to include only files starting with "HS"
filtered_files <- file_list[grep("^MM", basename(file_list))]

# Initialize an empty list to store the results
results <- list()

# Iterate over each file
for (file in filtered_files) {
  # Read the file using read.csv and assign it to the module3 variable
  module3 <- read.csv(file, header = FALSE, sep = "\t")
  
  # Rest of the code block
  contacts <- read.csv("~/Downloads/Enrichment_psychencode/HAR_contacts1.txt", header = FALSE)
  targets <- as.vector(contacts$V1)
  module3_genes <- as.vector(module3$V1)
  all_genes <- read.csv("~/Downloads/Enrichment_psychencode/Background_genes.txt", header = FALSE, sep = "\t")
  all_gene_list <- as.vector(all_genes$V1)
  obs <- length(intersect(targets, module3_genes))
  
  rem_genes <- length(all_gene_list) - length(targets)
  temp_results <- vector('list', 20000)
  for (i in 1:20000) {
    x <- split(sample(all_gene_list), rep(1:2, c(length(targets), rem_genes)))
    temp_results[[i]] <- length(intersect(x[[1]], module3_genes))
  }
  
  df <- data.frame(difs = unlist(temp_results))
  
  p <- sum(unlist(temp_results) >= obs) / 20000
  
  enrichment <- obs / mean(df$difs)
  
  # Store the results in the list
  results[[file]] <- data.frame(p = p, module3_length = length(module3_genes), enrichment=enrichment)
  print("Done")
}

