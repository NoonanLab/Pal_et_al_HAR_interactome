# Pal et al. 2024 
---

## Scripts accompanying the publication titled "Resolving the 3D Interactome of Human Accelerated Regions in Human and Chimpanzee Neurodevelopment"

### Step 1
`CHiC_probe_design.py` - Design CHi-C probes for HARs and HGEs in human and chimpanzee genomes.

### Step 2
`HiCUP_preprocessing.sh` - Preprocess CHi-C data for human and chimpanzee cell lines using the HiCUP pipelines.

### Step 3
`CHiCAGO_significant_interactions.Rmd` - Calling significant interactions using CHiCAGO for human and chimpanzee datasets.

### Step 4
`contact_profile.py` - Generate contact profiles for individual HARs and HGEs.

### Step 5
`gene_annotation.r` - Annotate genes targets of HARs and HGEs.

### Step 6
`species_specific_target_extension.py` - Modify the start and end positions of gene and transcript features in a GTF file to refine set of human and chimpanzee-specific gene targets.

### Step 7
`permutation_test_intervals.r` - Perform permutation tests to assess the significance of the overlap of HAR/HGE set compared to different background elements and genomic regions.

### Step 8
`permutation_test.r` - Perform permutation tests to assess the significance of the overlap of HAR/HGE gene targets compared to different functional gene sets.

### Step 9
`cut_and_run_preprocessing.sh` - Preprocess CUT&RUN data for human and chimpanzee datasets.

### Step 10
`gff_to_bed.py` - Convert GFF files to BED files for gene and transcript features.

### Step 11
`loop_type_call.py` - Identify activating (H3K27ac-marked) and repressed (H3K27me3-marked) loops in human and chimpazee NSCs and human neurons. Also identify the cases when the repressed loops in NSCs switch to activating loops in neurons.

### Step 12
`differetial_enrichment.r` - Identify differentially enriched regions between human and chimpanzee NSCs for either H3K27ac or H3K27me3 marking.

### Step 13
`rna_seq_analysis.sh` - Perform RNA-seq analysis for human and chimpanzee NSCs and neurons.

### Step 14
`DESeq2_RNA_seq.Rmd` - Perform differential gene expression analysis using DESeq2 between human and chimpanzee NSCs.

### Step 15
`HAR_contact_fetal_atlas_expression_bias.Rmd` - Assess the expression bias of HAR and HGE gene targets in different cell types of the human fetal brain.