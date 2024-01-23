# Load packages
required_packages <- c("dplyr", "lubridate", "tidyr", "ape", "seqinr")
suppressMessages(
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
)

fasta_files <- c("results/subsampled_Dengue_1.fasta", "results/subsampled_Dengue_2.fasta", 
                 "results/subsampled_Dengue_3.fasta", "results/subsampled_Dengue_4.fasta")

reference_sequences <- c("reference_genomes/sylvatic_dengue1.fasta", "reference_genomes/aligned_root_denv2.fasta", 
                         "reference_genomes/aligned_root_dengue_3.fasta", "reference_genomes/sylvatic_dengue4.fasta")

metadata_reference <- c("reference_genomes/sylvatic_dengue1.csv", "reference_genomes/aligned_root_denv2.csv", 
              "reference_genomes/aligned_root_dengue_3.csv", "reference_genomes/sylvatic_dengue4.csv")

metadata_sample <- c("results/subsampled_Dengue_1_infoTbl.csv", "results/subsampled_Dengue_2_infoTbl.csv", 
                     "results/subsampled_Dengue_3_infoTbl.csv", "results/subsampled_Dengue_4_infoTbl.csv")

gene_ranges <- list(
  "1" = c(935, 2419),
  "2" = c(937, 2421),
  "3" = c(935, 2413),
  "4" = c(939, 2423)
)

for (i in 1:length(fasta_files)) {
  fasta_file <- fasta_files[i]
  reference_file <- reference_sequences[i]
  
  subsampled_df <- read.dna(fasta_file, format = "fasta", as.matrix = TRUE)
  
  reference_df <- read.dna(reference_file,format = "fasta", as.matrix = TRUE)
  
  serotype <- gsub("results/subsampled_Dengue_([0-9]+).fasta", "\\1", fasta_file)
  
  start_pos <- gene_ranges[[serotype]][1]
  end_pos <- gene_ranges[[serotype]][2]
  trimmed_reference <- reference_df[, start_pos:end_pos]
  
  combined_df <- rbind(subsampled_df, trimmed_reference)
  
  outfile_combined <- paste0("results/Dengue_", serotype, "_combined.fasta")
  write.dna(combined_df, file=outfile_combined, format="fasta", nbcol=-1, colsep="")
  
  cat(paste0("Processing and combining fasta for Dengue serotype ", serotype, " completed.\n"))
}

#add metadata reference to sample metadata

for (i in 1:length(metadata_sample)) {
  sample_metadata <- read.csv(metadata_sample[i], stringsAsFactors = FALSE)
  reference_metadata <- read.csv(metadata_reference[i], stringsAsFactors = FALSE)
  
  sample_metadata <- rbind(sample_metadata, reference_metadata)
  sample_metadata$date <- as.Date(sample_metadata$date, format = "%Y-%m-%d")
  
  outfile <- paste0("results/Dengue_", i, "_combined_metadata.csv")
  write.csv(sample_metadata, file = outfile, row.names = FALSE)
  
  cat(paste0("Processing and combining metadata for Dengue serotype ", i, " completed.\n"))
}
