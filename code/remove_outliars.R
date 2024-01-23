if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(Biostrings)

# Define file paths for Dengue serotypes
tsv_files <- c("/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/dengue_1/outliers.tsv",
               "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/dengue_2/outliers.tsv",
               "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/dengue_3/outliers.tsv",
               "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/dengue_4/outliers.tsv")

rooting_metadata_files <- c("/Users/rhysinward/Documents/Dengue_anaysis/reference_genomes/sylvatic_dengue1.csv",
                            "/Users/rhysinward/Documents/Dengue_anaysis/reference_genomes/aligned_root_denv2.csv",
                            "/Users/rhysinward/Documents/Dengue_anaysis/reference_genomes/aligned_root_dengue_3.csv",
                            "/Users/rhysinward/Documents/Dengue_anaysis/reference_genomes/sylvatic_dengue4.csv")

fasta_files <- c("/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/Dengue_1_combined.fasta",
                 "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/Dengue_2_combined.fasta",
                 "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/Dengue_3_combined.fasta",
                 "/Users/rhysinward/Documents/Dengue_anaysis/results/treetime_round1/Dengue_4_combined.fasta")

output_files <- c("/Users/rhysinward/Documents/Dengue_anaysis/results/Dengue_1_cleaned.fasta",
                  "/Users/rhysinward/Documents/Dengue_anaysis/results/Dengue_2_cleaned.fasta",
                  "/Users/rhysinward/Documents/Dengue_anaysis/results/Dengue_3_cleaned.fasta",
                  "/Users/rhysinward/Documents/Dengue_anaysis/results/Dengue_4_cleaned.fasta")

# Process each serotype
for (i in 1:4) {
  # Read the TSV file to get the list of IDs to remove
  ids_to_remove <- read.table(tsv_files[i], header = TRUE, sep = "\t")$X
  
  # Add rooting sequence ID to the list of IDs to remove, if not already present
  rooting_seq_id <- read.csv(rooting_metadata_files[i])$name[1]
  if (!rooting_seq_id %in% ids_to_remove) {
    ids_to_remove <- c(ids_to_remove, rooting_seq_id)
  }
  
  # Read the FASTA file
  fasta_sequences <- readDNAStringSet(fasta_files[i])
  
  # Filter the sequences
  filtered_sequences <- fasta_sequences[!names(fasta_sequences) %in% ids_to_remove]
  
  # Write the filtered sequences to a new FASTA file
  writeXStringSet(filtered_sequences, filepath = output_files[i])
}

