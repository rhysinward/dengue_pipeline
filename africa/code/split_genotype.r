# List of required packages
required_packages <- c("optparse", "dplyr", "stringr")
library(Biostrings)
options(repos = c(CRAN = "http://cran.us.r-project.org"))

# Function to check and install packages
suppressMessages(
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
)

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type = "character",
                help = "Input CSV file containing cleaned metadata, with sequence identifiers matching those in the input FASTA file."),
    make_option(c("-f", "--fasta"), type = "character",
                help = "Input FASTA file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type = "character", default = "subsampled",
                help = "Base name for output files. Files will be named as '<outfile>_<genotype>.fasta', '<outfile>_<genotype>_infoTbl.tsv', and '<outfile>_<genotype>_infoTbl.csv'")
  )
)
opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

# Check if required arguments are provided
if (is.null(opt$metadata) | is.null(opt$fasta)) {
  cat("Both --metadata and --fasta arguments must be provided. Exiting now...\n")
  quit(status = 1)
}

# Read in input cleaned metadata file
metadata.df <- read.csv(opt$metadata, stringsAsFactors = FALSE)

# Verify that metadata contains the necessary columns
required_columns <- c("Sequence_name", "Major_Lineage")  # Adjust column names if necessary
missing_columns <- setdiff(required_columns, colnames(metadata.df))
if (length(missing_columns) > 0) {
  cat("The following required columns are missing in the metadata:", paste(missing_columns, collapse = ", "), "\n")
  quit(status = 1)
}

# Read in input FASTA file
fasta_sequences <- readDNAStringSet(opt$fasta)
fasta_ids <- names(fasta_sequences)

# Ensure that sequence IDs in metadata match those in FASTA
missing_ids <- setdiff(metadata.df$Sequence_name, fasta_ids)
if (length(missing_ids) > 0) {
  cat("Warning: The following sequence IDs are present in the metadata but missing in the FASTA file:\n")
  print(missing_ids)
  # Optionally, you can choose to remove these from metadata
  metadata.df <- metadata.df %>% filter(Sequence_name %in% fasta_ids)
}

# Get unique genotypes
genotypes <- unique(metadata.df$Major_Lineage)
cat("Identified genotypes:", paste(genotypes, collapse = ", "), "\n")

# Iterate over each genotype and write genotype-specific files
for (genotype in genotypes) {
  cat("Processing genotype:", genotype, "\n")
  
  # Subset metadata for the current genotype
  metadata_gen <- metadata.df %>% filter(Major_Lineage == genotype)
  
  # Get sequence IDs for the current genotype
  seq_ids_gen <- metadata_gen$Sequence_name
  
  # Subset FASTA sequences
  fasta_gen <- fasta_sequences[seq_ids_gen]
  
  # Define output paths
  fasta_output_path <- paste0(opt$outfile, genotype, ".fasta")
  metadata_tsv_output_path <- paste0(opt$outfile, genotype, "_infoTbl.tsv")
  metadata_csv_output_path <- paste0(opt$outfile, genotype, "_infoTbl.csv")
  
  # Write genotype-specific FASTA file
  writeXStringSet(fasta_gen, filepath = fasta_output_path)
  cat("Written FASTA to:", fasta_output_path, "\n")
  
  # Write genotype-specific metadata files
  write.table(metadata_gen,
              file = metadata_tsv_output_path,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  
  write.csv(metadata_gen,
            file = metadata_csv_output_path,
            row.names = FALSE)
  
  cat(paste0("Processing completed for genotype ", genotype, "\n"))
}

cat("All genotypes have been processed successfully!\n")
