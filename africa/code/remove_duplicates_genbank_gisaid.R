# Load required packages
required_packages <- c("optparse", "dplyr", "tidyr",
                       "readr", "ape", "seqinr", "stringr")

# Function to check and load packages
suppressMessages({
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
})

# Set up command-line options
option_list <- list(
  make_option(c("--metadata_gisaid"), type="character", help="Input TSV file containing metadata from GISAID, with sequence identifiers matching those in the input FASTA file."),
  make_option(c("--metadata_genbank"), type="character", help="Input TSV file containing metadata from GenBank, with sequence identifiers matching those in the input FASTA file."),
  make_option(c("--fasta_gisaid"), type="character", help="Input FASTA file from GISAID, with sequence identifiers matching those in the metadata file."),
  make_option(c("--fasta_genbank"), type="character", help="Input FASTA file from GenBank, with sequence identifiers matching those in the metadata file."),
  make_option(c("--outfile_fasta"), type="character", help="Output FASTA file."),
  make_option(c("--outfile_csv"), type="character", help="Output CSV file.")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read FASTA files using read.fasta from seqinr
if (!is.null(opt$fasta_genbank)) {
  genbank_seqs <- read.fasta(file = opt$fasta_genbank, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
} else {
  stop("Input GenBank FASTA file not provided. Exiting now...")
}

if (!is.null(opt$fasta_gisaid)) {
  gisaid_seqs <- read.fasta(file = opt$fasta_gisaid, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
} else {
  stop("Input GISAID FASTA file not provided. Exiting now...")
}

# Read metadata
if (!is.null(opt$metadata_genbank)) {
  genbank_metadata <- read.csv(opt$metadata_genbank)
} else {
  stop("Input GenBank metadata file not provided. Exiting now...")
}

if (!is.null(opt$metadata_gisaid)) {
  gisaid_metadata <- read.csv(opt$metadata_gisaid)
} else {
  stop("Input GISAID metadata file not provided. Exiting now...")
}

# Create data frames from sequences
genbank_seq_df <- data.frame(
  Sequence_name = names(genbank_seqs),
  Sequence = unlist(genbank_seqs),
  stringsAsFactors = FALSE
)

gisaid_seq_df <- data.frame(
  Sequence_name = names(gisaid_seqs),
  Sequence = unlist(gisaid_seqs),
  stringsAsFactors = FALSE
)

# Merge sequences with metadata
genbank_data <- genbank_seq_df %>%
  left_join(genbank_metadata, by = c("Sequence_name" = "Sequence_name"))

gisaid_data <- gisaid_seq_df %>%
  left_join(gisaid_metadata, by = c("Sequence_name" = "Sequence_name"))

# Combine datasets with source identifier
combined_data <- bind_rows(
  genbank_data %>% mutate(Source = "GenBank"),
  gisaid_data %>% mutate(Source = "GISAID")
)

# Create duplicate key
combined_data <- combined_data %>%
  mutate(
    Duplicate_Key = paste(Country, Date, Sequence, sep = "_")
  )

# Assign source preference
combined_data <- combined_data %>%
  mutate(
    Source_Preference = ifelse(Source == "GISAID", 1, 2)
  ) %>%
  arrange(Duplicate_Key, Source_Preference)

# Remove duplicates across sources
filtered_data <- combined_data %>%
  group_by(Duplicate_Key) %>%
  filter(Source_Preference == min(Source_Preference)) %>%
  ungroup()

# Prepare sequences for output
unique_seqs <- filtered_data %>%
  select(Sequence_name, Sequence) %>%
  distinct()

# Write sequences to FASTA file using write.fasta
write.fasta(sequences = as.list(unique_seqs$Sequence),
            names = unique_seqs$Sequence_name,
            file.out = opt$outfile_fasta)

# Save merged metadata
# Remove unnecessary columns
filtered_data <- filtered_data %>% select(-Duplicate_Key, -Source_Preference, -Sequence, -Source)
#change name of column

write.csv(filtered_data,
          file = opt$outfile_csv,
          row.names = FALSE)
