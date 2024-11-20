## load packages
required_packages <- c("optparse", "dplyr","lubridate","tidyr",
                       "readr","ape","seqinr","countrycode",
                       "ggplot2","purrr","zoo","rlang","wrswoR")
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
    make_option(c("-m", "--metadata"), type="character", help="Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type="character", help="Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-c", "--country"), type="character", help="Select country of interest"),
    make_option(c("-o", "--output_dir_fasta"), type="character", default="subsampled", help="Output directory for subsampled sequences and metadata"),
    make_option(c("-a", "--output_dir_csv"), type="character", default="subsampled", help="Output directory for subsampled sequences and metadata")
  )
)
opt = parse_args(opt_parser)

##########################################################
## read in input metadata file
##########################################################
setwd("~/rhys/Brazil_phylo")
if (!is.null(opt$metadata)) {
  metadata.df <- read.csv(opt$metadata)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

## read in input fasta file
if (!is.null(opt$fasta)) {
  seqs <- read.fasta(opt$fasta)
} else {
  cat("Input fasta file. Exiting now...")
  quit()
}

#select sequences with state information from your country of interest and put replace the country with the state

#if country = brazil and state = NA remove 

metadata_cleaned <- metadata.df %>%
  filter(!(Country == opt$country & (is.na(State) | State == "NA"))) %>%
  mutate(Country = ifelse(Country == opt$country, State, Country))

sequence_ids <- metadata_cleaned$name

matching_seq_ids <- sequence_ids %in% names(seqs)
if (any(!matching_seq_ids)) {
  warning(paste("The following GenBank_IDs have no matching sequences in the FASTA file:",
                paste(sequence_ids[!matching_seq_ids], collapse = ", ")))
}

filtered_seqs <- seqs[names(seqs) %in% sequence_ids]

cat("Total sequences in FASTA:", length(seqs), "\n")
cat("Sequences after filtering:", length(filtered_seqs), "\n")

Serotype <- unique(metadata.df$Serotype)

write.fasta(sequences=seqs, names=names(seqs),
            file.out=paste0(opt$output_dir_fasta))

write.csv(metadata.df,
          file = paste0(opt$output_dir_csv),
          row.names = FALSE,quote=FALSE)
