# List of required packages
required_packages <- c(
  "optparse", "dplyr", "lubridate", "tidyr",
  "readr", "ape", "seqinr", "countrycode",
  "ggplot2", "purrr", "zoo", "rlang", "wrswoR", "readr"
)

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type = "character", help = "Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--output_dir_fasta"), type = "character", default = "subsampled", help = "Output directory for subsampled sequences and metadata"),
    make_option(c("-a", "--output_dir_csv"), type = "character", default = "subsampled", help = "Output directory for subsampled sequences and metadata")
  )
)
opt <- parse_args(opt_parser)

##########################################################
## read in input metadata file
##########################################################

metadata_df <- safe_read_file_param(opt$metadata, read_csv, show_col_types = FALSE, required = TRUE)

## read in input fasta file
seqs <- safe_read_file_param(opt$fasta, read.fasta, required = TRUE)

names(seqs) <- gsub("(", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub(")", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub(",", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub("*", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub("'", "_", names(seqs), fixed = TRUE)

metadata_df$Sequence_name <- gsub("(", "_", metadata_df$Sequence_name, fixed = TRUE)
metadata_df$Sequence_name <- gsub(")", "_", metadata_df$Sequence_name, fixed = TRUE)
metadata_df$Sequence_name <- gsub(",", "_", metadata_df$Sequence_name, fixed = TRUE)
metadata_df$Sequence_name <- gsub("*", "_", metadata_df$Sequence_name, fixed = TRUE)
metadata_df$Sequence_name <- gsub("'", "_", metadata_df$Sequence_name, fixed = TRUE)

metadata_df <- metadata_df %>% rename(date = Decimal_Date)
metadata_df <- metadata_df %>% rename(name = Sequence_name)

metadata_df$Country <- gsub("Viet_Nam", "Vietnam", metadata_df$Country)


Serotype <- unique(metadata_df$Serotype)

write.fasta(sequences=seqs, names=names(seqs),
            file.out=paste0(opt$output_dir_fasta))

write.csv(metadata_df,
                    file = paste0(opt$output_dir_csv),
                    row.names = FALSE,quote=FALSE)
