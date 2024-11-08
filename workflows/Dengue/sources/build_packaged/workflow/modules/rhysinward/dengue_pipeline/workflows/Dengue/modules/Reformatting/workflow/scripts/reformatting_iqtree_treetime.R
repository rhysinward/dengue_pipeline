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
    make_option(c("-o", "--output_dir_fasta"), type="character", default="subsampled", help="Output directory for subsampled sequences and metadata"),
    make_option(c("-a", "--output_dir_csv"), type="character", default="subsampled", help="Output directory for subsampled sequences and metadata")
    
  )
)
opt = parse_args(opt_parser)

##########################################################
## read in input metadata file
##########################################################

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

names(seqs) <- gsub("(", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub(")", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub(",", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub("*", "_", names(seqs), fixed = TRUE)
names(seqs) <- gsub("'", "_", names(seqs), fixed = TRUE)

metadata.df$Sequence_name <- gsub("(", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub(")", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub(",", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub("*", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub("'", "_", metadata.df$Sequence_name, fixed = TRUE)

metadata.df <- metadata.df %>% rename(date = Decimal_Date)
metadata.df <- metadata.df %>% rename(name = Sequence_name)

metadata.df$Country <- gsub("Viet_Nam", "Vietnam", metadata.df$Country)


Serotype <- unique(metadata.df$Serotype)

write.fasta(sequences=seqs, names=names(seqs),
            file.out=paste0(opt$output_dir_fasta))

write.csv(metadata.df,
                    file = paste0(opt$output_dir_csv),
                    row.names = FALSE,quote=FALSE)
