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

names(seqs) <- gsub("_+", "_", names(seqs))
names(seqs) <- gsub("^_+|_+$", "", names(seqs))

metadata.df$Sequence_name <- gsub("(", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub(")", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub(",", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub("*", "_", metadata.df$Sequence_name, fixed = TRUE)
metadata.df$Sequence_name <- gsub("'", "_", metadata.df$Sequence_name, fixed = TRUE)

metadata.df$Sequence_name <- gsub("_+", "_", metadata.df$Sequence_name)
metadata.df$Sequence_name <- gsub("^_+|_+$", "", metadata.df$Sequence_name)

metadata.df <- metadata.df %>% rename(date = Decimal_Date)
metadata.df <- metadata.df %>% rename(name = Sequence_name)

metadata.df$Country <- gsub("Viet_Nam", "Vietnam", metadata.df$Country)

sequences_to_remove <- c(
  "EPI_ISL_14908847|Brazil|Sao_Paulo|Santa_Barbara_D_oeste|Dengue_1|2022-03-28|2022.23561643836",
  "EPI_ISL_14908874|Brazil|Sao_Paulo|Santa_Barbara_D_oeste|Dengue_1|2022-04-19|2022.29589041096",
  "EPI_ISL_19358701|Brazil|Sao_Paulo|Santa_Barbara_D_Oeste|Dengue_1|2024-04-02|2024.13698630137",
  "EPI_ISL_19358701|Brazil|Sao_Paulo|Santa_Barbara_D_Oeste|Dengue_1|2024-02-20|2024.13698630137",
  "EPI_ISL_17699762|Brazil|Sao_Paulo|Santa_Barbara_D_Oeste|Dengue_1|2023-02-15|2023.12328767123",
  "EPI_ISL_17699762|Brazil|Sao_Paulo|Santa_Barbara_D_Oeste|Dengue_1|2023-03-02|2023.12328767123",
  "EPI_ISL_19376613|Brazil|Rio_Grande_do_Sul|Sant_Ana_do_Livramento|Dengue_1|2024-03-10|2024.18630136986",
  "EPI_ISL_14908869|Brazil|Sao_Paulo|Santa_Barbara_D_oeste|Dengue_1|2022-04-18|2022.29315068493"
)

seqs <- seqs[!names(seqs) %in% sequences_to_remove]


Serotype <- unique(metadata.df$Serotype)

write.fasta(sequences=seqs, names=names(seqs),
            file.out=paste0(opt$output_dir_fasta))

write.csv(metadata.df,
          file = paste0(opt$output_dir_csv),
          row.names = FALSE,quote=FALSE)
