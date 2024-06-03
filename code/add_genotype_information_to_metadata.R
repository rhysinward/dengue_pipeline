# List of required packages
required_packages <- c("optparse", "dplyr", "lubridate", "tidyr",
                       "readr", "ape","seqinr", "stringr")

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
    make_option(c("-m", "--metadata"), type="character", help="Input csv file containing cleaned metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-b", "--genotype"), type="character", help="Input csv file containing typing data, with sequence identifiers matching those in the input metadata file."),
    make_option(c("-l", "--outfile_csv"), type="character", help="Outfile csv")
    
  )
)


opt = parse_args(opt_parser)

########################################################################
## main
########################################################################


## read in typing file
if (!is.null(opt$genotype)) {
  typing.df <- read.csv(opt$genotype, sep = ";")
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

#select genotype column 

typing.df <- dplyr:: select(typing.df,c("seqName","clade"))

#rename to genotype

typing.df <- rename(typing.df,genotype = clade)

# read in metadata file

if (!is.null(opt$metadata)) {
  metadata.df <- read.csv(opt$metadata)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

metadata.df <- left_join(metadata.df, typing.df, by = c("Sequence_name" = "seqName"))

#give any empty clades a value of unassigned

metadata.df$genotype <- sub("^$", "unassigned", metadata.df$genotype)

#write out metadata file

write.csv(metadata.df,
          file = opt$outfile_csv,
          row.names = FALSE)
cat(paste0("Adding serotype infomation complete","\n"))

