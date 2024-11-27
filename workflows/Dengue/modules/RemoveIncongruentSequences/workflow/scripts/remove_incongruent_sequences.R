# scripts/remove_incongruent_sequences.R

## Load packages
required_packages <- c("dplyr", "tidyr",
                       "seqinr", "ape", "optparse")

suppressMessages({
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
})

# Load optparse for command-line argument parsing
library(optparse)

# Define command-line options
option_list = list(
  make_option(c("-f", "--fasta"), type = "character", 
              help = "Input FASTA file."),
  make_option(c("-t", "--tree"), type = "character", 
              help = "Input tree file."),
  make_option(c("-e", "--metadata"), type = "character", 
              help = "Input TSV file containing metadata from GenBank, with sequence identifiers matching those in the input FASTA file."),
  make_option(c("-o", "--outfile"), type = "character", 
              default = "subsampled", 
              help = "Base name for output files.")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list = option_list))

##########################################################
## Main Execution
##########################################################

# Load FASTA
if (!is.null(opt$fasta)) {
  fasta <- read.fasta(opt$fasta)
} else {
  stop("Input FASTA file is missing. Exiting now...")
}

seq_names <- as.data.frame(as.matrix(attributes(fasta)$names), stringsAsFactors = FALSE)

# Load Tree
if (!is.null(opt$tree)) {
  tree <- read.tree(opt$tree)
} else {
  stop("Input tree file is missing. Exiting now...")
}

tree_names <- data.frame(name = tree$tip.label, stringsAsFactors = FALSE)
class(seq_names$V1)
# Filter FASTA based on tree tips
to_keep <- seq_names$V1 %in% tree_names$name
filtered_fasta <- fasta[to_keep]

# Load Metadata
if (!is.null(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
} else {
  stop("Input metadata file is missing. Exiting now...")
}

# Clean Metadata
metadata_cleaned <- filter(metadata, name %in% tree_names$name)

# Write Cleaned Metadata
write.csv(metadata_cleaned,
          file = paste0(opt$outfile, "_infoTbl.csv"),
          row.names = FALSE, quote=FALSE)


# Write Filtered FASTA
write.fasta(sequences = filtered_fasta, 
            names = seq_names$V1[to_keep], 
            file.out = paste0(opt$outfile, ".fasta"))
