# List of required packages
required_packages <- c("dplyr", "lubridate", "tidyr", "optparse", "ape", "seqinr", "readr")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-f", "--fasta"),
      type = "character",
      help = "Input fasta file, with sequence identifiers matching those in the metadata file."
    ),
    make_option(c("-w", "--WG_threshold"),
      type = "numeric", default = 0.31,
      help = "Set threshold for proportion of unknown bases in whole genome."
    ),
    make_option(c("-e", "--EG_threshold"),
      type = "numeric", default = 0.31,
      help = "Set threshold for proportion of unknown bases in E-gene."
    ),
    make_option(c("-k", "--outfile_fasta_eg"),
      type = "character",
      help = "Outfile fasta for E-gene"
    ),
    make_option(c("-d", "--outfile_fasta_wg"),
      type = "character",
      help = "Outfile fasta for wg"
    )
  )
)

opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

info_msg("Running with `WG_threshold`=", opt$WG_threshold, " and `EG_threshold`=", opt$EG_threshold)

# Read in input fasta file
fasta_df <- safe_read_file_param(opt$fasta, read.fasta, required = TRUE)

# Extract serotype number for file naming
seq_name <- as.data.frame(as.matrix(attributes(fasta_df)$names))
taxa_split <- data.frame(do.call("rbind", strsplit(as.character(seq_name$V1), "|", fixed = TRUE)))
serotype <- unique(taxa_split$X5)


# Get WG gene
seqs_wg <- rapply(fasta_df, function(x) ifelse(opt$WG_threshold * sum(table(x)) > sum(x == "-"), "good", x), how = "unlist")
vec.tokeep <- which(seqs_wg == "good")
info_msg("Number of good WG sequences for `", serotype, "`: ", length(vec.tokeep))
write.fasta(sequences = fasta_df[vec.tokeep], names = names(fasta_df)[vec.tokeep], file = opt$outfile_fasta_wg)

# Get EG gene
fasta_df <- safe_read_file_param(opt$fasta, read.dna, format = "fasta", as.matrix = TRUE, required = TRUE)

# Map serotypes to their respective column ranges for the gene wanted
eg_cut_ranges <- list(
  "Dengue_1" = c(935, 2419),
  "Dengue_2" = c(937, 2421),
  "Dengue_3" = c(935, 2413),
  "Dengue_4" = c(939, 2423)
)

# Extract the start and end positions for the current serotype
start_pos <- eg_cut_ranges[[serotype]][1]
end_pos <- eg_cut_ranges[[serotype]][2]

# Subset the dataframe for the current serotype's gene range
gene_wanted <- fasta_df[, start_pos:end_pos]

write.dna(gene_wanted, file = opt$outfile_fasta_eg, format = "fasta", nbcol = -1, colsep = "")

# Need to reload the data as different file format needed to remove those with a lot of -

fasta_df <- read.fasta(opt$outfile_fasta_eg)
seqs <- rapply(fasta_df, function(x) ifelse(opt$EG_threshold * sum(table(x)) > sum(x == "-"), "good", x), how = "unlist")
vec.tokeep <- which(seqs == "good")
write.fasta(sequences = fasta_df[vec.tokeep], names = names(fasta_df)[vec.tokeep], file = opt$outfile_fasta_eg)

fasta_df <- read.fasta(opt$outfile_fasta_eg)
seqs <- rapply(fasta_df, function(x) ifelse(opt$EG_threshold * sum(table(x)) > sum(x == "n"), "good", x), how = "unlist")
vec.tokeep <- which(seqs == "good")
write.fasta(sequences = fasta_df[vec.tokeep], names = names(fasta_df)[vec.tokeep], file = opt$outfile_fasta_eg)

info_msg("Processing for Dengue serotype `", serotype, "` completed")
