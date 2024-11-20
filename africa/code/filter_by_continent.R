# Load required packages
required_packages <- c("dplyr", "lubridate", "tidyr", "optparse", "ape", "seqinr", "countrycode")
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
    make_option(c("-e", "--metadata"), type = "character", help = "Input TSV file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input FASTA file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type = "character", default = "subsampled", help = "Base name for output files. Files will be named as '<outfile>.fasta'."),
    make_option(c("-c", "--continent"), type = "character", help = "Specify the continent to filter data.")
  )
)

opt = parse_args(opt_parser)

# Ensure inputs are provided
if (is.null(opt$metadata) || is.null(opt$fasta) || is.null(opt$continent)) {
  stop("Metadata file, FASTA file, and continent parameter are required. Exiting now...")
}

# Read in input metadata and FASTA files
metadata <- read.csv(opt$metadata)
seqs <- read.fasta(opt$fasta)

# Remove underscores and convert country names to continent
metadata <- metadata %>%
  mutate(
    Country = gsub("_", " ", Country),
    Continent = countrycode(Country, origin = "country.name", destination = "continent")
  )

# Filter metadata based on the specified continent
metadata_filtered <- filter(metadata, Continent == opt$continent)

# Filter sequences
seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))
keep <- seq_name %>%
  filter(V1 %in% metadata_filtered$Sequence_name)

taxa_split <- data.frame(do.call('rbind', strsplit(as.character(keep$V1), '|', fixed = TRUE)))
taxa_split$name <- keep$V1
species_to_keep <- taxa_split$name

vec_names <- unlist(lapply(strsplit(names(seqs), ";"), function(x) x[length(x)]))
vec_to_keep <- which(vec_names %in% species_to_keep)

# Write output files
write.fasta(sequences = seqs[vec_to_keep], names = names(seqs)[vec_to_keep],
            file.out = paste0(opt$outfile, ".fasta"))

write.csv(metadata_filtered,
          file = paste0(opt$outfile, "_infoTbl.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(metadata_filtered,
          file = paste0(opt$outfile, "_infoTbl.tsv"),
          row.names = FALSE, quote = FALSE)
