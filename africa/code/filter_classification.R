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
    make_option(c("-e", "--metadata"), type = "character", help = "Input TSV file containing metadata from GenBank."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input FASTA file."),
    make_option(c("-m", "--min_count"), type = "numeric", default = 100, help = "Minimum number of sequences required to retain [default= %default]."),
    make_option(c("-t", "--target"), type = "character", help = "Specify target region (country or continent)."),
    make_option(c("-c", "--classification"), type = "character", help = "Specify dengue lineage granularity to filter."),
    make_option(c("--fasta_outfile"), type = "character", default = NULL, help = "Output path for the filtered FASTA file."),
    make_option(c("--csv_outfile"), type = "character", default = NULL, help = "Output path for the filtered metadata CSV file.")
  )
)

opt <- parse_args(opt_parser)

########################################################################
## Main
########################################################################

# Read in input metadata file
if (!is.null(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
} else {
  cat("Input metadata file is required. Exiting now...")
  quit()
}

# Read in input FASTA file
if (!is.null(opt$fasta)) {
  seqs <- read.fasta(opt$fasta)
} else {
  cat("Input FASTA file is required. Exiting now...")
  quit()
}

# Clean and process metadata
metadata <- metadata %>%
  mutate(
    Country = gsub("_", " ", Country),
    Country = ifelse(Country == 'Micronesia', 'Micronesia, Federated States of', Country),
    Continent = countrycode(Country, origin = "country.name", destination = "continent"),
    Continent = ifelse(is.na(Continent), '_', Continent)
  )

# Filter by target region (country or continent)
if (!is.null(opt$target)) {
  metadata_filtered <- metadata %>%
    filter(Country == opt$target | Continent == opt$target)
}

# Filter genotypes by minimum count
genotype_counts <- metadata_filtered %>%
  group_by(.data[[opt$classification]]) %>%
  summarise(count = n())

allowed_genotypes <- genotype_counts %>%
  filter(count >= opt$min_count) %>%
  pull(.data[[opt$classification]])

metadata <- metadata %>%
  filter(.data[[opt$classification]] %in% allowed_genotypes)

# Filter sequences
seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))
keep <- seq_name %>%
  filter(V1 %in% metadata$Sequence_name)

taxa_split <- data.frame(do.call('rbind', strsplit(as.character(keep$V1), '|', fixed = TRUE)))
taxa_split$name <- keep$V1
species_to_keep <- taxa_split$name

vec_names <- unlist(lapply(strsplit(names(seqs), ";"), function(x) x[length(x)]))
vec_to_keep <- which(vec_names %in% species_to_keep)

# Write filtered FASTA
write.fasta(sequences = seqs[vec_to_keep], names = names(seqs)[vec_to_keep],
            file.out = opt$fasta_outfile)

# Write filtered metadata CSV
write.csv(metadata, file = opt$csv_outfile, row.names = FALSE, quote = FALSE)
