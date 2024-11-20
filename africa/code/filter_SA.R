## load packages
required_packages <- c("dplyr","lubridate","tidyr","optparse","ape",
                       "seqinr")
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
    make_option(c("-e", "--metadata"), type="character", help="Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type="character", help="Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type="character", default="subsampled", help="Base name for output files. Files will be named as '<outfile>.fasta'.")
  )
)

opt = parse_args(opt_parser)

########################################################################
## main
########################################################################


## read in input metadata file
if (!is.null(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
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

########################################################################


South_America <- c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
                   "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname",
                   "Uruguay", "Venezuela")

metadata_sea <- filter(metadata, Country %in% South_America)

#filter sequences

seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

keep <- seq_name %>%
  filter(V1 %in% metadata_sea$Sequence_name)

taxa_split <- data.frame(do.call('rbind',strsplit(as.character(keep$V1),'|',fixed = TRUE)))
taxa_split$name <- keep$V1
species.to.keep <- taxa_split$name

vec.names<-unlist(lapply(strsplit(names(seqs), ";"), function(x)x[length(x)]))

vec.tokeep <-which(vec.names %in%  species.to.keep)

#write files 

write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=paste0(opt$outfile, '.fasta'))

write.csv(metadata_sea,
          file = paste0(opt$outfile, "_infoTbl.tsv"),
          row.names = FALSE, quote=FALSE)

write.csv(metadata_sea,
          file = paste0(opt$outfile, "_infoTbl.csv"),
          row.names = FALSE, quote=FALSE)