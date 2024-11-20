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
option_list <- list(
  make_option(c("-f", "--fasta"), type="character", 
              help="Input fasta file, with sequence identifiers matching those in the metadata file."),
  make_option(c("-w", "--WG_threshold"), type="numeric", default=0.31, 
              help="Set threshold for proportion of unknown bases in whole genome."),
  make_option(c("-e", "--EG_threshold"), type="numeric", default=0.31, 
              help="Set threshold for proportion of unknown bases in E-gene."),
  make_option(c("-k", "--outfile_fasta_eg"), type="character", 
              help="Outfile fasta for E-gene"),
  make_option(c("-d", "--outfile_fasta_wg"), type="character", 
              help="Outfile fasta for wg")
)

opt_parser <- OptionParser(option_list=option_list)

# Parse command-line arguments
opt <- parse_args(opt_parser)
########################################################################
## main
########################################################################

## read in input fasta file
if (!is.null(opt$fasta)) {
  seqs <- read.fasta(opt$fasta)
} else {
  cat("Input fasta file. Exiting now...")
  quit()
}

seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name$V1),'|',fixed = TRUE)))
serotype <- unique(taxa_split$X5)

#get WG gene
seqs_wg <- rapply(seqs,function(x) ifelse(0.31*sum(table(x)) > sum(x == '-'),'good',x), how = "unlist")

vec.tokeep <-which(seqs_wg ==  'good')

length(vec.tokeep)

write.fasta(sequences = seqs[vec.tokeep], names = names(seqs)[vec.tokeep], file = opt$outfile_fasta_wg)

#get EG gene
fasta.df <- read.dna(opt$fasta, format = "fasta", as.matrix = TRUE)

# Map serotypes to their respective column ranges for the gene wanted

gene_ranges <- list(
  "Dengue_1" = c(935, 2419),
  "Dengue_2" = c(937, 2421), 
  "Dengue_3" = c(935, 2413),
  "Dengue_4" = c(939, 2423)  
)

  # Extract the start and end positions for the current serotype
  start_pos <- gene_ranges[[serotype]][1]
  end_pos <- gene_ranges[[serotype]][2]
  
  # Subset the dataframe for the current serotype's gene range
  gene_wanted <- fasta.df[, start_pos:end_pos]

write.dna(gene_wanted, file=opt$outfile_fasta_eg, format="fasta", nbcol=-1, colsep="")
  
#Need to reload the data as different file format needed to remove those with a lot of - 

fasta.df <- read.fasta(opt$outfile_fasta_eg)
seqs <- rapply(fasta.df,function(x) ifelse(opt[["EG_threshold"]]*sum(table(x)) > sum(x == '-'),'good',x), how = "unlist")
vec.tokeep <-which(seqs ==  'good')
write.fasta(sequences = fasta.df[vec.tokeep], names = names(fasta.df)[vec.tokeep], file = opt$outfile_fasta_eg)
fasta.df <- read.fasta(opt$outfile_fasta_eg)
seqs <- rapply(fasta.df,function(x) ifelse(opt[["EG_threshold"]]*sum(table(x)) > sum(x == 'n'),'good',x), how = "unlist")
vec.tokeep <-which(seqs ==  'good')
write.fasta(sequences = fasta.df[vec.tokeep], names = names(fasta.df)[vec.tokeep], file = opt$outfile_fasta_eg)
cat(paste0("Processing for Dengue serotype ", serotype, " completed.\n"))