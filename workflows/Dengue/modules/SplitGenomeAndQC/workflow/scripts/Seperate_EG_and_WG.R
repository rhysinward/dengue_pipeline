## load packages
required_packages <- c("dplyr","lubridate","tidyr","optparse","ape",
                       "seqinr")
suppressMessages(
  for (package in required_packages) {
    #if (!require(package, character.only = TRUE)) {
    #  install.packages(package, repos = "http://cran.us.r-project.org")
    #}
    library(package, character.only = TRUE)
  }
)

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(
        c("-w", "--WG_threshold"),
        type="numeric",
        default = 0.29,
        help="Set threshold for proportion of unknown bases in whole genome."
    ),
    make_option(
        c("-e", "--EG_threshold"),
        type="numeric",
        default = 0.05,
        help="Set threshold for proportion of unknown bases in E-gene."
    ),
    make_option(
        c("-i", "--input_dir"),
        type="character",
        default = "results",
        help="Set input directory."
    ),
    make_option(
        c("-o", "--output_dir"),
        type="character",
        default = "results",
        help="Set output directory."
    )
  )
)

opt = parse_args(opt_parser)
########################################################################
## main
########################################################################

## Define the array of serotype FASTA files
fasta_files <- c(
    paste0(opt$input_dir, "/Aligned_Dengue_1/nextalign.aligned.fasta"),
    paste0(opt$input_dir, "/Aligned_Dengue_2/nextalign.aligned.fasta"),
    paste0(opt$input_dir, "/Aligned_Dengue_3/nextalign.aligned.fasta"),
    paste0(opt$input_dir, "/Aligned_Dengue_4/nextalign.aligned.fasta")
)

## Loop over the FASTA files
for (fasta_file in fasta_files) {
  # Read in input FASTA file
  fasta.df <- read.fasta(fasta_file)
  # Extract serotype number for file naming
  serotype <- gsub(paste0(opt$input_dir, "/Aligned_Dengue_([0-9]+)/nextalign.aligned.fasta"), "\\1", fasta_file)
  
#get WG gene
seqs_wg <- rapply(fasta.df,function(x) ifelse(opt[["WG_threshold"]]*sum(table(x)) > sum(x == '-'),'good',x), how = "unlist")

vec.tokeep <-which(seqs_wg ==  'good')

length(vec.tokeep)

outfile_wg <- paste0(opt$output_dir, "/Dengue_", serotype, "_WG.fasta")

write.fasta(sequences = fasta.df[vec.tokeep], names = names(fasta.df)[vec.tokeep], file = outfile_wg)

#get EG gene
fasta.df <- read.dna(fasta_file, format = "fasta", as.matrix = TRUE)

# Map serotypes to their respective column ranges for the gene wanted

gene_ranges <- list(
  "1" = c(935, 2419),
  "2" = c(937, 2421),
  "3" = c(935, 2413),
  "4" = c(939, 2423)
)

# Loop through the serotypes
for (serotype in serotype) {
  # Extract the start and end positions for the current serotype
  start_pos <- gene_ranges[[serotype]][1]
  end_pos <- gene_ranges[[serotype]][2]
  
  # Subset the dataframe for the current serotype's gene range
  gene_wanted <- fasta.df[, start_pos:end_pos]
}
  
outfile_eg <- paste0(opt$output_dir, "/Dengue_", serotype, "_EG.fasta")
write.dna(gene_wanted, file=outfile_eg, format="fasta", nbcol=-1, colsep="")
  
  #Need to reload the data as different file format needed to remove those with a lot of - 

fasta.df <- read.fasta(outfile_eg)
seqs <- rapply(fasta.df,function(x) ifelse(opt[["EG_threshold"]]*sum(table(x)) > sum(x == '-'),'good',x), how = "unlist")
vec.tokeep <-which(seqs ==  'good')
outfile_eg <- paste0(opt$output_dir, "/Dengue_", serotype, "_EG.fasta")
write.fasta(sequences = fasta.df[vec.tokeep], names = names(fasta.df)[vec.tokeep], file = outfile_eg)
fasta.df <- read.fasta(outfile_eg)
seqs <- rapply(fasta.df,function(x) ifelse(opt[["EG_threshold"]]*sum(table(x)) > sum(x == 'n'),'good',x), how = "unlist")
vec.tokeep <-which(seqs ==  'good')
outfile_eg <- paste0(opt$output_dir, "/Dengue_", serotype, "_EG.fasta")
write.fasta(sequences = fasta.df[vec.tokeep], names = names(fasta.df)[vec.tokeep], file = outfile_eg)
cat(paste0("Processing for Dengue serotype ", serotype, " completed.\n"))
}
