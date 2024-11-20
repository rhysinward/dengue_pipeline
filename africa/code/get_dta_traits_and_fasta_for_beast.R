#Metadata for Treetime

###Code Start###


## load packages
required_packages <- c("dplyr","tidyverse","tidyr","lubridate","dplyr",
                       "tidygeocoder", "fuzzyjoin", "stringr", "seqinr",
                       "ape", "optparse")

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
    make_option(c("-w", "--WG_threshold"), type="numeric", default = 0.29,  help="Set threshold for proportion of unknown bases in whole genome."),
    make_option(c("-e", "--EG_threshold"), type="numeric", default = 0.05,  help="Set threshold for proportion of unknown bases in E-gene."),
    make_option(c("-f", "--fasta"), type="character", help="Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type="character", default="subsampled", help="Base name for output files. Files will be named as '<outfile>.fasta'.")
  )
)

########################################################################
## main
########################################################################

#load in fasta 
fasta <- seqinr :: read.fasta("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_github/results/subsampled_Dengue_1_cleaned.fasta")
seq_name <- as.data.frame(as.matrix(attributes(fasta)$names))
seq_name_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name$V1),'|',fixed = TRUE)))

#load in tree
tree <- read.nexus("/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_github/results/dengue_Dengue_1_timetree.nexus")
tree_name <- data.frame(tree$tip.label)
tree_name_split <- data.frame(do.call('rbind',strsplit(as.character(tree_name$tree.tip.label),'|',fixed = TRUE)))


#filter the fasta file to contain only names in the tree

vec.tokeep <-which(seq_name$V1 %in% tree_name$tree.tip.label)

# Write the sequences to a new FASTA file

write.fasta(sequences=fasta[vec.tokeep], names=names(fasta)[vec.tokeep],
            file.out=paste0('/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_github/results/dta.fasta'))

#metadata for DTA

tree_name_split$X2[tree_name_split$X2 == "Viet_Nam"] <- "Vietnam"

#create a dataframe with tree name and location

tree_name_split$traits <- tree_name$tree.tip.label

#rename x2 coloumn to location

colnames(tree_name_split)[2] <- "location"

metadata <- select(tree_name_split,c("traits","location"))

write.table(metadata, file = "/Users/rhysinward/Documents/grapevne_modules/dengue_pipeline_github/results/dta_metadata.txt", sep = "\t", row.names = FALSE)
