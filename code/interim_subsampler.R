## load packages
required_packages <- c("optparse", "dplyr","lubridate","tidyr",
                       "readr","ape","seqinr","countrycode")
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
    make_option(c("-c", "--country"), type="character", help="Country(ies) of Interest, the first country specified will determine local region"),
    make_option(c("-n", "--number_sequences"), type="character", help="Number of desired sequences (defult is 1:1 ratio with number of sequences from desired country(ies))")
  )
)

opt = parse_args(opt_parser)

##########################################################
## main
##########################################################

##########################################################
# Step 1: read in data
##########################################################
## read in input metadata file
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

##########################################################
# Step 2: processing of metadata
##########################################################

#do metadata and fasta file match? If not remove any in the metadata that do not match fasta file 

# Extract identifiers from FASTA file
seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

# Filter metadata to only include rows with identifiers that match the FASTA file
length <- nrow(metadata.df)
metadata.df <- metadata.df %>%
  filter(Sequence_name %in% seq_name$V1)

metadata.df <- metadata.df %>%
  mutate(Year = year(ymd(Date)))

# Check if any rows were removed
if (nrow(metadata.df) < length) {
  cat("Some identifiers in the metadata did not match the FASTA file. They have been removed.\n")
}

#assign region (SA, NA, Africa, Asia, Europe, Oceania)
metadata.df$Country <- ifelse(metadata.df$Country == "Micronesia", "Micronesia (Federated States of)", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "Borneo", "Indonesia", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "South_Korea", "Korea (Republic of)", metadata.df$Country)
metadata.df$Country <- ifelse(metadata.df$Country == "Republic_of_the_Congo", "Congo (Republic of)", metadata.df$Country)

remove_countries <- c("Pacific_Ocean","Saint_Martin")

metadata.df <- metadata.df %>%
  filter(!Country %in% remove_countries)

metadata.df$Continent <- countrycode(metadata.df$Country, origin = "country.name", destination = "continent")

#select all metadata from country of interest is ok if there are no metadata but number of required sequences needs to be specified
countries_of_interest <-   opt$country

countries_of_interest_first <- strsplit(countries_of_interest, ",")[[1]]

metadata_countries_of_interest <- metadata.df %>%
  filter(Country %in% countries_of_interest)

##############################################################################
# Step 1: Identify the region of the first country
##############################################################################

first_country <- countries_of_interest_first[1]
first_country_region <- countrycode(first_country, origin = "country.name", destination = "continent")

non_selected_country_metadata <- metadata.df %>%
  filter(!Country %in% countries_of_interest)

in_region_countries <- unique(non_selected_country_metadata$Country[non_selected_country_metadata$Continent == first_country_region])

out_of_region_countries <- unique(non_selected_country_metadata$Country[non_selected_country_metadata$Continent != first_country_region])

##############################################################################
# Step 2: Sample from region of and not of interest
##############################################################################

in_region_countries_metdata <- non_selected_country_metadata %>%
     dplyr :: filter(Country %in% in_region_countries) %>%
slice_sample (prop=1, replace = FALSE)

in_region_countries_metdata <- rbind(in_region_countries_metdata,metadata_countries_of_interest)
  
out_of_region_countries_metdata <- non_selected_country_metadata %>%
  dplyr :: filter(Country %in% out_of_region_countries) %>%
  slice_sample (n = round((0.4*nrow(in_region_countries_metdata))), replace = FALSE)

##############################################################################
# Step 3: Combine all sampled sequences into a single data frame
##############################################################################
final_sample <- rbind(out_of_region_countries_metdata,in_region_countries_metdata)
##############################################################################
# Step 4: Match sample to sequences and produce sampled metadata + sequences
##############################################################################
serotype <- unique(metadata.df$Serotype)

colnames(final_sample) <- c('GenBank_ID', "Country", "State",
                            "City", "Serotype", "Sequence_Type", "date","Decimal_Date","name","Year","Continent")

# final_sample <- select(final_sample,c(9,7,1,2,3,4,5,6,8,10,11))
# 
# write.table(final_sample,
#             file = paste0("results/subsampled_",serotype,'_infoTbl.txt'),
#             sep = "\t", 
#             col.names = TRUE, 
#             row.names = FALSE, 
#             quote = FALSE)
# 
# write.csv(final_sample,
#           file = paste0("results/subsampled_",serotype,'_infoTbl.csv'),
#           row.names = FALSE)

seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

keep <- seq_name %>%
  filter(V1 %in% final_sample$name)

taxa_split <- data.frame(do.call('rbind',strsplit(as.character(keep$V1),'|',fixed = TRUE)))
taxa_split$name <- keep$V1
species.to.keep <- taxa_split$name

vec.names<-unlist(lapply(strsplit(names(seqs), ";"), function(x)x[length(x)]))

vec.tokeep <-which(vec.names %in%  species.to.keep)

# Write the sequences to a new FASTA file


seq_name_kept <- as.data.frame(as.matrix(attributes(seqs[vec.tokeep])$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name_kept$V1),'|',fixed = TRUE)))
taxa_split$name <- seq_name_kept$V1
colnames(taxa_split) <- c('GenBank_ID', "Country", "State",
                          "City", "Serotype", "Sequence_Type", "date","Decimal_Date","name")

taxa_split <- select(taxa_split,c(9,7,1,2,3,4,5,6,8))
write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=paste0('results/subsampled_', serotype, '.fasta'))


write.csv(taxa_split,
            file = paste0("results/subsampled_",serotype,'_infoTbl.tsv'),
          row.names = FALSE, quote=FALSE)

write.csv(taxa_split,
          file = paste0("results/subsampled_",serotype,'_infoTbl.csv'),
          row.names = FALSE,quote=FALSE)
