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
  metadata.df <- read.csv("results/Dengue_1_infoTbl.csv")
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

metadata.df <- metadata.df %>%
  dplyr :: mutate(RowID = row_number())


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

##############################################################################
# Step 2: Input required number of sequences and calculate Sequences per Year
##############################################################################

#if opt$number_sequences is null (I.e. not inputed then automatically select 1:1 ratio) if not subtract 
#required number of sequences from number of sequences from countries of interest
if (is.null(3000)) {
  n_remainder <- as.numeric(nrow(metadata_countries_of_interest) * 2)
} else {
  n_remainder <- as.numeric(3000) - nrow(metadata_countries_of_interest) # Total number of sequences desired
}

if (n_remainder < 0) {
  stop("Error: The number of sequences specified (", opt$number_sequences, 
       ") is less than the number of sequences available in the country of interest (", 
       nrow(metadata_countries_of_interest), ").")
}

num_years <- length(unique(metadata.df$Year)) # Number of years

# In-region sequences per year
in_region_sequences <- round(0.8 * n_remainder)
# Out-of-region sequences per year
out_region_sequences <- round(0.2 * n_remainder)

##############################################################################
# Step 3: Sampling Sequences per Country within a Region
##############################################################################

non_selected_country_metadata <- metadata.df %>%
  filter(!Country %in% countries_of_interest)
  
in_region_countries <- unique(non_selected_country_metadata$Country[non_selected_country_metadata$Continent == first_country_region])

# Initialize variables
finished_sampling <- FALSE

# Create a function to sample one sequence for each entry in countries_vector
library(dplyr)

sample_sequences_in <- function(metadata, countries_vector,end_year, start_year, total_required) {
  # Initialize an empty dataframe to store the sampled sequences
  sampled_sequences <- data.frame()
  
  num_sequences_per_year <- in_region_sequences / (end_year - start_year + 1)  
  
  countries_with_data <- unique(metadata$Country[metadata$Year == year])
  countries_with_data <- intersect(countries_with_data, shuffled_countries)  
  countries_vector <- rep(countries_with_data, length.out = 500)
  
  for (year in start_year:end_year) {
    # Initialize a count for the number of sequences sampled for the current year
    sequences_sampled_this_year <- 0
    
    # Iterate over the countries vector
    for (country in countries_vector) {
      # Stop the inner loop if the required number of sequences for the year is reached
      if (sequences_sampled_this_year >= num_sequences_per_year) {
        break
      }
  
    # Check if the country has some data for the specified year
    if (country %in% countries_vector) {
      # Filter the metadata for the current country and year
      sequences <- metadata %>%
        filter(Country == country, Year == year) 
      
      # Sample 1 sequence if available
      if (nrow(sequences) >= 1) {
        sampled_sequence <- sequences %>% 
          slice_sample(n = 1)
        
        # Bind the sampled sequence to the sampled_sequences dataframe
        sampled_sequences <- rbind(sampled_sequences, sampled_sequence)
        metadata <- metadata[setdiff(1:nrow(metadata), which(metadata$Sequence_name == sampled_sequence$Sequence_name)),]
        sequences_sampled_this_year <- sequences_sampled_this_year + 1
        
        }
      # If no sequence is available, do nothing and move to the next country
    }
  }
  
#  if (nrow(sampled_sequences) < round(0.8 * total_required)) {
#    additional_needed <- round(0.8 * total_required) - nrow(sampled_sequences)
#    additional_sequences <- metadata %>%
#      filter(Country %in% countries_with_data) %>%
#      slice_sample(n = additional_needed, replace = FALSE)
#    sampled_sequences <- rbind(sampled_sequences, additional_sequences)
#  }
    
    if (nrow(sampled_sequences) >= total_required) {
      break
    }
  }
  
  return(sampled_sequences)
}


sampled_sequences_list <- list()

for (year in sort(unique(non_selected_country_metadata$Year))) {
  # Shuffle countries to give each one an equal chance in the sampling process
  shuffled_countries <- sample(in_region_countries)
  sampled_sequences <- sample_sequences(non_selected_country_metadata, shuffled_countries, max(non_selected_country_metadata$Year), min(non_selected_country_metadata$Year),in_region_sequences)
  sampled_sequences_list[[length(sampled_sequences_list) + 1]] <- sampled_sequences
}

sampled_sequences_within_region <- bind_rows(sampled_sequences_list)

##############################################################################
# Step 4: Sampling Sequences per Country within a Region
##############################################################################

out_of_region_countries <- unique(non_selected_country_metadata$Country[non_selected_country_metadata$Continent != first_country_region])

# Initialize variables
finished_sampling <- FALSE

# Initialize a list to store the sampled sequences
sampled_sequences_list <- list()

for (year in sort(unique(non_selected_country_metadata$Year))) {
  # Shuffle countries to give each one an equal chance in the sampling process
  shuffled_countries <- sample(out_of_region_countries)
  sampled_sequences <- sample_sequences(non_selected_country_metadata, shuffled_countries, max(non_selected_country_metadata$Year), min(non_selected_country_metadata$Year),out_region_sequences)
  sampled_sequences_list[[length(sampled_sequences_list) + 1]] <- sampled_sequences
}

sampled_sequences_out_of_region <- bind_rows(sampled_sequences_list)

##############################################################################
# Step 5: Combine all sampled sequences into a single data frame
##############################################################################
final_sample <- rbind(metadata_countries_of_interest,sampled_sequences_within_region,sampled_sequences_out_of_region)
##############################################################################
# Step 6: Match sample to sequences and produce sampled metadata + sequences
##############################################################################

row_ID <- as.numeric(final_sample$RowID)

keep <- metadata.df %>%
  filter(RowID %in% row_ID)



serotype <- unique(metadata.df$Serotype)

colnames(final_sample) <- c('GenBank_ID', "Country", "State",
                          "City", "Serotype", "Sequence_Type", "date","Decimal_Date","name","Year","Continent")

final_sample <- select(final_sample,c(9,7,1,2,3,4,5,6,8,10,11))

write.table(final_sample,
            file = paste0("results/subsampled_",serotype,'_infoTbl.txt'),
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

write.csv(final_sample,
          file = paste0("results/subsampled_",serotype,'_infoTbl.csv'),
          row.names = FALSE)

check <- read.csv(paste0("results/subsampled_",serotype,'_infoTbl.txt'),sep ="")
seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

not_in_seq_name <- anti_join(metadata.df, final_sample, by = c("GenBank_ID" = "GenBank_ID"))

keep <- metadata.df %>%
  filter(GenBank_ID %in% final_sample$GenBank_ID)

taxa_split <- data.frame(do.call('rbind',strsplit(as.character(keep$V1),'|',fixed = TRUE)))

species.to.keep <- not_in_seq_name$name

vec.names<-unlist(lapply(strsplit(names(seqs), ";"), function(x)x[length(x)]))

vec.tokeep <-which(vec.names %in%  species.to.keep)

# Write the sequences to a new FASTA file

write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=paste0('results/subsampled_', serotype, '.fasta'))

# Prepare and write the metadata information table
seq_name_kept <- as.data.frame(as.matrix(attributes(seqs[vec.tokeep])$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name_kept$V1),'|',fixed = TRUE)))
taxa_split$name <- seq_name_kept$V1
colnames(taxa_split) <- c('GenBank_ID', "Country", "State",
                          "City", "Serotype", "Sequence_Type", "date","Decimal_Date","name")

taxa_split <- select(taxa_split,c(9,7,1,2,3,4,5,6,8))

write.table(taxa_split,
            file = paste0("results/subsampled_",serotype,'_infoTbl.txt'),
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

write.csv(taxa_split,
          file = paste0("results/subsampled_",serotype,'_infoTbl.csv'),
          row.names = FALSE)
