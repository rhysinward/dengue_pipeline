## load packages
required_packages <- c("optparse", "dplyr","lubridate","tidyr",
                       "readr","ape","seqinr","countrycode","ggplot2",
                       "ggplot2","purrr")
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
    make_option(c("-n", "--number_sequences"), type="numeric", help="Number of desired sequences (defult is 1:1 ratio with number of sequences from desired country(ies))"),
    make_option(c("-e", "--prop_rd"), type="numeric",default = 0.8, help="Set proportion desired for region of interest (defult is 0.8)"),
    make_option(c("-w", "--prop_or"), type="numeric",default = 0.2, help="Set proprotion desired for outside region of interest (defult is 0.2)")
  )
  )

opt = parse_args(opt_parser)

##########################################################
## main
##########################################################

##########################################################
# Step 1: read in data
##########################################################
metadata_files <- strsplit(as.character(opt$metadata_files), ",")[[1]]
fasta_files <- strsplit(opt$fasta_files, ",")[[1]]

##########################################################
# Step 2: processing of metadata
##########################################################


process_files <- function(metadata_file, fasta_file) {
  # Read in metadata and FASTA files
  metadata.df <- read.csv(metadata_file)
  seqs <- read.fasta(fasta_file)
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

countries_of_interest <- opt$country

#opt$country

countries_of_interest_first <- strsplit(countries_of_interest, ",")[[1]]

metadata_countries_of_interest <- metadata.df %>%
  filter(Country %in% countries_of_interest) }

##############################################################################
# Step 3: Identify the region of the first country
##############################################################################

first_country <- countries_of_interest_first[1]
first_country_region <- countrycode(first_country, origin = "country.name", destination = "continent")

##############################################################################
# Step 4: Input required number of sequences and calculate Sequences per Year
##############################################################################


#if opt$number_sequences is null (I.e. not imputed then automatically select 1:1 ratio) if not subtract 
#required number of sequences from number of sequences from countries of interest
if (is.null(opt$number_sequences)) {
  n_remainder <- as.numeric(nrow(metadata_countries_of_interest))
} else {
  n_remainder <- as.numeric(opt$number_sequences) - nrow(metadata_countries_of_interest) # Total number of sequences desired
}

if (n_remainder < 0) {
  stop("Error: The number of sequences specified (", opt$number_sequences, 
       ") is less than the number of sequences available in the country of interest (", 
       nrow(metadata_countries_of_interest), ").")
}

total_proportion <- opt[["prop_rd"]] + opt[["prop_or"]]

# Check if the total proportion is not equal to 1
if (total_proportion != 1) {
  error_message <- sprintf("Error: The total proportion of 'prop_RD' and 'prop_OR' should be 1, but it is %.2f", total_proportion)
  stop(error_message)
}

# In-region sequences per year
in_region_sequences <- round(opt[["prop_rd"]] * n_remainder)
# Out-of-region sequences per year
out_region_sequences <- round(opt[["prop_or"]] * n_remainder)

##############################################################################
# Step 3: Sampling Sequences per Country within a Region
##############################################################################
  
in_region_countries <- unique(metadata.df$Country[metadata.df$Continent == first_country_region])

same_region_metadata <- metadata.df %>%
  filter(!Country %in% countries_of_interest & Country %in% in_region_countries)

#sample sequences from countries of interest

# Function to perform subsampler

subsampler <- function(df, max_per_year) {
  if (max_per_year < 1) {
    sampled_data <- data.frame()
} else {
  unique_years <- unique(df$Year)
  unique_countries <- unique(df$Country)
  sampled_data <- data.frame()
  
  for (year in unique_years) {
    year_data <- filter(df, Year == year)
    count <- 0
    
    while(count < max_per_year && nrow(year_data) > 0) {
      unique_countries <- sample(unique_countries)  # Randomize the order of countries so same countries aren't always sampled first
      for (country in unique_countries) {
        country_data <- filter(year_data, Country == country)
        if (nrow(country_data) > 0) {
          sampled_data <- rbind(sampled_data, slice_sample(country_data, n = 1))
          count <- count + 1
          if (count == max_per_year) break
        }
      }
      suppressMessages({
        year_data <- anti_join(year_data, sampled_data)
      })
          }
  }
  if (nrow(sampled_data) < max_per_year * length(unique_years)) {
    remaining <- max_per_year * length(unique_years) - nrow(sampled_data)
    additional_data <- df %>%
      filter(!Sequence_name %in% sampled_data$Sequence_name) %>%
      slice_sample(n = remaining)
    sampled_data <- rbind(sampled_data, additional_data)
  }
  sampled_data
}
}

# Apply advanced sampling to the data
sampling_region_of_interest <- subsampler(same_region_metadata, max_per_year = round(in_region_sequences/length(unique(same_region_metadata$Year))))

##############################################################################
# Step 4: Sampling Sequences per Country within a Region
##############################################################################
in_region_countries <- unique(metadata.df$Country[metadata.df$Continent == first_country_region])

out_of_region_countries <- unique(metadata.df$Country[metadata.df$Continent != first_country_region])

out_of_region_metadata <- metadata.df %>%
  filter(!Country %in% countries_of_interest & Country %in% out_of_region_countries)

# Apply advanced sampling to the data
sampling_out_of_region <- subsampler(out_of_region_metadata, max_per_year = round(out_region_sequences/length(unique(out_of_region_metadata$Year))))

##############################################################################
# Step 5: Combine all sampled sequences into a single data frame
##############################################################################
# Assuming metadata_countries_of_interest, sampling_out_of_region, and sampling_region_of_interest are already defined

final_sample <- metadata_countries_of_interest

#if metadata is empty do not join

if (nrow(sampling_out_of_region) == 0 && nrow(sampling_region_of_interest) == 0) {
} else {
  if (nrow(sampling_out_of_region) > 0) {
    final_sample <- rbind(final_sample, sampling_out_of_region)
  }
  if (nrow(sampling_region_of_interest) > 0) {
    final_sample <- rbind(final_sample, sampling_region_of_interest)
  }
}

##############################################################################
# Step 6: Combine all sampled sequences into a single data frame
##############################################################################

#get count of data 

final_sample_plot <- final_sample %>%
  mutate(Group = if_else(Country == first_country, first_country, as.character(Continent))) %>%
  group_by(Group, Year) %>%
  summarise(count = n()) %>%
  ungroup()

#plot data with good color scheme

ggplot(final_sample_plot, aes(x = Year, y = count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Year", y = "Number of Sequences", title = "Number of Sequences per Year") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "Set1")

#export to pdf

serotype <- unique(final_sample$Serotype)

ggsave(paste0('results/subsampled_', serotype, '_summary_plot.pdf'), width = 10, height = 10, units = "in")


##############################################################################
# Step 6: Match sample to sequences and produce sampled metadata + sequences
##############################################################################

colnames(final_sample) <- c('GenBank_ID', "Country", "State",
                            "City", "Serotype", "date","Decimal_Date","name","Year","Continent")

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
                          "City", "Serotype", "date","Decimal_Date","name")

taxa_split <- select(taxa_split,c(8,6,1,2,3,4,5,7))
write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=paste0('results/subsampled_', serotype, '.fasta'))

write.csv(taxa_split,
          file = paste0("results/subsampled_",serotype,'_infoTbl.tsv'),
          row.names = FALSE, quote=FALSE)

write.csv(taxa_split,
          file = paste0("results/subsampled_",serotype,'_infoTbl.csv'),
          row.names = FALSE,quote=FALSE)
