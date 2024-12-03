# List of required packages
required_packages <- c("optparse", "dplyr", "lubridate", "tidyr", "readr", "ape", "seqinr", "stringr", "magrittr")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type="character", help="Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type="character", help="Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-x", "--extra_metadata"), type="character", help="Option to add a csv file containing metadata from sequencing but not on Genbank"),
    make_option(c("-z", "--extra_fasta"), type="character", help="Option to add a fasta file containing metadata from sequencing but not on Genbank"),
    make_option(c("-j", "--outfile_tsv"), type="character", help="Outfile tsv"),
    make_option(c("-k", "--outfile_fasta"), type="character", help="Outfile fasta"),
    make_option(c("-l", "--outfile_csv"), type="character", help="Outfile csv"),
    make_option(c("-s", "--start_date"), type="character", default="2010-01-01", help="Start date for filtering (format YYYY-MM-DD)"),
    make_option(c("-e", "--end_date"), type="character", default="2023-12-31", help="End date for filtering (format YYYY-MM-DD)"),
    make_option(c("-H", "--host"), type="character", default="Homo sapiens", help="Host Type sample for sequencing was taken from")
  )
)

opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

## read in input metadata file from GenBank
metadata_df <- safe_read_file_param(opt$metadata, read_tsv, show_col_types = FALSE, required = TRUE) %>%
  select(-`Geographic State`)

## read in and combine extra metadata from sequencing
## Nb will need to be in a specific format to match GenBank metadata
## See Github for specific formatting requirements
extra_metadata <- safe_read_file_param(opt$extra_metadata, read_tsv, show_col_types = FALSE)
if (nrow(extra_metadata) == 0 || ncol(extra_metadata) == 0) {
  warn_msg("Extra metadata NOT included")
} else {
  metadata_df <- rbind(metadata_df, extra_metadata)
  info_msg("Extra metadata included")
}

# Process dates
process_date <- function(df) {
  date <- ifelse(nchar(df$"Isolate Collection date") == 4, paste(df$"Isolate Collection date", "06-15", sep = "-"),
    ifelse(nchar(df$"Isolate Collection date") == 7, paste(df$"Isolate Collection date", 15, sep = "-"), df$"Isolate Collection date")
  )
  df$Date <- as.Date(lubridate::parse_date_time(date, orders = c("mdy", "dmy", "myd", "y", "my", "m", "ymd", "ym")))

  df <- dplyr::select(df, -c("Isolate Collection date"))
}

metadata_df <- process_date(metadata_df)

# Filter host and date
metadata_df <- metadata_df %>%
  filter(`Host Name` == opt$host) %>%
  filter(Date >= as.Date(opt$start_date) & Date <= as.Date(opt$end_date))

# Extract state level information
metadata_df <- metadata_df %>%
  separate_wider_delim(
    cols = `Geographic Location`,
    delim = ":",
    names = c("Country", "State"),
    too_many = "error",
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    cols = State,
    delim = ",",
    names = c("State", "City"),
    too_many = "merge",
    too_few = "align_start"
  ) %>%
  mutate(
    City = str_split_i(City, ",", 1)
  ) %>%
  mutate(across(
    .cols = c(Country, State, City),
    .fns = ~ trimws(.x, which = "both")
  ))

# Select desired columns and remove any with NA in Date and Country
metadata_df <- metadata_df %>%
  select(Accession, `Virus Name`, Date, Country, State, City) %>%
  filter(!is.na(Date) & !is.na(Country))

# Match metadata to fasta file and rename sequences
## read in input fasta file
seqs <- safe_read_file_param(opt$fasta, read.fasta, required = TRUE)

extra_seqs <- safe_read_file_param(opt$extra_fasta, read.fasta)
if (is.null(extra_seqs)) {
  warn_msg("Extra fasta sequences are NOT included")
} else {
  seqs <- c(seqs, extra_seqs)
  info_msg("Extra fasta sequences are included")
}

taxa <- as.matrix(attributes(seqs)$names)
genbank_ID <- apply(taxa, 1, getEl, ind=1, sep=" ")

minds  <- match(genbank_ID, metadata.df$Accession)
dateTxt <- as.matrix(as.character(metadata.df$Date[minds]))
Virus <- as.matrix(as.character(metadata.df$`Virus Name`[minds]))
table(Virus)
decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
country <- as.matrix(metadata.df$Country[minds])
state   <- as.matrix(metadata.df$State[minds])
city   <- as.matrix(metadata.df$City[minds])

newTaxa <- paste(genbank_ID,country,state,city,Virus,dateTxt,decDate,sep="|")
newTaxa <- gsub(" ","_",newTaxa)
newTaxa <- gsub("\\(", "", newTaxa) 
newTaxa <- gsub("\\)", "", newTaxa)  
attributes(seqs)$names <- newTaxa

seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))
remove <- filter(seq_name,grepl('NA|NA|NA',V1,fixed = TRUE))
species.to.remove <- remove$V1
vec.names<-unlist(lapply(strsplit(names(seqs), ";"), function(x)x[length(x)]))
vec.tokeep <-which(! vec.names %in%  species.to.remove)

# Write the sequences to a new FASTA file
write.fasta(
  sequences = seqs[vec.tokeep], names = names(seqs)[vec.tokeep],
  file.out = opt$outfile_fasta
)

# Prepare and write the metadata information table
seq_name_kept <- as.data.frame(as.matrix(attributes(seqs[vec.tokeep])$names))
taxa_split <- data.frame(do.call("rbind", strsplit(as.character(seq_name_kept$V1), "|", fixed = TRUE)))
taxa_split$name <- seq_name_kept$V1
colnames(taxa_split) <- c(
  "GenBank_ID", "Country", "State",
  "City", "Virus_name", "Date", "Decimal_Date", "Sequence_name"
)

write.table(taxa_split,
  file = opt$outfile_tsv,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

write.csv(taxa_split,
  file = opt$outfile_csv,
  row.names = FALSE
)

info_msg("Processing metadata and sequences completed")
