## load packages
required_packages <- c("dplyr","tidyverse","tidyr","dplyr",
                        "seqinr","ape", "optparse")

suppressMessages(
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
    }
    library(package, character.only = TRUE)
  }
)

option_list = list(
  make_option(c("-f", "--fasta"), type="character", help="Input fasta file."),
  make_option(c("-t", "--tree"), type="character", help="Input tree file."),
  make_option(c("-m", "--metadata_china"), type="character", help="Input metadata file."),
  make_option(c("-v", "--metadata_vietnam"), type="character", help="Input metadata file."),
  make_option(c("-o", "--outfile"), type="character", default="subsampled", help="Base name for output files."),
  make_option(c("-l", "--DTA"), type="character", default="location", help="Name of the column to treat as location.")
)

opt = parse_args(OptionParser(option_list=option_list))

##########################################################
## main
##########################################################

# Load FASTA

if (!is.null(opt$fasta)) {
  fasta <- read.fasta(opt$fasta)
} else {
  cat("Input fasta file. Exiting now...")
  quit()
}

seq_name <- data.frame(name = sapply(fasta, function(x) attr(x, "name")))

# Load Tree

if (!is.null(opt$tree)) {
  tree <- read.nexus(opt$tree)
} else {
  cat("Input tree file. Exiting now...")
  quit()
}

tree_name <- data.frame(name = tree$tip.label)

#load granular Vietnam locations

if (!is.null(opt$metadata_vietnam)) {
  denv_vietnam <- read.csv(opt$metadata_vietnam)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

process_date <- function(x){
  date <- ifelse(nchar(x$'Collection_Date') == 4, paste(x$'Collection_Date', "06-15", sep = "-"),
                 ifelse(nchar(x$'Collection_Date') == 7, paste(x$'Collection_Date', 15, sep = "-"), x$'Collection_Date'))
  x$Date <- as.Date(parse_date_time(date, orders = c('mdy','dmy','myd','y','my','m','ymd','ym')))
  
  x <- dplyr :: select(x,-c('Collection_Date'))
}

metadata.df_vietnam <- process_date(denv_vietnam)

#select desired columns and remove any with NA in Date and Country

metadata.df_vietnam <- metadata.df_vietnam %>%
  dplyr:: select(Accession,Organism_Name, Region, Date) 

#merge different serotype naming schemes
metadata.df_vietnam <- metadata.df_vietnam %>%
  mutate(serotype = case_when(
    grepl("dengue virus 2|Dengue virus 2|dengue virus type 2", Organism_Name, ignore.case = TRUE) ~ "Dengue_2",
    grepl("dengue virus 3|Dengue virus 3|dengue virus type 3", Organism_Name, ignore.case = TRUE) ~ "Dengue_3",
    grepl("dengue virus 4|Dengue virus 4|dengue virus type 4", Organism_Name, ignore.case = TRUE) ~ "Dengue_4",
    grepl("dengue virus 1|Dengue virus 1|dengue virus type 1|dengue virus type I|Dengue virus", Organism_Name, ignore.case = TRUE) ~ "Dengue_1"))

metadata.df_complete_vietnam <- metadata.df_vietnam %>%
  filter(!Region %in% c("Imported case","Unknown"))

metadata.df_complete_vietnam <- metadata.df_complete_vietnam %>%  mutate(Region = case_when(
  grepl("Central|Central (maybe north)", Region, ignore.case = TRUE) ~ "Central_Vietnam",
  grepl("North (maybe central)|North", Region, ignore.case = TRUE) ~ "North_Vietnam",
  grepl("South", Region, ignore.case = TRUE) ~ "Southern_Vietnam"))

#load granular china locations

if (!is.null(opt$metadata_china)) {
  denv_china <- read.csv(opt$metadata_china)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

metadata.df_china <- process_date(denv_china)

#select desired columns and remove any with NA in Date and Country

metadata.df_china <- metadata.df_china %>%
  dplyr :: select(Accession,Organism_Name, Region, Date) 

#merge different serotype naming schemes
metadata.df_china <- metadata.df_china %>%
  mutate(serotype = case_when(
    grepl("dengue virus 2|Dengue virus 2|dengue virus type 2", Organism_Name, ignore.case = TRUE) ~ "Dengue_2",
    grepl("dengue virus 3|Dengue virus 3|dengue virus type 3", Organism_Name, ignore.case = TRUE) ~ "Dengue_3",
    grepl("dengue virus 4|Dengue virus 4|dengue virus type 4", Organism_Name, ignore.case = TRUE) ~ "Dengue_4",
    grepl("dengue virus 1|Dengue virus 1|dengue virus type 1|dengue virus type I|Dengue virus", Organism_Name, ignore.case = TRUE) ~ "Dengue_1"))

table(metadata.df_china$Region)

metadata.df_complete_china <- metadata.df_china %>%
  filter(!Region %in% c("Imported case","Unknown"))

metadata.df_complete_china$Region <- iconv(metadata.df_complete_china$Region, to = "UTF-8")

# Replace "\xcaHubei" with "Hubei" in the Region column
metadata.df_complete_china$Region <- gsub("\\\\xcaHubei", "Hubei", metadata.df_complete_china$Region)

# Replace "Yunnan\\xca" with "Yunnan" in the Region column
metadata.df_complete_china$Region <- gsub("Yunnan\\\\xca", "Yunnan", metadata.df_complete_china$Region)

#merge the two dataframes

metadata.df_complete <- rbind(metadata.df_complete_vietnam,metadata.df_complete_china)

metadata.df_complete <- dplyr:: select(metadata.df_complete,c("Accession","Region"))

# Filter FASTA based on tree tips
to_keep <- seq_name$name %in% tree_name$name
filtered_fasta <- fasta[to_keep]

# Metadata for DTA

fasta_split <- data.frame(do.call('rbind',strsplit(as.character(names(filtered_fasta)),'|',fixed = TRUE)))

fasta_split <- left_join(fasta_split,metadata.df_complete, by = c("X1" = "Accession"))

#if NA in region, place with value from X2

fasta_split$Region[is.na(fasta_split$Region)] <- fasta_split$X2[is.na(fasta_split$Region)]

#new dataframe for DTA metadata

dta_metadata <- data.frame(traits = names(filtered_fasta), location = fasta_split$Region)

#write table
write.table(dta_metadata, file = paste0(opt$outfile, "_metadata.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
# Write filtered FASTA
write.fasta(sequences = filtered_fasta, names = seq_name$name[to_keep], file.out = paste0(opt$outfile, "_filtered.fasta"))
