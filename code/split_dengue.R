# List of required packages
required_packages <- c("optparse", "dplyr", "lubridate", "tidyr", "readr", "ape", "seqinr", "readr")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type="character", help="Input csv file containing cleaned metadata, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type="character", help="Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type="character", default="subsampled", help="Base name for output files. Files will be named as '<outfile>_'serotype'.fasta', '<outfile>__'serotype'_infoTbl.tsv', and '<outfile>__'serotype'_infoTbl.csv'")
    
  )
)
opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

## read in input cleaned metadata file
metadata_df <- safe_read_file_param(opt$metadata, read_csv, show_col_types = FALSE, required = TRUE)

#merge different serotype naming schemes
metadata_df <- metadata_df %>%
  mutate(serotype = case_when(
    grepl("dengue_virus_2|Dengue_virus_2|dengue_virus_type_2", Virus_name, ignore.case = TRUE) ~ "Dengue_2",
    grepl("dengue_virus_3|Dengue_virus_3|dengue_virus_type_3", Virus_name, ignore.case = TRUE) ~ "Dengue_3",
    grepl("dengue_virus_4|Dengue_virus_4|dengue_virus_type_4", Virus_name, ignore.case = TRUE) ~ "Dengue_4",
    grepl("dengue_virus_1|Dengue_virus_1|dengue_virus_type_1|dengue_virus_type_I|Dengue_virus", Virus_name, ignore.case = TRUE) ~ "Dengue_1"
  ))


info_msg("Serotype table:")
table(metadata_df$serotype)

## read in input fasta file
seqs <- safe_read_file_param(opt$fasta, read.fasta, required = TRUE)

taxa <- as.matrix(attributes(seqs)$names)

#loop through each serotype to match, rename, and write to separate files
for (serotype in c("Dengue_1", "Dengue_2", "Dengue_3", "Dengue_4")) {
  # Subset the metadata for the current serotype
  serotype_metadata <- metadata.df %>% filter(serotype == !!serotype)
  minds  <- match(taxa, serotype_metadata$Sequence_name)
  genbank_ID <- as.matrix(serotype_metadata$GenBank_ID[minds])
  dateTxt <- as.matrix(as.character(serotype_metadata$Date[minds]))
  Virus <- serotype
  decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  country <- as.matrix(serotype_metadata$Country[minds])
  state   <- as.matrix(serotype_metadata$State[minds])
  city   <- as.matrix(serotype_metadata$City[minds])
  
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
  
  write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
              file.out=paste0(opt$outfile,serotype,".fasta"))

  # Prepare and write the metadata information table
  
  seq_name_kept <- as.data.frame(as.matrix(attributes(seqs[vec.tokeep])$names))
  taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name_kept$V1),'|',fixed = TRUE)))
  taxa_split$name <- seq_name_kept$V1
  colnames(taxa_split) <- c('GenBank_ID', "Country", "State",
                            "City", "Serotype", "Date","Decimal_Date","Sequence_name")
  
  write.table(taxa_split,
              file = paste0(opt$outfile,serotype,"_infoTbl.txt"),
              sep = "\t", 
              col.names = TRUE, 
              row.names = FALSE, 
              quote = FALSE)
  
  write.csv(taxa_split,
            file = paste0(opt$outfile,serotype,"_infoTbl.csv"),
            row.names = FALSE)
  cat(paste0("Processing completed for ", serotype, "\n"))
}
