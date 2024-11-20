# List of required packages
required_packages <- c("optparse", "dplyr", "lubridate", "tidyr",
                       "readr", "ape","seqinr", "stringr")

options(repos = c(CRAN = "http://cran.us.r-project.org"))

# Function to check and install packages
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
    make_option(c("-x", "--extra_metadata"), type="character", help="Option to add a csv file containing metadata from sequencing but not on Genbank"),
    make_option(c("-z", "--extra_fasta"), type="character", help="Option to add a fasta file containing metadata from sequencing but not on Genbank"),
    make_option(c("-j", "--outfile_tsv"), type="character", help="Outfile tsv"),
    make_option(c("-k", "--outfile_fasta"), type="character", help="Outfile fasta"),
    make_option(c("-l", "--outfile_csv"), type="character", help="Outfile csv"),
    make_option(c("-s", "--start_date"), type="character", default="2000-01-01", help="Start date for filtering (format YYYY-MM-DD)"),
    make_option(c("-e", "--end_date"), type="character", default="2024-05-31", help="End date for filtering (format YYYY-MM-DD)"),
    make_option(c("-H", "--host"), type="character", default="Homo sapiens", help="Host type sample for sequencing was taken from")
  )
)

opt = parse_args(opt_parser)

########################################################################
## main
########################################################################
## read in input metadata file from GenBank
if (!is.null(opt$metadata)) {
  metadata.df <- read_tsv(opt$metadata)
} else {
  cat("Input metadata file. Exiting now...")
  quit()
}

## read in and combine extra metadata from sequencing
## Nb will need to be in a specific format to match GenBank metadata
## See Github for specific formatting requirements

if (!is.null(opt$extra_metadata)) {
  metadata.extra <- read_tsv(opt$extra_metadata, show_col_types = FALSE)
  if (nrow(metadata.extra) == 0 || ncol(metadata.extra) == 0) {
    warning("Extra metadata NOT included")
  } else {
    metadata.df <- rbind(metadata.df, metadata.extra)
    print("Extra metadata included")
  }
}

#Process dates

process_date <- function(x){
  date <- ifelse(nchar(x$`Collection date`) == 4, paste(x$`Collection date`, "06-15", sep = "-"),
                 ifelse(nchar(x$`Collection date`) == 7, paste(x$`Collection date`, 15, sep = "-"), x$`Collection date`))
  x$Date <- as.Date(parse_date_time(date, orders = c('mdy','dmy','myd','y','my','m','ymd','ym')))
  
  x <- dplyr :: select(x,-c(`Collection date`))
}

metadata.df <- process_date(metadata.df)

#filter host and date

metadata.df <- metadata.df %>%
  filter(Host == opt[["host"]]) %>%
  filter(Date >= as.Date(opt[["start_date"]]) & Date <= as.Date(opt[["end_date"]]))

#extract state level information

metadata.df <- metadata.df %>%
  separate_wider_delim(
    cols = 'Location', 
    delim = "/", 
    names = c("Continent", "Country", "State", "City"), 
    too_many = "error", 
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    cols = City, 
    delim = ",", 
    names = c("CityPrimary", "CitySecondary"), 
    too_many = "merge", 
    too_few = "align_start"
  ) %>%
  mutate(
    City = coalesce(CityPrimary, CitySecondary)
  ) %>%
  select(-CityPrimary, -CitySecondary) %>%
  mutate(across(
    .cols = c(Continent, Country, State, City), 
    .fns = ~ str_trim(.x, side = "both")
  ))


#select desired columns and remove any with NA in Date and Country

metadata.df <- metadata.df %>%
  select(`Accession ID`, Serotype, Date, Country, State, City) %>% 
  filter(!is.na(Date) & !is.na(Country)) 

#functions

calcDecimalDate_fromTxt	<- function( dateTxt, sep="/", namedMonths=FALSE, dayFirst=FALSE) {
  els 	<- strsplit(dateTxt, sep)[[1]]
  if (dayFirst) {
    if (length(els) > 1) {
      els <- els[length(els):1]
    }
  }
  
  year 	<- as.integer(els[1])
  
  if (length(els)==1) {
    month <- 6  #7
    day	<- 15 #2
    decDate <- year + 0.5
  } else {
    
    if (length(els)==2) {
      if (nchar(els[2]) > 0) {
        if (namedMonths) {
          month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
        } else {
          month <- as.integer(els[2])
        }
        day	<- 15
        decDate <- calcDecimalDate(day,month,year)
      } else {
        month <- 6 #7
        day   <- 15 #2
        decDate <- year + 0.5
      }
    } else {
      if (namedMonths) {
        month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
      } else {
        month <- as.integer(els[2])
      }
      
      if (nchar(els[3]) > 0) {
        day 	<- as.integer(els[3])
      } else {
        day <- 15
      }
      decDate <- calcDecimalDate(day,month,year)
    }
  }
  
  
  return ( decDate )
}

calcDecimalDate_from_yymmdd	<- function( dateTxt, sep="/", ycutoff=15, defaultMonth=6, defaultDay=15 ) {
  els	<- strsplit(dateTxt, sep)[[1]]
  yy	<- as.integer(els[1])
  mm	<- as.integer(els[2])
  dd	<- as.integer(els[3])
  
  if (!is.finite(yy)) {
    return( -1 )
  } else {
    if (yy <= ycutoff) {
      yy <- yy+2000
    }
    if ((yy > ycutoff) & (yy < 99)) {
      yy <- yy+1900
    }
    
    if (!is.finite(mm)) {
      mm <- 0
    }
    if (!is.finite(dd)) {
      dd <- 0
    }
    return ( calcDecimalDate( dd, mm, yy, defaultMonth=defaultMonth, defaultDay=defaultDay ) )
  }
  
}

getEl	<- function(line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
  els	<- strsplit(line, sep)[[1]]
  
  if (ind[1] != -1) {
    if (fromEnd) {
      ind <- length(els)-(ind-1)
    }
  }
  
  if (final) {
    return( els[length(els)] )
  } else {
    
    if (reconstruct) {
      if (ex[1] > 0) {
        if (fromEnd) {
          ex <- length(els)-(ex-1)
        }
        ind <- setdiff((1:length(els)),ex)
      }
      
      newLine <- els[ind[1]]
      if (length(ind) > 1) {
        for (i in 2:length(ind)) {
          newLine <- paste(newLine, els[ind[i]], sep=sep)
        }
      }
      return ( newLine )
    } else {
      if ( ind[1] == -1 ) {
        return( els )
      } else {
        return( els[ind] )
      }
    }
  }
}

calcDecimalDate	<- function(day, month, year, defaultMonth=6, defaultDay=15) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334)
  
  if (month==0) {
    if (defaultMonth >= 1) {
      month <- defaultMonth
    } else {
      month	<- ceiling(runif(1)*12)
    }
  }
  
  if (day==0) {
    if (defaultDay >= 1) {
      day	<- defaultDay
    } else {
      day	<- ceiling(runif(1)*30)
    }
  }
  
  dd	<- cd[month] + day - 1
  
  decDate <- year + (dd/365)
  
  return ( decDate )
}

invertDecimalDate <- function( decDate, formatAsTxt=FALSE, ddmmyy=FALSE ) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
  fractD<- cd/365
  
  year		<- floor(decDate)
  fractYear 	<- decDate-year
  month		<- which(fractD >= fractYear)[1]-1
  
  if (month > 0) {
    fractMonth  <- fractYear-fractD[month]
    day		<- round((fractMonth*365)+1)
  } else {
    month <- 1
    day   <- 1
  }
  
  res <- c(day,month,year)
  
  if (formatAsTxt) {
    if (month < 10) {
      mm  <- paste("0",month,sep="")
    } else {
      mm <- month
    }
    res <- paste(year,mm,day,sep="-")
  }
  
  if (ddmmyy) {
    months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    if (day < 10) {
      dd <- paste("0",day,sep="")
    } else {
      dd <- day
    }
    res <- paste(dd,months[month],year,sep="-")
  }
  return( res )
  
}

#match metadata to fasta file and rename sequences

## read in input fasta file

if (!is.null(opt$fasta)) {
  seqs <- read.fasta(opt$fasta)
} else {
  cat("Input fasta file. Exiting now...")
  quit()
}

if (!is.null(opt$extra_fasta)) {
  safe_read_fasta <- purrr::safely(\(x) read.fasta(x))
  seqs_extra <- safe_read_fasta(opt$extra_fasta)
  
  if (is.null(seqs_extra$error)) {
    seqs <- c(seqs, seqs_extra$result)
    print("Extra fasta data included")
  } else {
    warning("Extra fasta data NOT included")
  }
} 

taxa <- as.matrix(attributes(seqs)$names)
genbank_ID <- apply(taxa, 1, getEl, ind=2, sep="\\|")
minds  <- match(genbank_ID, metadata.df$`Accession ID`)
dateTxt <- as.matrix(as.character(metadata.df$Date[minds]))
Virus <- as.matrix(as.character(metadata.df$Serotype[minds]))
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

write.fasta(sequences=seqs[vec.tokeep], names=names(seqs)[vec.tokeep],
            file.out=opt$outfile_fasta)

# Prepare and write the metadata information table

seq_name_kept <- as.data.frame(as.matrix(attributes(seqs[vec.tokeep])$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(seq_name_kept$V1),'|',fixed = TRUE)))
taxa_split$name <- seq_name_kept$V1
colnames(taxa_split) <- c('GenBank_ID', "Country", "State",
                          "City", "Virus_name", "Date","Decimal_Date","Sequence_name")

write.table(taxa_split,
            file = opt$outfile_tsv,
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

write.csv(taxa_split,
          file = opt$outfile_csv,
          row.names = FALSE)
cat(paste0("Processing metadata and sequences completed","\n"))
