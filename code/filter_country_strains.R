## load packages
required_packages <- c(
  "dplyr", "tidyverse", "tidyr", "dplyr",
  "seqinr", "ape", "optparse", "readr"
)

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

opt_parser <- OptionParser(
    option_list = list(
        make_option(c("-c", "--csv"), type = "character", help = "Input strain `infoTbl` csv file."),
        make_option(c("-u", "--country"), type = "character", help = "Country of interest to filter strains."),
        make_option(c("-o", "--outfile"), type = "character", help = "Base name for output files.")
    )
)

opt <- parse_args(opt_parser)

##########################################################
## main
##########################################################

# Load infoTbl csv
info_tbl <- safe_read_file_param(opt$csv, read_csv, required = TRUE)

# filter VN strains
vn_strains <- info_tbl %>%
  filter(Country == opt$country) %>%
  pull(name)

# save filtered results
write_lines(vn_strains, opt$outfile)

info_msg(sprintf("Finished filtering strains from `%s`", opt$country))