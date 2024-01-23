#!/bin/bash

R_SCRIPT="code/Clean_metadata_and_fasta.R"

# Metadata file
METADATA_FILE="data/metadata.tsv"

# Fasta file
FASTA_FILE="data/ncbi_dataset/data/genomic.fna"

# Start date
START_DATE="2000-01-01"

# End date
END_DATE="2023-12-24"

# Run the R script with the specified parameters
Rscript $R_SCRIPT --metadata $METADATA_FILE --fasta $FASTA_FILE --start-date $START_DATE --end-date $END_DATE
