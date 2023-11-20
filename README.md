# Dengue pipeline 

This respoitory contains all the code and data used in routinely processing and analysing dengue data from GenBank.

# Step by step process of pipeline 

## Step 1: Download genomic and metadata from Genbank

- Download data for serotypes from GenBank

```
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

conda activate ncbi_datasets

datasets download virus genome taxon "Dengue Virus"  --filename virus.zip

datasets summary virus genome taxon "Dengue Virus" --as-json-lines | dataformat tsv virus-genome > metadata.tsv

unzip /Users/rhysinward/Documents/Dengue_anaysis/virus.zip
```
- Here we downloaded both the fasta  and metadata for all sequences

## Step 2: Process and clean metadata and output fasta

- Load the sequences and metadata and using the R script xxxx to clean the metadata and output fasta sequences with correct naming

```
Rscript /Users/rhysinward/Documents/Dengue_anaysis/Code/Clean_metadata_and_fasta.R --metadata /Users/rhysinward/Documents/Dengue_anaysis/metadata.tsv --fasta /Users/rhysinward/Documents/Dengue_anaysis/virus/ncbi_dataset/data/genomic.fna
```
- It also splits serotypes into Dengue 1,2,3,4 and Unknown
- This code can also be used to select the required dates of the data

## Step 3: Checking serotypes and Genotypes 

- Here we would ideally have a command line tool to double check and designte serotypes and genotypes independantly to GenBank due to errors in naming and number of Unknown Sequences
- At the moment to confirm serotype we are using using [Genome Detective](http://krisp.ukzn.ac.za/app/typingtool/dengue/))
  - Note that the typing tool as a maximum of 100,000 sequences it can run for free (should be enough for all runs pretty much)
  - Code to process the outputs of Genomic Detective and produce the final datasets can be found at xxxxxx
- This is the final processing step




