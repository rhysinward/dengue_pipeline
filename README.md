# Dengue pipeline 

This repository hosts the comprehensive suite of code and datasets utilized for the routine processing and analysis of Dengue virus data sourced from GenBank.

# Pipeline Workflow 

## Step 1: Acquisition of Genomic Data and Metadata from GenBank

Process:
- Acquire data corresponding to Dengue virus serotypes directly from GenBank.

Commands:
First, set up the required environment and tools:

```
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
```

Next, download the Dengue virus dataset and metadata:

```
conda activate ncbi_datasets

# Download the virus dataset
datasets download virus genome taxon "Dengue Virus"  --filename virus.zip

# Download and format metadata
datasets summary virus genome taxon "Dengue Virus" --released-after 2000-01-01 --as-json-lines | dataformat tsv virus-genome > metadata.tsv

# Unzip the dataset
unzip /Users/rhysinward/Documents/Dengue_anaysis/virus.zip
```

Highlights:

- This step involves downloading both FASTA format sequences and associated metadata for all Dengue virus sequences.

## Step 2: Processing and Cleaning Metadata and Output FASTA

Process:

- Utilize the Clean_metadata_and_fasta.R script to refine the metadata and generate appropriately named FASTA sequence files.

Commands:

- Run the R script with necessary arguments:

```
Rscript /Users/rhysinward/Documents/Dengue_anaysis/Code/Clean_metadata_and_fasta.R --metadata /Users/rhysinward/Documents/Dengue_anaysis/metadata.tsv --fasta /Users/rhysinward/Documents/Dengue_anaysis/virus/ncbi_dataset/data/genomic.fna --start-date 2000-01-01 --end-date 2023-12-01
```

Highlights:
- The script separates sequences into serotypes Dengue 1, 2, 3, and 4.
- It allows for the selection of data based on specified date ranges.
- Includes functionality to integrate additional metadata and sequences not currently available on GenBank.
- Harmonizes varying formats for dates and other metadata fields to maintain consistency across the dataset.

## Step 3: Verifying Serotypes and Genotypes

Objective:

- Implement a method for independently verifying serotypes and genotypes, as GenBank entries may contain inconsistencies or unknown sequences.

Current Approach:

- For serotype confirmation, we utilize the [Dengue Typing Tool](http://krisp.ukzn.ac.za/app/typingtool/dengue/)
- Note: This tool is capable of processing up to 100,000 sequences at no cost, which should suffice for most datasets.

Ideal Approch: 
- Develop a robust command-line tool to independently verify and assign serotypes and genotypes and be easily be integrated within Dengue Pipelines. This tool aims to address and rectify potential errors in naming and the presence of unknown sequences often encountered in GenBank entries.
 
## Step 4: Sequence Alignment

Process:

- We employ nextalign, part of the Nextstrain suite, for aligning our sequences.

Setup Instructions:

- Detailed installation guidelines for Nextstrain are available [here](https://docs.nextstrain.org/en/latest/install.html). We recommend using Nextstrain within Docker for optimal performance.
- To run this please use the bash script found here bash/align_sequences.sh

Execution:

- To align sequences, use the following bash script: bash/align_sequences.sh.

```
bash bash/align_sequences.sh
```

## Step 5: Segregating E Gene and Whole Genomes (WG) with Custom Thresholds

Purpose:

- This step involves dividing aligned genomes into Envelope gene (EG) segments and Whole Genomes (WG).

Criteria (default):

- Whole Genomes (WG) with more than 30% missing bases are excluded. This threshold is set considering Yale's sequencing criteria (70% completeness).
- For the Envelope gene (EG) segments, a stricter criterion is applied: sequences with more than 5% missing bases are removed. This percentage is adjustable based on project requirements.

Implementation Details:

- The EG positions are determined according to each serotype's genemap.
- To exucute this process please use the Rscript - Seperate_EG_and_WG.R

```
rscript Code/Seperate_EG_and_WG.R
```

## Step 6: Sub-sampler (WORK IN PROGRESS INTERIM SUB-SAMPLER INCLUDED)


![Sampling Pipeline (1)_page-0001](https://github.com/rhysinward/dengue_pipeline/assets/67955642/06b5a02e-9b14-4324-91e1-79b54d8a1682)
**Figure 1.** Subsampler pipeline.


## Step 7: Tree Building

```
nohup iqtree2 -m TIM2+F+R4 -s x.fasta
```

## Step 8: Time-Scaling

```
treetime --tree x.treefile --aln x.fasta --dates x.csv  --clock-filter 4 --confidence
```

## Step 9: Mugration Anaysis 

```
treetime mugration --tree x.treefile --states x.csv --attribute country
```

## Step 10: Visulisation within Nextstrain - to come

## Step 11: Transmission Lineages - to come




