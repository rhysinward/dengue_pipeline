# Dengue pipeline 

This repository hosts the comprehensive suite of code and datasets utilized for the routine processing and analysis of Dengue virus data sourced from GenBank.

# Unified Pipeline for Dengue Sequence Analysis Using Snakemake

## Overview

This pipeline is designed for the retrieval, processing, subsampling, and phylogenetic analysis of the latest dengue virus sequences. It leverages Snakemake, a powerful workflow management system, to ensure reproducibility and efficiency in bioinformatics analyses.

## Prerequisites

Before proceeding, ensure you have the following prerequisites installed:

- Git (for cloning this repository)
- Mambaforge (for managing environments and dependencies)

## Installation

### Step 1: Clone the Repository

First, clone this repository to your local machine:

```
git clone git@github.com:rhysinward/dengue_pipeline.git
```

### Step 2: Install Mambaforge

If Mambaforge is not already installed, execute the following commands:

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

### Step 3: Install Snakemake

Install Snakemake using Mamba, which facilitates the installation in an isolated environment:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### Step 4: Activate Snakemake Environment

Activate the Snakemake environment and verify the installation:

```
mamba activate snakemake
snakemake --help
```

## Usage

### Running the Pipeline

Within the cloned repository, execute the following command to run the pipeline:

For Mac:

```
CONDA_SUBDIR=osx-64 snakemake --use-conda --cores 4
```

For Linux:

```
CONDA_SUBDIR=linux-64 snakemake --use-conda --cores 4
```

Adjust the --cores parameter based on your system's capabilities.

### Visualizing Outputs

The pipeline outputs can be visualized locally or online using [Auspice](https://auspice.us/):

```
nextstrain view auspice/
```

# Step by step of the pipeline

# Pipeline Workflow - Summary of each of the steps

## Step 1: Acquisition of Genomic Data and Metadata from GenBank

- This step involves downloading both FASTA format sequences and associated metadata for all Dengue virus sequences.
- You are able to change the date from which you want to obtain sequences from

## Step 2: Clean metadata and FASTA files 

- This step processes and cleans the data data fownloaded from NCBI
- It allows for the selection of data based on specified date ranges and host type.
- Includes functionality to integrate additional metadata and sequences not currently available on GenBank.
- Harmonizes varying formats for dates and other metadata fields to maintain consistency across the dataset.

## Step 3: filter for sequences from SEA

- Filter data for countries in SEA and select only sequences from china and Vietnam with known geo-coded sequences
- This is quite a specific step for our anaysis, can be removed to make the pipeline more generalisable

## Step 4: Split into serotype, add serotypes to sequence name and generate sequence specific metadata

- The script separates sequences into serotypes Dengue 1, 2, 3, and 4.
- Generates serotype specific metadata

## Step 5: (Future step not currently implemented) Verifying Serotypes and Genotypes
 
Objective:

- Implement a method for independently verifying serotypes and genotypes, as GenBank entries may contain inconsistencies or unknown sequences.

Current Approach:

- For serotype confirmation, we utilize the [Dengue Typing Tool](http://krisp.ukzn.ac.za/app/typingtool/dengue/)
- Note: This tool is capable of processing up to 100,000 sequences at no cost, which should suffice for most datasets.

Ideal Approch: 
- Develop a robust command-line tool to independently verify and assign serotypes and genotypes and be easily be integrated within Dengue Pipelines. This tool aims to address and rectify potential errors in naming and the presence of unknown sequences often encountered in GenBank entries.

## Step 6: Sequence alignment 

- Sequences are aligned using nextalign

## Step 7: Segregating E gene and Whole Genomes and performing quaility control

- Segregating E gene and whole genomes from aligned Dengue virus sequences and performing quality control.
- Can set the threshold in which an ambiguous number of bases is acceptable within the data
- Here both Whole Genomes (WG) and E genes with more than 31% missing bases are excluded. This threshold is set considering the [Grubaugh Lab](https://grubaughlab.com/) in Yale's sequencing criteria (69% completeness).

## Step 8: Subsampler

The selection of sequences was performed using a weighted random sampling technique. 

### Even Weighting Scheme

One approch used we coin an 'Even Weighted Sub-Sampling' approach was chosen to ensure a representative and equitable selection of sequences across different locations, collected over various time periods.

Each sample from a geographical location collected in a specifiable time period was assigned a weight according to the formula:

Where:
- `1/x_ij` represents the number of sequences of location `i` collected during time period `j`.

This weighting scheme inversely correlates the weight of a sequence with the abundance of sequences from the same location in a specifiable time period. Therefore, locations with fewer sequences in a given time period are given higher selection probabilities, ensuring a balanced representation of geographies in the final dataset used for analysis.

### Proportional Weighting Scheme 

The 'Proportional Weighted Sub-sampling' method is designed to ensure that the selection of sub-samples is representative of a particular variable of interest, such as the number of cases, mobility rates, etc., across various geographical locations and over different time periods. This approach allows for a more nuanced analysis that takes into account the variable's distribution, ensuring that the sub-sample accurately reflects the broader dataset's characteristics.

In this approach, each sample is associated with a specific geographical location and falls within a definable time period. To ensure representativeness, each sample is assigned a weight based on the following formula:

- `x_ij`  is the value of the variable of interest for location `i` during time period `j`.

The weighting scheme is at the heart of the 'Proportional Weighted Sub-sampling' method, ensuring that each sample's likelihood of selection is directly correlated with the chosen variable of interest within a given time frame. This correlation means that areas or time periods with higher values of the variable of interest will have a proportionally greater influence on the sub-sample. 

# Parameters Description

The script accepts a range of command-line options to customize the input, output, and behavior of the sampling process. Below is a detailed description of each option available in the `option_list`.

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| `-m`, `--metadata` | character | | Input TSV file containing metadata from GenBank, with sequence identifiers matching those in the input FASTA file. |
| `-f`, `--fasta` | character | | Input FASTA file, with sequence identifiers matching those in the metadata file. |
| `-c`, `--location_local` | character | | CSV file specifying granularity and location of interest for the analysis. |
| `-x`, `--location_background` | character | | CSV file to specify the locations for even sampling or to conduct proportional sampling. Please include granularity, location of background, and number of sequences desired (only relevant for proportional sampling). |
| `-t`, `--time_interval` | character | `Year` | Select the sampling interval from Year, Month, or Week. This option determines the temporal granularity of the analysis. |
| `-n`, `--number_sequences_local` | numeric | | Number of desired sequences from the location of interest. |
| `-e`, `--number_sequences_background` | numeric | | Number of desired sequences from the background locations. |
| `-w`, `--sampling_method` | character | `Even` | Choose between even or proportional sampling methods. |
| `-o`, `--outfile` | character | `subsampled` | Base name for output files. Files will be named as `<outfile>_fasta.fasta`, `<outfile>_infoTbl.tsv`, and `<outfile>_infoTbl.csv`. |

## Step 9: Correct metadata and fasta files into the correct format for iqtree and treetime  

- Correct metadata and fasta files into the correct format for iqtree and treetime

## Step 10: ML-Treebuilding

- Uses IQTREE2

## Step 11: Build time-calibrated trees

- Inferring time-calibrated trees for each Dengue virus serotype using treetime

# Step 12 - 15 all utilise the nextstrain suite of tools available [here](https://docs.nextstrain.org/en/latest/install.html)

## Step 12: Infer "ancestral" mutations across the tree

## Step 13: Translate sequences

## Step 14: Discrete trait reconstruction

## Step 15: Export for visualisation in Auspice

## Step 16: Extract annotated tree from nextstrain JSON format 

- Extract annotated tree from the JSON file produced by the nextstrain suite of tools

## Step 17: Extract information from tree 

- Extract annotations from the tree

## 18: Quantify number of exports and imports from desired country



# Miscellaneous 

## Add Rooting Sequences to Sub-sampled Datasets

Objective:

To enhance phylogenetic analysis, this optional step involves appending specific rooting sequences to each sub-sampled Dengue virus serotype dataset.

```
rscript Code/add_rooting_sequence.R
```
## Removal of Outliers Based on Clock Filter

This step involves the removal of outliers identified using the clock filter from TreeTime. You have two options: either completely remove these sequences and re-estimate the tree using IQ-TREE 2, or remove the tips from the tree using gotree. We will include the rooting sequences within the outliars so they are removed before we have the final tree. 

### Option 1: Remove Sequences and Re-estimate Tree Using IQ-TREE 2

* Identify Outliers: Use the clock filter from TreeTime to identify outlier sequences.
* Remove Sequences: Exclude these outlier sequences from your dataset.
* Re-estimate Tree (following step 7)

### Option 2: Remove Tips from Tree Using gotree
* Identify Outliers: Use the clock filter from TreeTime to pinpoint the outlier tips in your tree.
* Remove Tips Using gotree: Use the gotree tool to remove these tips from your ML tree and repeat step 8. This step does not require re-estimating the entire tree.

example command (can't get it to run outside docker atm): 

```
docker run --platform linux/amd64 -v $PWD:$PWD -w $PWD -i -t evolbioinfo/gotree:v0.2.8b prune -f dengue_1/outliers.tsv -i Dengue_1_combined_trimmed.newick --format newick > Dengue_1_combined_trimmed.newick
```




