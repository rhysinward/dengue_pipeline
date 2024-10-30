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

## Step 1a: Acquisition of Genomic Data and Metadata from GenBank

- This step involves downloading both FASTA format sequences and associated metadata for all Dengue virus sequences.
- You are able to change the date from which you want to obtain sequences from

There is a scipt in the next step that allows for the addition of metadata + sequences outside of GenBank e.g. GISAID

## Step 2a: Clean metadata and FASTA files 

- This step processes and cleans the data data downloaded from NCBI
- It allows for the selection of data based on specified date ranges and host type.
- Includes functionality to integrate additional metadata and sequences not currently available on GenBank.
- Harmonizes varying formats for dates and other metadata fields to maintain consistency across the dataset.

Please see the schema for the metadata and naming of fasta file below if you want to merge datasets taken from multiple sources

## Metadata Schema and FASTA File Naming Convention

To ensure consistency when merging datasets from multiple sources, please adhere to the following metadata schema and FASTA file naming conventions.

## Metadata Fields

| Field          | Description                                         | Format                             | Example                        | Note                         |
|----------------|-----------------------------------------------------|------------------------------------|--------------------------------|-------------------------------|
| GenBank_ID     | Unique identifier assigned by GenBank to each sequence | Alphanumeric string                | PP773768.1                     | Required                      |
| Country        | Country where the sample was collected              | Full country name                  | Thailand                       | Required       |
| State          | State or province where the sample was collected    | Full state/province name           | Chanthaburi_Province                             | Use NA if not available       |
| City           | City where the sample was collected                 | Full city name                     | Chanthaburi                             | Use NA if not available       |
| Serotype       | Serotype of the virus                               | Alphanumeric string without spaces | Dengue_1                       |         Required         |
| Date           | Date when the sample was collected                  | YYYY-MM-DD (ISO 8601 format)       | 2018-11-26                     |  Required    |
| Decimal_Date   | Collection date as a decimal for computational analyses | Decimal number with up to 14 decimal places | 2018.90136986301 |    Required   |
| Sequence_name  | Concatenated string of key metadata fields          | GenBank_ID\|Country\|State\|City\|Serotype\|Date\|Decimal_Date | PP773768.1\|Thailand\|NA\|NA\|Dengue_1\|2018-11-26\|2018.90136986301 | Required |

## FASTA File Naming Convention

Each sequence in the FASTA file should use the `Sequence_name` as the header line.

### Example of a FASTA entry:

```
>PP773768.1|Thailand|Chanthaburi_Province|Chanthaburi|Dengue_1|2018-11-26|2018.90136986301
ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC...
```

## Step 2b (Optional): Processing of Metadata from Outside of Genbank (example of data taken from GISAID - note that this data is not in the repro and if you don't want to do these steps please delete from snakefile)

- You will need your own personal account to download data from [GISAID](https://gisaid.org/).
- This step harmonizes metadata and FASTA files to match the schema format.

## Step 2c (Optional): De-duplication of Sequences and Metadata from GenBank and GISAID

- Since some sequences are present in both sources, it's necessary to remove duplicates before creating a joint dataset.
- Due to low sampling of dengue, identical sequences across datasets are expected to be rare.
- We identify duplicates by comparing sequences with the same country and collection date; if sequences are identical, we de-duplicate them.

## Step 2d (Optional): Concatenate Sequences from GenBank and GISAID (applicable to any metadata + FASTA following the same naming convention)

- Merge FASTA and metadata from both sources.

## Step 3: filter for sequences from SEA

- Filter data for countries in SEA and select only sequences from china and Vietnam with known geo-coded sequences
- This is quite a specific step for our anaysis, can be removed to make the pipeline more generalisable

## Step 4: Split into serotype, add serotypes to sequence name and generate sequence specific metadata

- The script separates sequences into serotypes Dengue 1, 2, 3, and 4.
- Generates serotype specific metadata

## Step 5: Asign Serotypes and Genotypes (Future Step, Not Currently Implemented in Pipeline)

- Assignation of genotype and serotype against reference genomes.

### Future Step:

- Incorporate Verity Hill's [dengue lineage classification](https://pmc.ncbi.nlm.nih.gov/articles/PMC11118645/) for more detailed lineage information.

## Step 6: Sequence alignment 

- Sequences are aligned using nextalign

## Step 7: Segregating E gene and Whole Genomes and performing quaility control

- Segregating E gene and whole genomes from aligned Dengue virus sequences and performing quality control.
- Can set the threshold in which an ambiguous number of bases is acceptable within the data
- Here both Whole Genomes (WG) and E genes with more than 31% missing bases are excluded. This threshold is set considering the [Grubaugh Lab](https://grubaughlab.com/) in Yale's sequencing criteria (69% completeness).

**Note:** All current steps are performed separately for Whole Genomes (WG) and E-genes.

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

# Steps 10-19: Maximum Likelihood + Nextstrain Analysis Workflow

## Step 10: Maximum Likelihood (ML) Treebuilding

- This step constructs a Maximum Likelihood (ML) phylogenetic tree
- Uses IQ-TREE2, which allows users to choose the model best suited for their dataset

## Step 11: Build time-calibrated trees

- Inferring time-calibrated trees for each Dengue virus serotype using treetime

## Steps 12 - 15: Nextstrain Analysis

All of the following steps utilize the Nextstrain suite of tools, which provides a comprehensive pipeline for pathogen genomic analysis. Instructions for installing Nextstrain are available [here](https://docs.nextstrain.org/en/latest/install.html).

### Step 12: Infer "Ancestral" Mutations Across the Tree

- Identifies and maps mutations at each ancestral node within the phylogenetic tree

### Step 13: Translate Sequences

- Converts nucleotide sequences to amino acid sequences

### Step 14: Discrete Trait Reconstruction

- Performs discrete trait reconstruction to infer the geographic or epidemiological traits at each node in the phylogenetic tree

### Step 15: Export for Visualization in Auspice

- Exports the processed and annotated data for visualization in Auspice

## Step 16: Extract Annotated Tree from Nextstrain JSON Format

- Extracts the fully annotated phylogenetic tree from the JSON file produced by Nextstrain

## Step 17: Plot Tree

- Plot tree

## Step 18: Extract Information from Tree

- Retrieves specific annotations such as the number of imports/exports

## Step 19: Quantify number of exports and imports from desired country 

- Plots the number of imports/exports

# Steps 20-24: BEAST Analysis Workflow

## Step 20: Prepare Filtered FASTA and Metadata for BEAST

- Generate updated FASTA files and metadata tailored for BEAST analysis using pruned data from TreeTime.

## Step 21: Generate BEAST XML Configuration Using BEAUti

- Create a BEAST XML configuration file using BEAUti

## Step 22: Run BEAST Analysis

- Perform Bayesian phylogenetic analysis using BEAST.

## 23: Analyze BEAST Output with Tracer

- Assess convergence and performance of BEAST MCMC analysis using Tracer.

## 24: LogCombiner (only use if running multiple runs)

- Merge multiple BEAST runs.

## 24: Generate Maximum Clade Credibility (MCC) Tree with TreeAnnotator

- Summarise BEAST trees into an MCC tree for interpretation.

## 
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




