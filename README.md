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

- Whole Genomes (WG) with more than 29% missing bases are excluded. This threshold is set considering the [Grubaugh Lab](https://grubaughlab.com/) in Yale's sequencing criteria (69% completeness).
- For the Envelope gene (EG) segments, a stricter criterion is applied: sequences with more than 5% missing bases are removed.
- Both thresholds are adjustable based on project requirements.

Implementation Details:

- The EG positions are determined according to each serotype's genemap.
- To exucute this process please use the Rscript - Seperate_EG_and_WG.R

```
rscript Code/Seperate_EG_and_WG.R --WG_threshold 0.29 --EG_threshold 0.05
```

## Step 6: Sub-sampler

![Sampling Pipeline (1)_page-0001](https://github.com/rhysinward/dengue_pipeline/assets/67955642/06b5a02e-9b14-4324-91e1-79b54d8a1682)
**Figure 1.** Subsampler pipeline.

To use subsampler please use the following line of code:

```
rscript Code/subsampler.R --metadata results/Dengue_1_infoTbl.csv --fasta results/Dengue_1_EG.fasta --country Colombia  --number_sequences 3000 --prop_rd 0.8 --prop_or 0.2
```

## Step 6b (Optional): Add Rooting Sequences to Sub-sampled Datasets

Objective:

To enhance phylogenetic analysis, this optional step involves appending specific rooting sequences to each sub-sampled Dengue virus serotype dataset.

```
rscript Code/add_rooting_sequence.R
```

## Step 7: Tree Building

```
nohup /home/zoo/zool2561/rhys/dengue_pipeline/iqtree-2.2.6-Linux/bin/iqtree2 -m MFP -s Dengue_1_combined.fasta -o 'EF457905.1_Dengue_virus_type1_isolate_P72-1244_complete_genome'
```

```
nohup /home/zoo/zool2561/rhys/dengue_pipeline/iqtree-2.2.6-Linux/bin/iqtree2 -m MFP -s Dengue_2_combined.fasta -o 'EU003591.1_Dengue_virus_type_2_isolate_IBH11234_polyprotein_gene_complete_cds'
```

```
nohup /home/zoo/zool2561/rhys/dengue_pipeline/iqtree-2.2.6-Linux/bin/iqtree2 -m MFP -s Dengue_3_combined.fasta -o 'KU050695.1_Dengue_virus_3_complete_genome'
```

```
nohup /home/zoo/zool2561/rhys/dengue_pipeline/iqtree-2.2.6-Linux/bin/iqtree2 -m MFP -s Dengue_4_combined.fasta -o 'JF262780.1_Dengue_virus_4_isolate_P73-1120_complete_genome'
```

## Step 8: Time-Scaling

```
treetime \
  --tree Dengue_1_combined.fasta.treefile \
  --aln Dengue_1_combined.fasta \
  --dates Dengue_1_combined_metadata.csv \
  --clock-filter 3 \
  --stochastic-resolve \
  --reroot EF457905.1_Dengue_virus_type1_isolate_P72-1244_complete_genome \
  --outdir dengue_1
```

```
treetime \
  --tree Dengue_2_combined.fasta.treefile \
  --aln Dengue_2_combined.fasta \
  --dates Dengue_2_combined_metadata.csv \
  --clock-filter 3 \
  --stochastic-resolve \
  --reroot EU003591.1_Dengue_virus_type_2_isolate_IBH11234_polyprotein_gene_complete_cds \
  --outdir dengue_2
```

```
treetime \
  --tree Dengue_3_combined.fasta.treefile \
  --aln Dengue_3_combined.fasta \
  --dates Dengue_3_combined_metadata.csv \
  --clock-filter 3 \
  --stochastic-resolve \
  --reroot KU050695.1_Dengue_virus_3_complete_genome \
  --outdir dengue_3
```

```
treetime \
  --tree Dengue_4_combined.fasta.treefile \
  --aln Dengue_4_combined.fasta \
  --dates Dengue_4_combined_metadata.csv \
  --clock-filter 3 \
  --stochastic-resolve \
  --reroot JF262780.1_Dengue_virus_4_isolate_P73-1120,_complete_genome \
  --outdir dengue_4
```

## Step 9: Removal of Outliers Based on Clock Filter

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

## Step 9b: Extract upto date meta data from tree tip name 


## Step 10: Mugration Anaysis 

```
treetime mugration --tree x.treefile --states x.csv --attribute country
```

## Step 11: Visulisation within Nextstrain - to come

## Step 12: Transmission Lineages - to come




