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

datasets summary virus genome taxon "Dengue Virus" --released-after 2000-01-01 --as-json-lines | dataformat tsv virus-genome > metadata.tsv

unzip /Users/rhysinward/Documents/Dengue_anaysis/virus.zip
```

- Here we downloaded both the fasta  and metadata for all sequences

## Step 2: Process and clean metadata and output fasta

- Load the sequences and metadata and using the R script code/Clean_metadata_and_fasta.R to clean the metadata and output fasta sequences with correct naming

```
Rscript /Users/rhysinward/Documents/Dengue_anaysis/Code/Clean_metadata_and_fasta.R --metadata /Users/rhysinward/Documents/Dengue_anaysis/metadata.tsv --fasta /Users/rhysinward/Documents/Dengue_anaysis/virus/ncbi_dataset/data/genomic.fna
```

- It also splits serotypes into Dengue 1,2,3,4 and Unknown
- This code can also be used to select the required dates of the data
- Need to add arguments to set date within the terminal wasn't currently working and not sure why so date in set in the rscript for now

## Step 3: Checking serotypes and Genotypes 

- Here we would ideally have a command line tool to double check and designte serotypes and genotypes independantly to GenBank due to errors in naming and number of Unknown Sequences
- At the moment to confirm serotype we are using using [Dengue Typing Tool](http://krisp.ukzn.ac.za/app/typingtool/dengue/)
  - Note that the typing tool as a maximum of 100,000 sequences it can run for free (should be enough for all runs pretty much)
  - Code to process the outputs of Genomic Detective and produce the final datasets can be found at xxxxxx
 
## Step 4: Alignment

- Here we are aligning our sequences using nextalign which is wrapped in the nextstrain pipeline
- Please see here about installing [nextstrain](https://docs.nextstrain.org/en/latest/install.html). N.b we are using nextstrain within docker
- To run this please use the bash script found here bash/align_sequences.sh

```
bash bash/align_sequences.sh
```

## Step 5: Seperating E gene and Whole Genomes (WG)

- Here we are splitting the aligned genomes into E gene (EG) segments and Whole genomes (WG)
- We say that if more than 30% of a WG bases (ATCG) are missing then we remove (70% set due to Yale sequencing thresholds)
- We say that is more the 5% of an EG bases are missing then we remove (stricter criteria is arbitrary can be changed)
- The EG position for each serotype is based on the genemap
- To run this please use the rscript found here Code/Seperate_EG_and_WG.R (N.b code is slightly cluncky but works)

```
rscript Code/Seperate_EG_and_WG.R
```

## Step 6: Sub-sampler (WORK IN PROGRESS INTERIM SUB-SAMPLER INCLUDED)


![Sampling Pipeline (1)_page-0001](https://github.com/rhysinward/dengue_pipeline/assets/67955642/06b5a02e-9b14-4324-91e1-79b54d8a1682)
**Figure 1.** Subsampler pipeline.


## Step 7: Tree Building - classic iqtree command



## Step 8: Time-Scaling - Problem keeps infering into the future

## Step 9: Mugration Anaysis 

## Step 10: Visulisation within Nextstrain 

## Step 11: Transmission Lineages




