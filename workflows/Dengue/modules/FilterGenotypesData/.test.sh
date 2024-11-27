#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/metadata" ]; then
    ln -s ../../AddGenotypeInformation/results/out results/metadata
fi
if [ ! -d "results/metadata" ]; then
    echo "No input data found. Please run the AddGenotypeInformation pipeline first."
    exit 1
fi
if [ ! -d "results/fasta" ]; then
    ln -s ../../SplitGenomeAndQC/results/out results/fasta
fi
if [ ! -d "results/fasta" ]; then
    echo "No input data found. Please run the SplitGenomeAndQC pipeline first."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/dengue_target.yaml target
